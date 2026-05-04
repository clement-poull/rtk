#pragma once

#include "scion/scion.hpp"

#include "scion/scion.hpp"

#include "rtk/core/projection.hpp"

#include "rtk/dif/utils.hpp"

#include "rtk/core/vertices.hpp"

#include "rtk/model/pifs.hpp"

#include "rtk/core/math.hpp"

#include "Eigen/Dense"

#include <span>
#include <optional>
#include <expected>
#include <ranges>
#include <cmath>
#include <algorithm>

/// \file rtk/dif/dcf.hpp
/// \brief provides a tool to compute the Surface Differential Characteristic Function of a tensor product of two contractive transformations
/// \doi 10.5220/0013191700003912

namespace rtk::dif {
    template <scion::usize D>
    auto sdcf(rtk::projection const& projection, std::span<scion::f32 const> const evals_s, std::span<rtk::vecxf const> const evecs_s, rtk::vecxf const& starting_point_s, std::span<scion::f32 const> evals_t, std::span<rtk::vecxf const> evecs_t, rtk::vecxf starting_point_t, scion::f32 const s, scion::f32 const t) -> std::expected<rtk::vertex<scion::f32, D, D>, std::monostate> {
        { // check preconditions
            if (evals_s.size() < 2) {
                throw std::logic_error("width_proj_s must be at least 2");
            }

            if (evals_t.size() < 2) {
                throw std::logic_error("width_proj_t must be at least 2");
            }

            if (projection.width() != evals_s.size() * evals_t.size()) {
                throw std::logic_error("inconsistent width_proj");
            }

            if (evals_s.size() != evecs_s.size()) {
                throw std::logic_error("inconsistent width_proj_s");
            }

            if (evals_t.size() != evecs_t.size()) {
                throw std::logic_error("inconsistent width_proj_s");
            }

            if (!rtk::pifs<>::validate_evals(evals_s).has_value()) {
                return std::unexpected(std::monostate{});
            }
            if (!rtk::pifs<>::validate_evecs(evecs_s).has_value()) {
                return std::unexpected(std::monostate{});
            }

            if (!rtk::pifs<>::validate_evals(evals_t).has_value()) {
                return std::unexpected(std::monostate{});
            }
            if (!rtk::pifs<>::validate_evecs(evecs_t).has_value()) {
                return std::unexpected(std::monostate{});
            }
        }

        // computes the coordinates of s_starting_point in the eigenbasis s_eigenvecs
        std::vector<scion::f32> xis_s = rtk::dif::compute_xis<scion::f32>(evals_s.size(), evecs_s, starting_point_s);
        xis_s[0] = 1.0f;
        // computes the coordinates of t_starting_point in the eigenbasis t_eigenvecs
        std::vector<scion::f32> xis_t = rtk::dif::compute_xis<scion::f32>(evals_t.size(), evecs_t, starting_point_t);
        xis_t[0] = 1.0f;

        // computes the s_alphas
        std::vector<scion::f32> ais_s = rtk::dif::compute_ais(evals_s.size(), evals_s);
        ais_s[0] = 0.0f;
        ais_s[1] = 1.0f;
        // computes the t_alphas
        std::vector<scion::f32> ais_t = rtk::dif::compute_ais(evals_t.size(), evals_t);
        ais_t[0] = 0.0f;
        ais_t[1] = 1.0f;

        // computes the position in projective space
        rtk::vecxf const position_bn = std::ranges::fold_left_first(std::views::cartesian_product(std::views::iota(0UZ, evals_s.size()), std::views::iota(0UZ, evals_t.size())) | std::views::transform([&](auto const& indices) -> rtk::vecxf { return xis_s[std::get<0>(indices)] * xis_t[std::get<1>(indices)] * rtk::vector_product_tensor_flatten<scion::f32>(evecs_s[std::get<0>(indices)], evecs_t[std::get<1>(indices)]) * std::pow(s, ais_s[std::get<0>(indices)]) * std::pow(t, ais_t[std::get<1>(indices)]); } ), [](rtk::vecxf const& a, rtk::vecxf const& b) -> rtk::vecxf { return a + b; }).value();

        // computes the position in projected space
        rtk::vecf<D> const position_rn = projection.project<D>(position_bn).value();

        // computes the cross tangents in projective space
        rtk::vecxf const tan_s_bn = (std::ranges::fold_left_first(std::views::cartesian_product(std::views::iota(1UZ, evals_s.size()), std::views::iota(0UZ, evals_t.size())) | std::views::transform([&](auto const& indices) -> rtk::vecxf { return ais_s[std::get<0>(indices)] * xis_s[std::get<0>(indices)] * xis_t[std::get<1>(indices)] * rtk::vector_product_tensor_flatten<scion::f32>(evecs_s[std::get<0>(indices)], evecs_t[std::get<1>(indices)]) * std::pow(s, ais_s[std::get<0>(indices)] - 1.0f) * std::pow(t, ais_t[std::get<1>(indices)]); } ), [](rtk::vecxf const& a, rtk::vecxf const& b) -> rtk::vecxf { return a + b; })).value();
        rtk::vecxf const tan_t_bn = (std::ranges::fold_left_first(std::views::cartesian_product(std::views::iota(0UZ, evals_s.size()), std::views::iota(1UZ, evals_t.size())) | std::views::transform([&](auto const& indices) -> rtk::vecxf { return ais_t[std::get<1>(indices)] * xis_s[std::get<0>(indices)] * xis_t[std::get<1>(indices)] * rtk::vector_product_tensor_flatten<scion::f32>(evecs_s[std::get<0>(indices)], evecs_t[std::get<1>(indices)]) * std::pow(s, ais_s[std::get<0>(indices)]) * std::pow(t, ais_t[std::get<1>(indices)] - 1.0f); } ), [](rtk::vecxf const& a, rtk::vecxf const& b) -> rtk::vecxf { return a + b; })).value();


        // computes the cross tangents in projected space
        rtk::vecf<D> const tan_s_rn = projection.project<D>(tan_s_bn).value();
        rtk::vecf<D> const tan_t_rn = projection.project<D>(tan_t_bn).value();

        // computes the normal in projected space
        rtk::vecf<D> const normal_rn = tan_t_rn.cross(tan_s_rn).normalized();

        return rtk::vertex<scion::f32, D, D>{position_rn, normal_rn};
    }

    template <scion::usize D>
    auto sdcf(rtk::projection const& projection, std::span<scion::f32 const> const evals_s, std::span<rtk::vecxf const> const evecs_s, rtk::vecxf const& starting_point_s, std::span<scion::f32 const> evals_t, std::span<rtk::vecxf const> evecs_t, rtk::vecxf starting_point_t, scion::f32 const sa, scion::f32 const sz, scion::f32 const ta, scion::f32 const tz, scion::usize const count) -> std::expected<std::vector<rtk::vertex<scion::f32, 3, 3>>, std::monostate> {
        std::vector<rtk::vertex<scion::f32, D, D>> sdcf;
        sdcf.reserve(count * count);

        for (scion::usize i_s = 0; i_s < count; ++i_s) {
            scion::f32 const s = sa + (sz - sa) * (static_cast<scion::f32>(i_s) / static_cast<scion::f32>(count - 1));

            for (scion::usize i_t = 0; i_t < count; ++i_t) {
                scion::f32 const t = ta + (tz - ta) * (static_cast<scion::f32>(i_t) / static_cast<scion::f32>(count - 1));

                if (auto result = rtk::dif::sdcf<D>(projection, evals_s, evecs_s, starting_point_s, evals_t, evecs_t, starting_point_t, s, t); result.has_value()) {
                    sdcf.emplace_back(result.value());
                } else {
                    return std::unexpected(std::monostate{});
                }
            }
        }

        return sdcf;
    }

    /// \brief Computes the gaussian curvature of the SDCF as (s, t) approaches (0, 0) if it exists.
    /// \param projection Projection of the model.
    /// \param s_width_proj Number of control points of the s-curve.
    /// \param s_eigenvals Eigenvalues of the s-curve.
    /// \param s_eigenvecs Eigenvectors of the s-curve.
    /// \param s_starting_point Starting point of the s-curve.
    /// \param t_width_proj Number of control points of the t-curve.
    /// \param t_eigenvals Eigenvalues of the t-curve.
    /// \param t_eigenvecs Eigenvectors of the t-curve.
    /// \param t_starting_point Starting point of the t-curve.
    /// \param transformation transformation to apply to the vectors. Use identity if none.
    /// \return gaussian curvature
    // template <scion::usize D, std::floating_point S>
    // auto sdcf_curvature_gaussian(rtk::projection const& projection, scion::usize s_width_proj, std::span<const S> s_eigenvals, std::span<const rtk::vec<S>> s_eigenvecs, rtk::vec<S> s_starting_point, scion::usize t_width_proj, std::span<const S> t_eigenvals, std::span<const rtk::vec<S>> t_eigenvecs, rtk::vec<S> t_starting_point, rtk::mat<S> const& transformation) -> std::expected<S, std::monostate> {
    //     { // check preconditions
    //         // projection is the projection for the tensor product, so its width is the product of the projection with of the two curves
    //         scion::assert_debug(projection.size() == s_width_proj * t_width_proj);
    //
    //         // checks that there is exactly the right number of eigenvalues and eigenvectors
    //         scion::assert_debug(s_width_proj == s_eigenvals.size());
    //         scion::assert_debug(s_width_proj == s_eigenvecs.size());
    //         // checks that the eigenvectors are all the correct size
    //         scion::assert_debug(std::ranges::all_of(s_eigenvecs, [s_width_proj](auto const& eigenvec) { return s_width_proj == eigenvec.size(); }));
    //         // checks that the eigenvectors and eigenvalues don't contain NaN or Inf
    //         scion::assert_debug(std::ranges::all_of(s_eigenvals, [](auto const& eigenval) { return !std::isnan(eigenval) && !std::isinf(eigenval); }));
    //         scion::assert_debug(std::ranges::all_of(s_eigenvecs, [](auto const& eigenvec) { return !eigenvec.hasNaN() && eigenvec.allFinite(); }));
    //
    //         // checks that there is exactly the right number of eigenvalues and eigenvectors
    //         scion::assert_debug(t_width_proj == t_eigenvals.size());
    //         scion::assert_debug(t_width_proj == t_eigenvecs.size());
    //         // checks that the eigenvectors are all the correct size
    //         scion::assert_debug(std::ranges::all_of(t_eigenvecs, [t_width_proj](auto const& eigenvec) { return t_width_proj == eigenvec.size(); }));
    //         // checks that the eigenvectors and eigenvalues don't contain NaN or Inf
    //         scion::assert_debug(std::ranges::all_of(t_eigenvals, [](auto const& eigenval) { return !std::isnan(eigenval) && !std::isinf(eigenval); }));
    //         scion::assert_debug(std::ranges::all_of(t_eigenvecs, [](auto const& eigenvec) { return !eigenvec.hasNaN() && eigenvec.allFinite(); }));
    //     }
    //
    //     // computes the coordinates of s_starting_point in the eigenbasis s_eigenvecs
    //     std::vector<S> const s_xis = rtk::details::dif::compute_xis(s_width_proj, s_eigenvals, s_eigenvecs, s_starting_point);
    //     // computes the coordinates of t_starting_point in the eigenbasis t_eigenvecs
    //     std::vector<S> const t_xis = rtk::details::dif::compute_xis(t_width_proj, t_eigenvals, t_eigenvecs, t_starting_point);
    //
    //     // computes the s_alphas
    //     std::vector<S> const s_ais = rtk::details::dif::compute_ais(s_width_proj, s_eigenvals);
    //     // computes the t_alphas
    //     std::vector<S> const t_ais = rtk::details::dif::compute_ais(t_width_proj, t_eigenvals);
    //
    //     auto const pd = [&projection, &s_xis, &t_xis, &s_eigenvecs, &t_eigenvecs, &transformation](scion::usize const i, scion::usize const j) -> rtk::vec<S, D> {
    //         return projection.project(transformation * (s_xis[i] * t_xis[j] * rtk::vector_product_tensor_flatten<S>(s_eigenvecs[i], t_eigenvecs[j]))).value();
    //     };
    //
    //     rtk::vec<S, D> normal = pd(1, 0).cross(pd(0, 1)).normalized();
    //
    //     S const i_det = pd(1, 0).dot(pd(1, 0)) * pd(0, 1).dot(pd(0, 1)) - std::pow(pd(1, 0).dot(pd(0, 1)), 2.0);
    //     S const m_cst = std::pow(pd(1, 1).dot(normal), 2.0);
    //
    //     S const k_constant = m_cst / i_det;
    //
    //     if (s_ais[2] > 2 && t_ais[2] > 2) {
    //         return std::copysign(std::numeric_limits<S>::infinity(), pd(2, 0).dot(normal) * pd(0, 2).dot(normal) / i_det);
    //     }
    //
    //     if (s_ais[2] > 2 && t_ais[2] == 2) {
    //         return std::copysign(std::numeric_limits<S>::infinity(), pd(2, 0).dot(normal) * pd(0, 2).dot(normal) / i_det);
    //     }
    //
    //     if (s_ais[2] > 2 && t_ais[2] < 2) {
    //         return std::unexpected{std::monostate{}};
    //     }
    //
    //     if (s_ais[2] == 2 && t_ais[2] > 2) {
    //         return std::copysign(std::numeric_limits<S>::infinity(), pd(2, 0).dot(normal) * pd(0, 2).dot(normal) / i_det);
    //     }
    //
    //     if (s_ais[2] == 2 && t_ais[2] == 2) {
    //         return (2.0 * pd(2, 0).dot(normal) * 2.0 * pd(0, 2).dot(normal) - m_cst) / i_det;
    //     }
    //
    //     if (s_ais[2] == 2 && t_ais[2] < 2) {
    //         return -k_constant;
    //     }
    //
    //     if (s_ais[2] < 2 && t_ais[2] > 2) {
    //         return std::unexpected{std::monostate{}};
    //     }
    //
    //     if (s_ais[2] < 2 && t_ais[2] == 2) {
    //         return -k_constant;
    //     }
    //
    //     if (s_ais[2] < 2 && t_ais[2] < 2) {
    //         return -k_constant;
    //     }
    //
    //     return std::unexpected{std::monostate{}};
    // }

    // template <scion::usize D, std::floating_point S>
    // auto sdcf_curvature_gaussian_attractor(rtk::projection<D> const& projection, scion::usize s_width_proj, std::span<const S> s_eigenvals, std::span<const rtk::vec<S>> s_eigenvecs, rtk::vec<S> s_starting_point, scion::usize t_width_proj, std::span<const S> t_eigenvals, std::span<const rtk::vec<S>> t_eigenvecs, rtk::vec<S> t_starting_point, S) -> std::expected<rtk::vertex<S, 3>, std::monostate> {
    //     scion::assert_debug(projection.size() == s_width_proj * t_width_proj);
    //
    //     scion::assert_debug(s_width_proj == s_eigenvals.size());
    //     scion::assert_debug(s_width_proj == s_eigenvecs.size());
    //     scion::assert_debug(std::ranges::all_of(s_eigenvecs, [s_width_proj](auto const& eigenvec) { return s_width_proj == eigenvec.size(); }));
    //
    //     scion::assert_debug(std::ranges::all_of(s_eigenvals, [](auto const& eigenval) { return !std::isnan(eigenval) && !std::isinf(eigenval); }));
    //     scion::assert_debug(std::ranges::all_of(s_eigenvecs, [](auto const& eigenvec) { return !eigenvec.hasNaN() && eigenvec.allFinite(); }));
    //
    //     scion::assert_debug(t_width_proj == t_eigenvals.size());
    //     scion::assert_debug(t_width_proj == t_eigenvecs.size());
    //     scion::assert_debug(std::ranges::all_of(t_eigenvecs, [t_width_proj](auto const& eigenvec) { return t_width_proj == eigenvec.size(); }));
    //
    //     scion::assert_debug(std::ranges::all_of(t_eigenvals, [](auto const& eigenval) { return !std::isnan(eigenval) && !std::isinf(eigenval); }));
    //     scion::assert_debug(std::ranges::all_of(t_eigenvecs, [](auto const& eigenvec) { return !eigenvec.hasNaN() && eigenvec.allFinite(); }));
    //
    //     auto const s_width_proj_i = static_cast<scion::isize>(s_width_proj);
    //     auto const t_width_proj_i = static_cast<scion::isize>(t_width_proj);
    //
    //     return std::unexpected<std::monostate>{{}};
    // }
}
