#pragma once

#include "scion/scion.hpp"

#include "rtk/core/projection.hpp"

#include "rtk/dif/utils.hpp"

#include "rtk/core/vertices.hpp"

#include "rtk/model/pifs.hpp"

#include "Eigen/Dense"

#include <span>
#include <optional>
#include <expected>
#include <ranges>
#include <cmath>
#include <algorithm>

/// \file rtk/dif/dcf.hpp
/// \brief provides a tool to compute the Differential Characteristic Function of a contractive transformation
/// \doi 10.5220/0012574800003660

namespace rtk::dif {
    /// \brief compute the DCF and its normal vector at t
    /// \param projection projection to use
    /// \param evals evals of the t-curve.
    /// \param evecs evecs of the t-curve.
    /// \param starting_point starting point of the t-curve
    /// \param t parameter
    /// \return coordinates of the point at (t) and its normal
    template <scion::usize D, bool StrictlySorted>
    auto dcf(rtk::projection const& projection, std::span<scion::f32 const> const evals, std::span<rtk::vecxf const> const evecs, rtk::vecxf const& starting_point, scion::f32 const t) -> std::expected<rtk::vertex<scion::f32, D, D>, std::monostate> {
        { // check preconditions
            if (projection.width() < 2) {
                throw std::logic_error("width_proj must be at least 2");
            }

            if (projection.width() != evals.size() || projection.width() != evecs.size()) {
                throw std::logic_error("inconsistent width_proj");
            }

            if (!rtk::pifs<StrictlySorted>::validate_evals(evals).has_value()) {
                return std::unexpected(std::monostate{});
            }
            if (!rtk::pifs<StrictlySorted>::validate_evecs(evecs).has_value()) {
                return std::unexpected(std::monostate{});
            }
        }

        // computes the coordinates of s_starting_point in the eigenbasis eigenvecs
        std::vector<scion::f32> xis = rtk::dif::compute_xis<scion::f32>(projection.width(), evecs, starting_point);
        xis[0] = 1.0f;

        // computes the alphas
        std::vector<scion::f32> ais = rtk::dif::compute_ais(projection.width(), evals);
        ais[0] = 0.0f;
        ais[1] = 1.0f;

        // computes the position in projective space
        rtk::vecxf const position_bn = std::ranges::fold_left_first(std::views::iota(0UZ, projection.width()) | std::views::transform([&](scion::usize const& index) -> rtk::vecxf { return xis[index] * evecs[index] * std::pow(t, ais[index]); } ), [](rtk::vecxf const& l, rtk::vecxf const& r) -> rtk::vecxf { return l + r; }).value();

        // computes the cross normals in projective space
        // rtk::vecxf const normal_s_bn = (std::ranges::fold_left_first(std::views::cartesian_product(std::views::iota(1u, s_width_proj), std::views::iota(0u, t_width_proj)) | std::views::transform([&](auto const& indices) -> rtk::vec<S> { return s_ais[std::get<0>(indices)] * s_xis[std::get<0>(indices)] * t_xis[std::get<1>(indices)] * rtk::vector_product_tensor_flatten<S>(s_eigenvecs[std::get<0>(indices)], t_eigenvecs[std::get<1>(indices)]) * std::pow(s, s_ais[std::get<0>(indices)] - 1.0f) * std::pow(t, t_ais[std::get<1>(indices)]); } ), [](rtk::vec<S> const& a, rtk::vec<S> const& b) -> rtk::vec<S> { return a + b; })).value();

        // computes the position in projected space
        rtk::vecf<D> const position_rn = projection.project<D>(position_bn).value();

        // computes the cross normals in projected space
        // rtk::vecf<D> const normal_s_rn = projection.project(normal_s_bn).value();

        // computes the normal in projected space
        // rtk::vecf<D> const normal_rn = normal_s_rn.cross(normal_t_rn).normalized();

        return rtk::vertex<scion::f32, D, D>{position_rn, position_rn};
    }

    template <scion::usize D, bool StrictlySorted>
    auto dcf(rtk::projection const& projection, std::span<scion::f32 const> const evals, std::span<rtk::vecxf const> const evecs, rtk::vecxf const& starting_point, scion::f32 const ta, scion::f32 const tz, scion::usize const count) -> std::expected<std::vector<rtk::vertex<scion::f32, D, D>>, std::monostate> {
        std::vector<rtk::vertex<scion::f32, D, D>> dcf;
        dcf.reserve(count);

        for (scion::usize i = 0; i < count; ++i) {
            scion::f32 const t = ta + (tz - ta) * (static_cast<scion::f32>(i) / static_cast<scion::f32>(count - 1));

            if (auto result = rtk::dif::dcf<D, StrictlySorted>(projection, evals, evecs, starting_point, t); result.has_value()) {
                dcf.emplace_back(result.value());
            } else {
                return std::unexpected(std::monostate{});
            }
        }

        return dcf;
    }

    // template <scion::usize D, std::floating_point S>
    // auto dcf_curvature(rtk::projection const& projection, Eigen::MatrixX<S> const& transform, std::span<const S> eigenvals, std::span<const Eigen::VectorX<S>> eigenvecs, Eigen::VectorX<S> starting_point, S const t) -> std::expected<S, std::monostate> {
    //     scion::assert_debug(projection.width() == transform.cols());
    //     scion::assert_debug(projection.width() == transform.rows());
    //     scion::assert_debug(projection.width() == eigenvals.size());
    //     scion::assert_debug(projection.width() == eigenvecs.size());
    //     scion::assert_debug(std::ranges::all_of(eigenvecs, [projection](auto const& eigenvec) { return projection.width() == eigenvec.size(); }));
    //
    //     scion::assert_debug(std::ranges::all_of(eigenvals, [](auto const& eigenval) { return !std::isnan(eigenval) && !std::isinf(eigenval); }));
    //     scion::assert_debug(std::ranges::all_of(eigenvecs, [](auto const& eigenvec) { return !eigenvec.hasNaN() && eigenvec.allFinite(); }));
    //
    //     return std::unexpected(std::monostate{});
    // }

//    template <typename S>
//    auto dcf(scion::usize const width_proj, Eigen::MatrixX<S> const& transformation, Eigen::VectorX<S> const& starting_point, S const t) -> std::expected<std::vector<Eigen::VectorX<S>>, std::monostate> {
//        scion::assert_debug(width_proj == transformation.cols());
//        scion::assert_debug(width_proj == transformation.rows());
//
//        // TODO check for nan and inf
//
//        std::vector<S> eigenvals;
//        eigenvals.reserve(width_proj);
//
//        std::vector<Eigen::VectorX<S>> eigenvecs;
//        eigenvecs.reserve(width_proj);
////
//        Eigen::EigenSolver<Eigen::MatrixX<S>> solver{transformation};
//        Eigen::VectorX<S> solver_output_eigenvalues = solver.eigenvalues().real();
//        Eigen::MatrixX<S> solver_output_eigenvectors = solver.eigenvectors().real();

//        if (solver_output_eigenvalues.size() != width_proj) {
//            return std::unexpected<std::monostate>{{}};
//        }
//
//        if (solver_output_eigenvectors.size() != width_proj) {
//            return std::unexpected<std::monostate>{{}};
//        }
//
//        for(scion::usize i = 0; i < width_proj; ++i) {
//            eigenvals.push_back(solver_output_eigenvalues[i]);
//            eigenvecs.push_back(solver_output_eigenvectors.col(i));
//        }
//
//        std::ranges::sort(std::ranges::zip_view(eigenvals, eigenvecs), [](auto const& lambda_1, auto const& lambda_2) -> bool { return std::abs(lambda_1) > std::abs(lambda_2); }, [](auto const& eigenpair) -> decltype(auto) { return std::get<0>(eigenpair); });
//
//        return rtk::dif::dcf(width_proj, transformation, eigenvals, eigenvecs, starting_point, t);
//    }
}
