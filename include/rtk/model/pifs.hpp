#pragma once

#include "scion/core/core.hpp"

#include "rtk/core/types.hpp"

#include "rtk/util/assert.hpp"

#include "rtk/core/projection.hpp"

#include "rtk/dif/utils.hpp"

#include <expected>
#include <iostream>
#include <memory>
#include <ranges>
#include <span>
#include <stdexcept>
#include <vector>

/// \file rtk/model/pifs.hpp
/// \brief provides a class to represent a Projected Iterated Function System

namespace rtk
{
    /// \brief represents a projected-IFS
    /// \tparam StrictlySorted specifies if the eigenvalues are to be sorted in strict descending order of modulus, or in descending order of modulus
    template <bool StrictlySorted = true>
    class pifs {
    private: // static private members

    protected: // static protected members

    public: // static public members
        /// \brief validates an eval, returning a short description of the error if it is invalid
        /// \param eval eval to validate
        /// \return short description of the error if the eval is invalid
        [[nodiscard]] static inline auto validate_eval(scion::f32 const eval) -> std::expected<void, std::string> {
            if (!rtk::is_val_valid(eval)) {
                return std::unexpected<std::string>("eval is not valid");
            }

            if (!rtk::is_val_eigenval_contraction(eval)) {
                return std::unexpected<std::string>("eval is not a contraction");
            }

            return {};
        }

        /// \brief validates an evals, returning a short description of the error if it is invalid
        /// \param evals eval to validate
        /// \return short description of the error if the evals is invalid
        [[nodiscard]] static inline auto validate_evals(std::span<scion::f32 const> const evals) -> std::expected<void, std::string> {
            auto result = evals | std::views::transform([](scion::f32 const eval) -> std::expected<void, std::string> {
                return rtk::pifs<StrictlySorted>::validate_eval(eval);
            }) | std::views::filter([](std::expected<void, std::string> const & expected) -> bool {
                return !expected.has_value();
            });
            if (auto it = result.begin(); it != result.end()) {
                return *it;
            }

            if (evals[0] < 1.0f - rtk::EPS_MP<scion::f32> || 1.0f + rtk::EPS_MP<scion::f32> < evals[0]) {
                return std::unexpected<std::string>("first eval is not 1");
            }

            if constexpr (StrictlySorted) {
                for (scion::usize i = 1; i < evals.size(); i++) {
                    if (std::abs(evals[i - 1]) <= std::abs(evals[i])) {
                        return std::unexpected<std::string>("evals is not sorted");
                    }
                }
            } else {
                for (scion::usize i = 1; i < evals.size(); i++) {
                    if (std::abs(evals[i - 1]) < std::abs(evals[i])) {
                        return std::unexpected<std::string>("evals is not sorted");
                    }
                }
            }

            if (scion::usize const width = evals.size(); width < 2) {
                return std::unexpected<std::string>("evals width must be at least 2");
            }

            return {};
        }

        /// \brief validates an evec, returning a short description of the error if it is invalid
        /// \param evec evecs to validate
        /// \return short description of the error if the evec is invalid
        [[nodiscard]] static inline auto validate_evec(rtk::vecxf const & evec) -> std::expected<void, std::string> {
            if (!rtk::is_vec_valid(evec)) {
                return std::unexpected<std::string>("evec is not valid");
            }

            return {};
        }

        /// \brief validates an evec, returning a short description of the error if it is invalid
        /// \param evec evecs to validate
        /// \return short description of the error if the evec is invalid
        [[nodiscard]] static inline auto validate_evec_ptn(rtk::vecxf const & evec) -> std::expected<void, std::string> {
            if (!rtk::is_vec_pnt(evec)) {
                return std::unexpected<std::string>("evec does not describe a point");
            }

            return {};
        }

        /// \brief validates an evec, returning a short description of the error if it is invalid
        /// \param evec evecs to validate
        /// \return short description of the error if the evec is invalid
        [[nodiscard]] static inline auto validate_evec_vec(rtk::vecxf const & evec) -> std::expected<void, std::string> {
            if (!rtk::is_vec_vec(evec)) {
                return std::unexpected<std::string>("evec does not describe a vector");
            }

            return {};
        }

        /// \brief validates an evecs, returning a short description of the error if it is invalid
        /// \param evecs evecs to validate
        /// \return short description of the error if the evecs is invalid
        [[nodiscard]] static inline auto validate_evecs(std::span<rtk::vecxf const> const evecs) -> std::expected<void, std::string> {
            scion::usize const width_proj = evecs.size();
            scion::isize const width_proj_i = static_cast<scion::isize>(width_proj);

            {
                auto result = evecs | std::views::transform([](rtk::vecxf const & evec) -> std::expected<void, std::string> {
                    return rtk::pifs<StrictlySorted>::validate_evec(evec);
                }) | std::views::filter([](std::expected<void, std::string> const & expected) -> bool {
                    return !expected.has_value();
                });
                if (auto it = result.begin(); it != result.end()) {
                    return *it;
                }
            }

            {
                if (auto const result = rtk::pifs<StrictlySorted>::validate_evec_ptn(evecs[0]); !result.has_value()) {
                    return result;
                }
            }

            {
                auto result = evecs | std::views::drop(1) | std::views::transform([](rtk::vecxf const & evec) -> std::expected<void, std::string> {
                    return rtk::pifs<StrictlySorted>::validate_evec_vec(evec);
                }) | std::views::filter([](std::expected<void, std::string> const & expected) -> bool {
                    return !expected.has_value();
                });
                if (auto it = result.begin(); it != result.end()) {
                    return *it;
                }
            }

            {
                if (width_proj < 2) {
                    return std::unexpected<std::string>("evecs width must be at least 2");
                }

                auto result = evecs | std::views::filter([width_proj_i](rtk::vecxf const & evec) {
                    return width_proj_i != evec.size();
                });
                if (auto it = result.begin(); it != result.end()) {
                    return std::unexpected<std::string>("evecs has inconsistent width");
                }
            }

            {
                rtk::matxf V = rtk::matxf::Zero(width_proj_i, width_proj_i);

                for (scion::isize d = 0; d < width_proj_i; ++d) {
                    V.col(d) = evecs[d];
                }

                Eigen::FullPivLU<rtk::matxf> lu(V);
                if (lu.rank() != V.rows()) {
                    return std::unexpected<std::string>("evecs does not form an invertible matrix");
                }
            }

            return {};
        }

        /// \brief validates a tran, returning a short description of the error if it is invalid
        /// \param tran tran to validate
        /// \return short description of the error if the tran is invalid
        [[nodiscard]] static inline auto validate_tran(rtk::matxf const & tran) -> std::expected<void, std::string> {
            if (!rtk::is_mat_valid(tran)) {
                return std::unexpected<std::string>("tran has invalid value");
            }

            if (!rtk::is_mat_semi_right_stochastic(tran)) {
                return std::unexpected<std::string>("tran is not semi-right-stochastic");
            }

            return {};
        }

        /// \brief computes a transformation from its eigenvalues and eigenvectors
        /// \param evals eigenvalues to use
        /// \param evecs eigenvectors to use
        /// \return transformation
        /// \errors incoherent projection width
        /// \errors first eigenvalue is not 1
        /// \errors an eigenvalue is not valid
        /// \errors eigenvalues are not sorted
        /// \errors an eigenvector is not valid
        /// \errors first eigenvector does not represent a point
        /// \errors a subsequent eigenvector is not a vector
        /// \errors eigenvector matrix is not invertible
        ///
        /// The returned matrix is a semi-right-stochastic matrix.
        [[nodiscard]] static auto compute_transformation(std::span<scion::f32 const> const evals, std::span<rtk::vecxf const> const evecs) -> std::expected<rtk::matxf, std::invalid_argument> {
            scion::usize const width_proj = evals.size();
            scion::isize const width_proj_i = static_cast<scion::isize>(width_proj);

            if (width_proj < 2) {
                return std::unexpected<std::invalid_argument>("projection width must be at least 2");
            }

            if (evals.size() != evecs.size()) {
                return std::unexpected<std::invalid_argument>("incoherent projection width");
            }

            if (auto result = rtk::pifs<StrictlySorted>::validate_evals(evals); !result.has_value()) {
                return std::unexpected<std::invalid_argument>(std::move(result.error()));
            }

            if (auto result = rtk::pifs<StrictlySorted>::validate_evecs(evecs); !result.has_value()) {
                return std::unexpected<std::invalid_argument>(std::move(result.error()));
            }

            rtk::matxf L = rtk::matxf::Zero(width_proj_i, width_proj_i);
            rtk::matxf V = rtk::matxf::Zero(width_proj_i, width_proj_i);

            for (scion::isize d = 0; d < width_proj_i; ++d) {
                L(d, d) = evals[d];
                V.col(d) = evecs[d];
            }

            rtk::matxf tran = V * L * V.inverse();

            if (auto result = rtk::pifs<StrictlySorted>::validate_tran(tran); !result.has_value()) {
                return std::unexpected<std::invalid_argument>(result.error());
            }

            return tran;
        }

        /// \brief computes evals and evecs from a transformation
        /// \param tran transformation to use
        /// \return evals and evecs
        static auto compute_eigenproperties(rtk::matxf const & tran) -> std::expected<std::pair<std::vector<scion::f32>, std::vector<rtk::vecxf>>, std::invalid_argument> {
            scion::usize const width_proj = tran.cols();
            scion::isize const width_proj_i = static_cast<scion::isize>(width_proj);

            if (tran.cols() != tran.rows()) {
                return std::unexpected<std::invalid_argument>("width_proj is inconsistent");
            }

            if (width_proj < 2) {
                return std::unexpected<std::invalid_argument>("width_proj must be at least 2");
            }

            if (auto result = rtk::pifs<StrictlySorted>::validate_tran(tran); !result.has_value()) {
                return std::unexpected<std::invalid_argument>(result.error());
            }

            std::vector<scion::f32> evals;
            std::vector<rtk::vecxf> evecs;

            Eigen::EigenSolver<rtk::matxf> const solver{tran};

            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("eigensolver failed");
            }

            rtk::vecxf solver_output_eigenvals = solver.eigenvalues().real();
            rtk::matxf solver_output_eigenvecs = solver.eigenvectors().real();

            std::vector<scion::isize> indices{};
            for (scion::isize index = 0; index < width_proj_i; ++index) indices.emplace_back(index);

            std::ranges::sort(indices, std::greater{}, [&solver_output_eigenvals](scion::isize const index) -> scion::f32 {
                return std::abs(solver_output_eigenvals(index));
            });

            for (scion::isize index = 0; index < width_proj_i; ++index) {
                evals.push_back(solver_output_eigenvals[indices[index]]);
                evecs.push_back(solver_output_eigenvecs.col(indices[index]));
            }

            {
                scion::f32 acc = 0.0f;

                for (scion::isize index_proj = 0; index_proj < width_proj_i; ++index_proj) {
                    acc += evecs[0](index_proj);
                }

                if (acc == 0.0f) {
                    return std::unexpected<std::invalid_argument>("first eigenvector does not describe a point");
                }

                evals[0] = 1.0f;
                evecs[0] /= acc;
            }

            if (evals.size() != width_proj) {
                throw std::logic_error("incoherent projection width");
            }

            if (auto result = rtk::pifs<StrictlySorted>::validate_evals(evals); !result.has_value()) {
                return std::unexpected<std::invalid_argument>(result.error());
            }

            if (auto result = rtk::pifs<StrictlySorted>::validate_evecs(evecs); !result.has_value()) {
                return std::unexpected<std::invalid_argument>(result.error());
            }

            return std::make_pair(std::move(evals), std::move(evecs));
        }

        /// \brief makes an empty pifs with the specified width_tran and width_proj
        /// \tparam Ptr std::unique_ptr or std::shared_ptr if you want a dynamically allocated pifs
        /// \param width_tran width_tran to use
        /// \param width_proj width_proj to use
        /// \return newly created pifs, or error
        template <template<typename> typename Ptr = std::type_identity>
            requires (std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::type_identity<rtk::pifs<StrictlySorted>>> || std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::unique_ptr<rtk::pifs<StrictlySorted>>> || std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::shared_ptr<rtk::pifs<StrictlySorted>>>)
        static auto make(scion::usize const width_tran, scion::usize const width_proj) -> std::expected<std::conditional_t<std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::type_identity<rtk::pifs<StrictlySorted>>>, rtk::pifs<StrictlySorted>, Ptr<rtk::pifs<StrictlySorted>>>, std::invalid_argument> {
            if (width_tran < 2) {
                return std::unexpected(std::invalid_argument{"transformation width must be greater or equal than 2"});
            }

            if (width_proj < 2) {
                return std::unexpected(std::invalid_argument{"projection width must be greater or equal than 2"});
            }

            std::vector<std::vector<scion::f32>> eigenvals = std::vector(width_tran, std::vector(width_proj, 1.0f));
            std::vector<std::vector<rtk::vecxf>> eigenvecs = std::vector<std::vector<rtk::vecxf>>(width_tran, std::vector<rtk::vecxf>(width_proj, rtk::vecxf::Zero(width_proj)));
            std::vector<rtk::matxf> trans = std::vector<rtk::matxf>(width_tran, rtk::matxf::Identity(width_proj, width_proj));

            if constexpr (std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::type_identity<rtk::pifs<StrictlySorted>>>) {
                return rtk::pifs<StrictlySorted>{
                    M{
                        width_tran,
                        width_proj,
                        std::move(eigenvals),
                        std::move(eigenvecs),
                        std::move(trans),
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::vector<rtk::vecxf>{},
                        0UZ,
                        false,
                        std::vector<rtk::vecxf>{},
                        4,
                        0,
                        {},
                        {},
                    }
                };
            } else {
                return Ptr<rtk::pifs<StrictlySorted>>(new rtk::pifs<StrictlySorted>{
                    M{
                        width_tran,
                        width_proj,
                        std::move(eigenvals),
                        std::move(eigenvecs),
                        std::move(trans),
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::vector<rtk::vecxf>{},
                        0UZ,
                        false,
                        std::vector<rtk::vecxf>{},
                        4,
                        0,
                        {},
                        {},
                    }
                });
            }
        }

        /// \brief makes a new pifs with the specified eigenvals and eigenvecs
        /// \tparam Ptr std::unique_ptr or std::shared_ptr if you want a dynamically allocated pifs
        /// \param eigenvals eigenvals to use
        /// \param eigenvecs eigenvecs to use
        /// \return newly created pifs, or error
        /// \throws std::logic_error if eigenvals and eigenvecs are particularly unsanitized
        template <template<typename> typename Ptr = std::type_identity>
            requires (std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::type_identity<rtk::pifs<StrictlySorted>>> || std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::unique_ptr<rtk::pifs<StrictlySorted>>> || std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::shared_ptr<rtk::pifs<StrictlySorted>>>)
        static auto make(std::vector<std::vector<scion::f32>> && eigenvals, std::vector<std::vector<rtk::vecxf>> && eigenvecs) -> std::expected<std::conditional_t<std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::type_identity<rtk::pifs<StrictlySorted>>>, rtk::pifs<StrictlySorted>, Ptr<rtk::pifs<StrictlySorted>>>, std::invalid_argument> {
            scion::usize const width_tran = eigenvecs.size();
            { // ensures coherent width_tran
                if (eigenvals.size() != eigenvecs.size()) {
                    return std::unexpected<std::invalid_argument>("inconsistent width_tran");
                }

                if (width_tran < 2) {
                    return std::unexpected<std::invalid_argument>("width_tran must be at least 2");
                }
            }

            scion::usize const width_proj = eigenvecs[0].size();
            { // ensures coherent width_proj
                for (scion::usize index_tran = 0; index_tran < width_tran; ++index_tran) {
                    if (eigenvals[index_tran].size() != width_proj || eigenvecs[index_tran].size() != width_proj) {
                        return std::unexpected<std::invalid_argument>("inconsistent width_proj");
                    }
                }

                if (width_proj < 2) {
                    return std::unexpected<std::invalid_argument>("width_proj must be at least 2");
                }
            }

            { // ensures consistent evec direction
                rtk::dif::compute_xis_force_positive<scion::f32>(width_proj, eigenvecs.front(), eigenvecs.back().front());
                rtk::dif::compute_xis_force_positive<scion::f32>(width_proj, eigenvecs.back(), eigenvecs.front().front());
            }

            auto range = std::views::zip_transform([](std::span<scion::f32 const> const evals, std::span<rtk::vecxf const> const evecs) -> std::expected<rtk::matxf, std::invalid_argument> {
                return rtk::pifs<StrictlySorted>::compute_transformation(evals, evecs);
            }, eigenvals, eigenvecs);

            std::vector<rtk::matxf> trans{};
            for (const auto & result : range) {
                if (!result.has_value()) {
                    return std::unexpected<std::invalid_argument>(result.error());
                }
                trans.emplace_back(result.value());
            }

            if constexpr (std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::type_identity<rtk::pifs<StrictlySorted>>>) {
                return rtk::pifs<StrictlySorted>{
                    {
                        width_tran,
                        width_proj,
                        std::move(eigenvals),
                        std::move(eigenvecs),
                        std::move(trans),
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::vector<rtk::vecxf>{},
                        0UZ,
                        true,
                        std::vector<rtk::vecxf>{},
                        4,
                        0,
                        {},
                        {},
                    }
                };
            } else {
                return Ptr<rtk::pifs<StrictlySorted>>(new rtk::pifs<StrictlySorted>{
                    {
                        width_tran,
                        width_proj,
                        std::move(eigenvals),
                        std::move(eigenvecs),
                        std::move(trans),
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::vector<rtk::vecxf>{},
                        0UZ,
                        true,
                        std::vector<rtk::vecxf>{},
                        4,
                        0,
                        {},
                        {},
                    }
                });
            }
        }

        /// \brief makes a new pifs with the specified trans
        /// \tparam Ptr std::unique_ptr or std::shared_ptr if you want a dynamically allocated pifs
        /// \param trans trans to use
        /// \return newly created pifs, or error
        /// \throws std::logic_error if trans is particularly unsanitized
        template <template<typename> typename Ptr = std::type_identity>
            requires (std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::type_identity<rtk::pifs<StrictlySorted>>> || std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::unique_ptr<rtk::pifs<StrictlySorted>>> || std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::shared_ptr<rtk::pifs<StrictlySorted>>>)
        static auto make(std::vector<rtk::matxf> && trans) -> std::expected<std::conditional_t<std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::type_identity<rtk::pifs<StrictlySorted>>>, rtk::pifs<StrictlySorted>, Ptr<rtk::pifs<StrictlySorted>>>, std::invalid_argument> {
            scion::usize const width_tran = trans.size();
            scion::usize const width_proj = trans[0].cols();

            if (width_tran < 2) {
                return std::unexpected<std::invalid_argument>("width_tran must be at least 2");
            }

            if (width_proj < 2) {
                return std::unexpected<std::invalid_argument>("width_proj must be at least 2");
            }

            std::vector<std::vector<scion::f32>> eigenvals;
            std::vector<std::vector<rtk::vecxf>> eigenvecs;

            for (auto const & tran : trans) {
                auto [evals, evecs] = rtk::pifs<StrictlySorted>::compute_eigenproperties(tran).value();

                eigenvals.emplace_back(std::move(evals));
                eigenvecs.emplace_back(std::move(evecs));
            }

            rtk::dif::compute_xis_force_positive<scion::f32>(width_proj, eigenvecs.front(), eigenvecs.back().front());
            rtk::dif::compute_xis_force_positive<scion::f32>(width_proj, eigenvecs.back(), eigenvecs.front().front());

            if constexpr (std::is_same_v<Ptr<rtk::pifs<StrictlySorted>>, std::type_identity<rtk::pifs<StrictlySorted>>>) {
                return rtk::pifs<StrictlySorted>{
                    {
                        width_tran,
                        width_proj,
                        std::move(eigenvals),
                        std::move(eigenvecs),
                        std::move(trans),
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::vector<rtk::vecxf>{},
                        0UZ,
                        true,
                        std::vector<rtk::vecxf>{},
                        4,
                        0,
                        {},
                        {},
                    }
                };
            } else {
                return Ptr<rtk::pifs<StrictlySorted>>(new rtk::pifs<StrictlySorted>{
                    {
                        width_tran,
                        width_proj,
                        std::move(eigenvals),
                        std::move(eigenvecs),
                        std::move(trans),
                        std::nullopt,
                        std::nullopt,
                        std::nullopt,
                        std::vector<rtk::vecxf>{},
                        0UZ,
                        true,
                        std::vector<rtk::vecxf>{},
                        4,
                        0,
                        {},
                        {},
                    }
                });
            }
        }

    private: // static private functions

    protected: // static protected functions

    public: // static public functions

    private: // private members

    protected: // protected members
        /// \brief member data
        struct M {
            /// \brief number of transformations of the PIFS
            ///
            /// \def width_tran: number of transformations of the PIFS
            scion::usize width_tran;

            /// \brief size of the projection
            ///
            /// \def width_proj: size of the projection
            scion::usize width_proj;

            /// \brief set of all eigenvalues of all transformations
            ///
            /// \def eigenvals: set of all eigenvalues of all transformations
            /// \def evals:     set of all eigenvalues of a transformations
            /// \def eval:      a single eigenvalue
            std::vector<std::vector<scion::f32>> eigenvals;

            /// \brief set of all eigenvectors of all transformations
            ///
            /// \def eigenvecs: set of all eigenvectors of all transformations
            /// \def evecs:     set of all eigenvectors of a transformations
            /// \def evec:      a single eigenvector
            std::vector<std::vector<rtk::vecxf>> eigenvecs;

            /// \brief set of all transformations
            ///
            /// \def trans: set of all transformations
            /// \def tran:  a single transformation
            std::vector<rtk::matxf> trans;

            /// \brief projection to use for a modeling space in arbitrary dimension
            std::optional<std::shared_ptr<rtk::projection>> projection;

            /// \brief projection to use for $\mathbb{R}^2$
            std::optional<std::shared_ptr<rtk::projection>> projection_d2;

            /// \brief projection to use for $\mathbb{R}^3$
            std::optional<std::shared_ptr<rtk::projection>> projection_d3;

            /// \brief attractor in barycentric space
            std::vector<rtk::vecxf> attractor_bn;

            /// \brief attractor in modeling space
            ///
            /// scion::usize              attractor has not yet been projected, will be computed with the specified dimension if a projection is available
            /// std::vector<rtk::vecf<0>> attractor was projected using m.projection
            /// std::vector<rtk::vecf<2>> attractor was projected using m.projection_d2
            /// std::vector<rtk::vecf<3>> attractor was projected using m.projection_d3
            std::variant<scion::usize, std::vector<rtk::vecf<0>>, std::vector<rtk::vecf<2>>, std::vector<rtk::vecf<3>>> attractor_rn;

            /// \brief flag for valid/invalid state
            bool state;

            /// \brief primitive to use for computing the attractor in barycentric space
            std::vector<vecxf> primitive;

            /// \brief number of iteration to apply when computing the attractor
            scion::usize iterations = 4;

            /// \brief state_counter of the projection that was used to compute the attractor in modeling space
            scion::usize projection_counter = 0;

            /// \brief aggregate of the update requested
            struct {
                /// \brief true if the transformations are to be recomputed
                bool transformations = false;

                /// \brief true if the eigenproperties are to be recomputed
                bool eigenproperties = false;

                /// \brief true if the attractor in barycentric space is to be recomputed
                bool attractor_bn = false;

                /// \brief true if the attractor in modeling space is to be recomputed
                bool attractor_rn = false;
            } update_requests;

            /// \brief aggregate of the update processed
            struct {
                /// \brief true if the transformations update was processed
                bool transformations = false;

                /// \brief true if the eigenproperties update was processed
                bool eigenproperties = false;

                /// \brief true if the attractor in barycentric space update was processed
                bool attractor_bn = false;

                /// \brief true if the attractor in modeling space update was processed
                bool attractor_rn = false;
            } update_processed;
        } m;

    public: // public members

    private: // private methods

    protected: // protected methods
        /// \brief constructor for a PIFS
        /// \param m member data
        explicit pifs(rtk::pifs<StrictlySorted>::M m) noexcept : m{std::move(m)} {
            rtk::log::trace("init", "ctor for rtk::pifs@{}", static_cast<void *>(this));
        }

        /// \brief computes trans from eigenvals and eigenvecs
        void update_trans() {
            rtk::log::debug("update", "trans of {}", static_cast<void *>(this));

            for (scion::usize index_tran = 0; index_tran < this->width_tran(); ++index_tran) {
                if (auto result = rtk::pifs<StrictlySorted>::compute_transformation(this->evals(index_tran), this->evecs(index_tran)); result.has_value()) {
                    m.trans[index_tran] = result.value();
                } else {
                    rtk::log::error("update", "could not update trans of {}: {}", static_cast<void *>(this), result.error().what());

                    this->invalidate();

                    return;
                }
            }

            m.update_processed.transformations = true;
        }

        /// \brief computes the eigenvalues and eigenvectors from the transformations
        void update_eigen() {
            rtk::log::debug("update", "eigen of {}", static_cast<void *>(this));

            for (scion::usize index_tran = 0; index_tran < this->width_tran(); ++index_tran) {
                if (auto result = rtk::pifs<StrictlySorted>::compute_eigenproperties(this->tran(index_tran)); result.has_value()) {
                    auto eigen = result.value();

                    m.eigenvals[index_tran] = eigen.first;
                    m.eigenvecs[index_tran] = eigen.second;
                } else {
                    rtk::log::error("update", "could not update eigen of {}: {}", static_cast<void *>(this), result.error().what());

                    this->invalidate();

                    return;
                }
            }

            m.update_processed.eigenproperties = true;
        }

        /// \brief computes the attractor in barycentric space
        void update_attractor_bn() {
            rtk::log::debug("update", "attractor_bn of {}", static_cast<void *>(this));

            m.attractor_bn.clear();

            if (m.primitive.empty()) {
                rtk::log::warn("update", "attractor_bn of {}: empty primitive", static_cast<void *>(this));

                m.update_processed.attractor_bn = true;

                return;
            }

            auto solver = [this](auto && solver, rtk::matxf const & transformation, scion::usize const iteration) {
                if (iteration == 0) {
                    for (rtk::vecxf const & prim : m.primitive) {
                        m.attractor_bn.emplace_back(transformation * prim);
                    }

                    return;
                }

                for (rtk::matxf const & matrix : m.trans) {
                    solver(solver, transformation * matrix, iteration - 1);
                }
            };

            solver(solver, rtk::matxf::Identity(this->width_proj(), this->width_proj()), m.iterations);

            rtk::log::debug("update", "attractor_bn of {}: {} vertices", static_cast<void *>(this), m.attractor_bn.size());

            m.update_processed.attractor_bn = true;
        }

        /// \brief computes the attractor in modeling space using the projection for the appropriate dimension
        void update_attractor_rn() {
            rtk::log::debug("update", "attractor_rn of {}", static_cast<void *>(this));

            switch (this->dimension()) {
            case 0: {
                if (std::holds_alternative<std::vector<rtk::vecf<0>>>(m.attractor_rn)) {
                    std::get<std::vector<rtk::vecf<0>>>(m.attractor_rn).clear();
                }

                if (this->has_projection<0>()) {
                    if (auto result = m.projection.value()->template project<0>(m.attractor_bn); result.has_value()) {
                        m.attractor_rn.template emplace<std::vector<rtk::vecxf>>(std::move(result.value()));
                    } else {
                        rtk::log::error("update", "could not update attractor_rn of {}: {}", static_cast<void *>(this), result.error().what());

                        this->invalidate();

                        return;
                    }

                    m.projection_counter = this->projection<0>()->counter();
                }

                break;
            }
            case 2: {
                if (std::holds_alternative<std::vector<rtk::vecf<2>>>(m.attractor_rn)) {
                    std::get<std::vector<rtk::vecf<2>>>(m.attractor_rn).clear();
                }

                if (this->has_projection<2>()) {
                    if (auto result = this->projection<2>()->template project<2>(m.attractor_bn); result.has_value()) {
                        m.attractor_rn.template emplace<std::vector<rtk::vec2f>>(std::move(result.value()));
                    } else {
                        rtk::log::error("update", "could not update attractor_rn of {}: {}", static_cast<void *>(this), result.error().what());

                        this->invalidate();

                        return;
                    }

                    m.projection_counter = this->projection<2>()->counter();
                }

                break;
            }
            case 3: {
                if (std::holds_alternative<std::vector<rtk::vecf<3>>>(m.attractor_rn)) {
                    std::get<std::vector<rtk::vecf<3>>>(m.attractor_rn).clear();
                }

                if (this->has_projection<3>()) {
                    if (auto result = this->projection<3>()->template project<3>(m.attractor_bn); result.has_value()) {
                        m.attractor_rn.template emplace<std::vector<rtk::vec3f>>(std::move(result.value()));
                    } else {
                        rtk::log::error("update", "could not update attractor_rn of {}: {}", static_cast<void *>(this), result.error().what());

                        this->invalidate();

                        return;
                    }

                    m.projection_counter = this->projection<3>()->counter();
                }

                break;
            }
            default:
                throw std::logic_error("dimension not supported");
            }

            m.update_processed.attractor_rn = true;
        }

    public: // public methods
        /// \brief copy constructor
        /// \param other other model to copy
        pifs(rtk::pifs<StrictlySorted> const & other): m{other.m} {
            rtk::log::trace("copy", "ctor for rtk::pifs@{} (from rtk::pifs@{})", static_cast<void *>(this), static_cast<void const *>(&other));
        }

        /// \brief move constructor
        /// \param other other model to move
        pifs(rtk::pifs<StrictlySorted> && other) noexcept: m{std::move(other.m)} {
            rtk::log::trace("move", "ctor for rtk::pifs@{} (from rtk::pifs@{})", static_cast<void *>(this), static_cast<void *>(&other));
        }

        /// \brief copy assignment operator
        /// \param other other model to copy
        /// \return instance to which it was copied
        rtk::pifs<StrictlySorted> & operator=(rtk::pifs<StrictlySorted> const & other) {
            this->m = other.m;

            rtk::log::trace("copy", "otor for rtk::pifs@{} (from rtk::pifs@{})", static_cast<void *>(this), static_cast<void const *>(&other));

            return *this;
        }

        /// \brief move assignment operator
        /// \param other other model to move
        /// \return instance to which it was moved
        rtk::pifs<StrictlySorted> & operator=(rtk::pifs<StrictlySorted> && other) noexcept {
            this->m = std::move(other.m);

            rtk::log::trace("move", "otor for rtk::pifs@{} (from rtk::pifs@{})", static_cast<void *>(this), static_cast<void *>(&other));

            return *this;
        }

        /// \brief destructor for a pifs
        virtual ~pifs() noexcept {
            rtk::log::trace("drop", "dtor for rtk::pifs@{}", static_cast<void *>(this));
        }

        /// \brief flags the transformations for update
        ///
        /// Transformations require an update if the eigenvalues or eigenvectors were changed.
        void request_update_transformations() noexcept {
            m.update_requests.transformations = true;
        }

        /// \brief flags the eigenproperties for update
        ///
        /// Eigenproperties require an update if the transformations were changed.
        void request_update_eigenproperties() noexcept {
            m.update_requests.eigenproperties = true;
        }

        /// \brief flags the attractor in barycentric space for update
        ///
        /// Barycentric attractor requires an update if the transformations were changed.
        void request_update_attractor_bn() noexcept {
            m.update_requests.attractor_bn = true;
        }

        /// \brief flags the attractor in modeling space for update
        ///
        /// Modeling attractor requires an update if the barycentric attractor was changed, or if the projection was changed.
        void request_update_attractor_rn() noexcept {
            m.update_requests.attractor_rn = true;
        }

        /// \brief updates the fields that were flagged for an update, propagating the updates as necessary
        virtual void update() {
            { // resets the update_processed
                m.update_processed.transformations = false;
                m.update_processed.eigenproperties = false;
                m.update_processed.attractor_bn = false;
                m.update_processed.attractor_rn = false;
            }

            { // check for projection-based update
                if (this->dimension() == 0) {
                    if (this->has_projection<0>() && this->projection<0>()->counter() != m.projection_counter) {
                        this->request_update_attractor_rn();
                    }
                }

                if (this->dimension() == 2) {
                    if (this->has_projection<2>() && this->projection<2>()->counter() != m.projection_counter) {
                        this->request_update_attractor_rn();
                    }
                }

                if (this->dimension() == 3) {
                    if (this->has_projection<3>() && this->projection<3>()->counter() != m.projection_counter) {
                        this->request_update_attractor_rn();
                    }
                }
            }

            this->validate<rtk::log::debug_t>();

            { // processes update_requests
                if (this->valid() && m.update_requests.transformations) {
                    rtk::dif::compute_xis_force_positive<scion::f32>(this->width_proj(), this->evecs_mut(0), this->evec(this->width_tran() - 1, 0));
                    rtk::dif::compute_xis_force_positive<scion::f32>(this->width_proj(), this->evecs_mut(this->width_tran() - 1), this->evec(0, 0));

                    this->update_trans();
                } else if (this->valid() && m.update_requests.eigenproperties) {
                    this->update_eigen();
                }

                if (this->valid() && (m.update_requests.attractor_bn || m.update_processed.transformations || m.update_processed.eigenproperties)) {
                    this->update_attractor_bn();
                }

                if (this->valid() && (m.update_requests.attractor_rn || m.update_processed.attractor_bn)) {
                    this->update_attractor_rn();
                }
            }

            { // resets the update_requests
                m.update_requests.transformations = false;
                m.update_requests.eigenproperties = false;
                m.update_requests.attractor_bn = false;
                m.update_requests.attractor_rn = false;
            }
        }

        /// \brief returns true if the queried projection is populated
        /// \tparam Dim dimension of the projection to query (0 for dynamic)
        /// \return true if the queried projection is populated
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        [[nodiscard]] inline auto has_projection() const noexcept -> bool {
            if constexpr (Dim == 0) {
                return m.projection.has_value();
            } else if constexpr (Dim == 2) {
                return m.projection_d2.has_value();
            } else if constexpr (Dim == 3) {
                return m.projection_d3.has_value();
            } else {
                static_assert(false);

                std::unreachable();
            }
        }

        /// \brief returns the queried projection
        /// \tparam Dim dimension of the projection to get (0 for dynamic)
        /// \return queried projection
        /// \throws std::bad_variant_access if the queried projection is not populated
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        [[nodiscard]] inline auto projection() const noexcept -> std::shared_ptr<const rtk::projection> {
            if constexpr (Dim == 0) {
                return m.projection.value();
            } else if constexpr (Dim == 2) {
                return m.projection_d2.value();
            } else if constexpr (Dim == 3) {
                return m.projection_d3.value();
            } else {
                static_assert(false);

                std::unreachable();
            }
        }

        /// \brief sets a new projection
        /// \tparam Dim dimension of the projection to ste (0 for dynamic)
        /// \param projection projection to use
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        void set_projection(std::shared_ptr<rtk::projection> const & projection) {
            if (Dim != 0 && projection->dimension() != Dim) {
                throw std::logic_error("dimension mismatch");
            }

            if (this->width_proj() != projection->width()) {
                throw std::logic_error("width mismatch");
            }

            if constexpr (Dim == 0) {
                m.projection = projection;
                m.projection_counter = projection->counter();

                this->request_update_attractor_rn();
            } else if constexpr (Dim == 2) {
                m.projection_d2 = projection;
                m.projection_counter = projection->counter();

                this->request_update_attractor_rn();
            } else if constexpr (Dim == 3) {
                m.projection_d3 = projection;
                m.projection_counter = projection->counter();

                this->request_update_attractor_rn();
            } else {
                static_assert(false);

                std::unreachable();
            }
        }

        /// \brief returns the width of the transformations
        /// \return width of the transformations
        [[nodiscard]] inline auto width_tran() const noexcept -> scion::usize {
            return m.width_tran;
        }

        /// \brief returns the width of the projections
        /// \return width of the projections
        [[nodiscard]] inline auto width_proj() const noexcept -> scion::usize {
            return m.width_proj;
        }

        /// \brief returns the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \return tran at index index_tran
        /// \throws std::logic_error if no such tran exists
        [[nodiscard]] inline auto tran(scion::usize const index_tran) const -> rtk::matxf const & {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            return m.trans[index_tran];
        }

        /// \brief returns the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \return tran at index index_tran
        /// \throws std::logic_error if no such tran exists
        /// \note client code is responsible to maintain the invariants of the object if changes are made through the returned reference
        [[nodiscard]] inline auto tran_mut(scion::usize const index_tran) noexcept -> rtk::matxf & {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            return m.trans[index_tran];
        }

        /// \brief returns the evals of the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \return evals of the tran at index index_tran
        /// \throws std::logic_error if no such tran exists
        [[nodiscard]] inline auto evals(scion::usize const index_tran) const -> std::span<const scion::f32> {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            return m.eigenvals[index_tran];
        }

        /// \brief returns the evals of the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \return evals of the tran at index index_tran
        /// \throws std::logic_error if no such tran exists
        /// \note client code is responsible to maintain the invariants of the object if changes are made through the returned reference
        [[nodiscard]] inline auto evals_mut(scion::usize const index_tran) -> std::span<scion::f32> {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            return m.eigenvals[index_tran];
        }

        /// \brief returns the evecs of the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \return evecs of the tran at index index_tran
        /// \throws std::logic_error if no such tran exists
        [[nodiscard]] inline auto evecs(scion::usize const index_tran) const -> std::span<const rtk::vecxf> {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            return m.eigenvecs[index_tran];
        }

        /// \brief returns the evecs of the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \return evecs of the tran at index index_tran
        /// \throws std::logic_error if no such tran exists
        /// \note client code is responsible to maintain the invariants of the object if changes are made through the returned reference
        [[nodiscard]] inline auto evecs_mut(scion::usize const index_tran) -> std::span<rtk::vecxf> {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            return m.eigenvecs[index_tran];
        }

        /// \brief returns the eval at index index_proj of the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \param index_proj index of the eval to query
        /// \return eval at index index_proj of the tran at index index_tran
        /// \throws std::logic_error if no such tran or eval exist
        [[nodiscard]] inline auto eval(scion::usize const index_tran, scion::usize const index_proj) const -> scion::f32 const & {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            if (this->width_proj() <= index_proj) {
                throw std::logic_error("no such eval");
            }

            return m.eigenvals[index_tran][index_proj];
        }

        /// \brief returns the eval at index index_proj of the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \param index_proj index of the eval to query
        /// \return eval at index index_proj of the tran at index index_tran
        /// \throws std::logic_error if no such tran or eval exist
        /// \note client code is responsible to maintain the invariants of the object if changes are made through the returned reference
        [[nodiscard]] inline auto eval_mut(scion::usize const index_tran, scion::usize const index_proj) -> scion::f32 & {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            if (this->width_proj() <= index_proj) {
                throw std::logic_error("no such eval");
            }

            return m.eigenvals[index_tran][index_proj];
        }

        /// \brief returns the evec at index index_proj of the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \param index_proj index of the evec to query
        /// \return evec at index index_proj of the tran at index index_tran
        /// \throws std::logic_error if no such tran or evec exist
        [[nodiscard]] inline auto evec(scion::usize const index_tran, scion::usize const index_proj) const -> rtk::vecxf const & {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            if (this->width_proj() <= index_proj) {
                throw std::logic_error("no such evec");
            }

            return m.eigenvecs[index_tran][index_proj];
        }

        /// \brief returns the evec at index index_proj of the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \param index_proj index of the evec to query
        /// \return evec at index index_proj of the tran at index index_tran
        /// \throws std::logic_error if no such tran or evec exist
        /// \note client code is responsible to maintain the invariants of the object if changes are made through the returned reference
        [[nodiscard]] inline auto evec_mut(scion::usize const index_tran, scion::usize const index_proj) -> rtk::vecxf & {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            if (this->width_proj() <= index_proj) {
                throw std::logic_error("no such evec");
            }

            return m.eigenvecs[index_tran][index_proj];
        }

        /// \brief returns the targeted iteration level
        /// \return targeted iteration level
        [[nodiscard]] inline auto iterations() const noexcept -> scion::usize {
            return m.iterations;
        }

        /// \brief sets the targeted iteration level
        /// \param iterations new targeted iteration level
        /// \throws std::logic_error if targeted iteration level is greater than 16
        inline void set_iterations(scion::usize const iterations) {
            if (iterations > 16) {
                throw std::logic_error("too many iterations");
            }

            m.iterations = iterations;

            this->request_update_attractor_bn();
        }

        /// \brief returns the attractor in barycentric space
        /// \return attractor in barycentric space
        [[nodiscard]] inline auto attractor_bn() const -> std::span<vecxf const> {
            return m.attractor_bn;
        }

        /// \brief returns the attractor in modeling space using the specified dimension
        /// \return attractor in modeling space using the specified dimension
        /// \throws std::bad_variant_access if it was not computed for the specified dimension
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        [[nodiscard]] inline auto attractor_rn() const -> std::span<vecf<Dim> const> {
            return std::get<std::vector<rtk::vecf<Dim>>>(m.attractor_rn);
        }

        /// \brief sorts the evals and evecs of the tran at index index_tran
        /// \param index_tran index of the tran to query
        /// \throws std::logic_error if no such tran or evec exist
        void sort_eigen(scion::usize const index_tran) {
            if (this->width_tran() <= index_tran) {
                throw std::logic_error("no such tran");
            }

            std::vector<scion::f32> evals = m.eigenvals[index_tran];
            std::vector<rtk::vecxf> evecs = m.eigenvecs[index_tran];

            std::vector<scion::usize> indices{};
            for (scion::usize index = 0; index < this->width_proj(); ++index) indices.emplace_back(index);

            std::ranges::sort(indices, std::greater{}, [&evals](scion::usize const index) -> scion::f32 {
                return std::abs(evals[index]);
            });

            for (scion::usize index = 0; index < this->width_proj(); ++index) {
                m.eigenvals[index_tran][index] = evals[indices[index]];
                m.eigenvecs[index_tran][index] = evecs[indices[index]];
            }
        }

        /// \brief returns the dimension of the attractor in modeling space, current or targeted
        /// \return dimension of the attractor in modeling space, current or targeted
        [[nodiscard]] inline auto dimension() const -> scion::usize {
            if (std::holds_alternative<scion::usize>(m.attractor_rn)) {
                return std::get<scion::usize>(m.attractor_rn);
            }

            if (std::holds_alternative<std::vector<rtk::vecf<0>>>(m.attractor_rn)) {
                return 0;
            }

            if (std::holds_alternative<std::vector<rtk::vecf<2>>>(m.attractor_rn)) {
                return 2;
            }

            if (std::holds_alternative<std::vector<rtk::vecf<3>>>(m.attractor_rn)) {
                return 3;
            }

            throw std::logic_error("dimension not supported");
        }

        /// \brief sets the primitive
        /// \param primitive primitive
        void set_primitive(std::vector<rtk::vecxf> && primitive) {
            m.primitive = std::move(primitive);
        }

        /// \brief sets the dimension of the modeling space in which to compute the attractor
        /// \tparam Dim dimension of the modeling space in which to compute the attractor
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        inline void set_dimension() noexcept {
            if (!std::holds_alternative<std::vector<rtk::vecf<Dim>>>(m.attractor_rn)) {
                m.attractor_rn = Dim;
            }

            this->request_update_attractor_rn();
        }

        /// \brief returns true if the PIFS is in a valid state
        /// \return true if the PIFS is in a valid state
        /// \note the returned value is only relevant after a call to this->update() and before any mutation of its state
        [[nodiscard]] inline auto valid() const noexcept -> bool {
            return m.state;
        }

        /// \brief marks the state as invalid
        inline void invalidate() noexcept {
            rtk::log::warn("validate", "pifs:{} is invalidated", static_cast<void *>(this));

            m.state = false;
        };

        /// \brief validates the eval at index index_proj of the tran at index index_tran
        /// \tparam Log one of the functors in rtk::log to log the error, if any, or rtk::log::silent_t for no log
        /// \param index_tran index of the tran to validate
        /// \param index_proj index of the eval to validate
        /// \return true if the eval is valid
        template <rtk::log::functor Log = rtk::log::silent_t>
        [[nodiscard]] inline auto validate_eval(scion::usize const index_tran, scion::usize const index_proj) const -> bool {
            scion::f32 const eval = this->eval(index_tran, index_proj);

            if (!rtk::is_val_valid(eval)) {
                if constexpr (!std::is_same_v<Log, rtk::log::silent_t>) {
                    Log{}("validate", "eval {}:[{}, {}] is not valid", static_cast<void const *>(this), index_tran, index_proj);
                }

                return false;
            }

            if (!rtk::is_val_eigenval_contraction(eval)) {
                if constexpr (!std::is_same_v<Log, rtk::log::silent_t>) {
                    Log{}("validate", "eval {}:[{}, {}] is not a contraction", static_cast<void const *>(this), index_tran, index_proj);
                }

                return false;
            }

            return true;
        }

        /// \brief validates the evals of the tran at index index_tran
        /// \tparam Log one of the functors in rtk::log to log the error, if any, or rtk::log::silent_t for no log
        /// \param index_tran index of the tran to validate
        /// \return true if the evals is valid
        template <rtk::log::functor Log = rtk::log::silent_t>
        [[nodiscard]] inline auto validate_evals(scion::usize const index_tran) const -> bool {
            std::span<scion::f32 const> const evals = this->evals(index_tran);

            if (!std::ranges::all_of(std::views::iota(0UZ, this->width_proj()), [this, index_tran](scion::usize const index_proj) -> bool {
                return this->validate_eval<Log>(index_tran, index_proj);
            })) {
                return false;
            }

            if (evals[0] < 1.0f - rtk::EPS_MP<scion::f32> || 1.0f + rtk::EPS_MP<scion::f32> < evals[0]) {
                if constexpr (!std::is_same_v<Log, rtk::log::silent_t>) {
                    Log{}("validate", "eval {}:[{}, 0] is not 1", static_cast<void const *>(this), index_tran);
                }

                return false;
            }

            if (!std::ranges::is_sorted(evals | std::views::transform([](scion::f32 const a) -> scion::f32 {
                return std::abs(a);
            }), std::ranges::greater{})) {
                if constexpr (!std::is_same_v<Log, rtk::log::silent_t>) {
                    Log{}("validate", "evals {}:[{}] is not sorted", static_cast<void const *>(this), index_tran);
                }

                return false;
            }

            return true;
        }

        /// \brief validates the eigenvals
        /// \tparam Log one of the functors in rtk::log to log the error, if any, or rtk::log::silent_t for no log
        /// \return true if the eigenvals is valid
        template <rtk::log::functor Log = rtk::log::silent_t>
        [[nodiscard]] inline auto validate_eigenvals() const -> bool {
            return std::ranges::all_of(std::views::iota(0UZ, this->width_tran()), [this](scion::usize const index_tran) -> bool {
                return this->validate_evals<Log>(index_tran);
            });
        }

        /// \brief validates the evec at index index_proj of the tran at index index_tran
        /// \tparam Log one of the functors in rtk::log to log the error, if any, or rtk::log::silent_t for no log
        /// \param index_tran index of the tran to validate
        /// \param index_proj index of the evec to validate
        /// \return true if the evec is valid
        template <rtk::log::functor Log = rtk::log::silent_t>
        [[nodiscard]] inline auto validate_evec(scion::usize const index_tran, scion::usize const index_proj) const -> bool {
            rtk::vecxf const & evec = this->evec(index_tran, index_proj);

            if (!rtk::is_vec_valid(evec)) {
                if constexpr (!std::is_same_v<Log, rtk::log::silent_t>) {
                    Log{}("validate", "evec {}:[{}, {}] is not valid", static_cast<void const *>(this), index_tran, index_proj);
                }

                return false;
            }

            return true;
        }

        /// \brief validates the evecs of the tran at index index_tran
        /// \tparam Log one of the functors in rtk::log to log the error, if any, or rtk::log::silent_t for no log
        /// \param index_tran index of the tran to validate
        /// \return true if the evecs is valid
        template <rtk::log::functor Log = rtk::log::silent_t>
        [[nodiscard]] inline auto validate_evecs(scion::usize const index_tran) const -> bool {
            std::span<rtk::vecxf const> const evecs = this->evecs(index_tran);

            if (!std::ranges::all_of(std::views::iota(0UZ, this->width_proj()), [this, index_tran](scion::usize const index_proj) -> bool {
                return this->validate_evec<Log>(index_tran, index_proj);
            })) {
                return false;
            }

            {
                if (!rtk::is_vec_pnt(evecs[0])) {
                    if constexpr (!std::is_same_v<Log, rtk::log::silent_t>) {
                        Log{}("validate", "evec {}:[{}, 0] does not describe a point", static_cast<void const *>(this), index_tran);
                    }

                    return false;
                }
            }

            {
                for (scion::usize index_proj = 1; index_proj < this->width_proj(); ++index_proj) {
                    if (!rtk::is_vec_vec(evecs[index_proj])) {
                        if constexpr (!std::is_same_v<Log, rtk::log::silent_t>) {
                            Log{}("validate", "evec {}:[{}, {}] does not describe a vector", static_cast<void const *>(this), index_tran, index_proj);
                        }

                        return false;
                    }
                }
            }

            {
                rtk::matxf V = rtk::matxf::Zero(static_cast<scion::isize>(evecs.size()), static_cast<scion::isize>(evecs.size()));

                for (scion::isize d = 0; d < static_cast<scion::isize>(evecs.size()); ++d) {
                    V.col(d) = evecs[d];
                }

                Eigen::FullPivLU<rtk::matxf> lu(V);
                if (lu.rank() != V.rows()) {
                    if constexpr (!std::is_same_v<Log, rtk::log::silent_t>) {
                        Log{}("validate", "evecs {}:[{}] does not construct an invertible matrix", static_cast<void const *>(this), index_tran);
                    }

                    return false;
                }
            }

            return true;
        }

        /// \brief validates the eigenvecs
        /// \tparam Log one of the functors in rtk::log to log the error, if any, or rtk::log::silent_t for no log
        /// \return true if the eigenvecs is valid
        template <rtk::log::functor Log = rtk::log::silent_t>
        [[nodiscard]] inline auto validate_eigenvecs() const -> bool {
            return std::ranges::all_of(std::views::iota(0UZ, this->width_tran()), [this](scion::usize const index_tran) -> bool {
                return this->validate_evecs<Log>(index_tran);
            });
        }

        /// \brief validates the tran at index index_tran
        /// \tparam Log one of the functors in rtk::log to log the error, if any, or rtk::log::silent_t for no log
        /// \param index_tran index of the tran to validate
        /// \return true if the tran is valid
        template <rtk::log::functor Log = rtk::log::silent_t>
        [[nodiscard]] inline auto validate_tran(scion::usize const index_tran) const -> bool {
            rtk::matxf const & tran = this->tran(index_tran);

            if (!rtk::is_mat_valid(tran)) {
                if constexpr (!std::is_same_v<Log, rtk::log::silent_t>) {
                    Log{}("validate", "tran {}:[{}] has an invalid value", static_cast<void const *>(this), index_tran);
                }

                return false;
            }

            if (!rtk::is_mat_semi_right_stochastic(tran)) {
                if constexpr (!std::is_same_v<Log, rtk::log::silent_t>) {
                    Log{}("validate", "tran {}:[{}] is not semi-right-stochastic", static_cast<void const *>(this), index_tran);
                }

                return false;
            }

            return true;
        }

        /// \brief validates the trans
        /// \tparam Log one of the functors in rtk::log to log the error, if any, or rtk::log::silent_t for no log
        /// \return true if the trans is valid
        template <rtk::log::functor Log = rtk::log::silent_t>
        [[nodiscard]] inline auto validate_trans() const -> bool {
            return std::ranges::all_of(std::views::iota(0UZ, this->width_tran()), [this](scion::usize const index_tran) -> bool {
                return this->validate_tran<Log>(index_tran);
            });
        }

        /// \brief validates the state of the pifs and returns true if it is valid
        /// \return true if the pifs was validated
        /// \note this function starts by invalidating the PIFS, then proceed to do the validation
        template <rtk::log::functor Log = rtk::log::silent_t>
        inline auto validate() -> bool {
            bool const state = this->m.state;

            this->m.state = false;

            if (!this->validate_eigenvals<Log>()) {
                if (state) {
                    rtk::log::warn("validate", "pifs:{} is invalidated", static_cast<void *>(this));
                }

                return false;
            }

            if (!this->validate_eigenvecs<Log>()) {
                if (state) {
                    rtk::log::warn("validate", "pifs:{} is invalidated", static_cast<void *>(this));
                }

                return false;
            }

            if (!this->validate_trans<Log>()) {
                if (state) {
                    rtk::log::warn("validate", "pifs:{} is invalidated", static_cast<void *>(this));
                }

                return false;
            }

            rtk::log::trace("validate", "pifs:{} was validated", static_cast<void *>(this));

            this->m.state = true;

            return this->valid();
        };
    };
}
