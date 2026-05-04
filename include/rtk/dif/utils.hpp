#pragma once

#include "scion/core/core.hpp"

#include <ranges>

/// \file rtk/dif/utils.hpp
/// \brief provides internal utilities for the dif module

namespace rtk::dif {
    /// \brief computes the coordinates of the starting point in the specified eigensystem
    /// \param width_proj Number of control points.
    /// \param evecs eigenvectors to use
    /// \param starting_point
    /// \return $x_i$.
    template<std::floating_point S>
    auto compute_xis(scion::usize const width_proj, std::span<const rtk::vecx<S>> evecs, rtk::vecx<S> const& starting_point)->std::vector<S> {
        rtk::matx<S> T = rtk::matx<S>::Identity(static_cast<scion::isize>(width_proj), static_cast<scion::isize>(width_proj));

        for (scion::isize v_i = 0; v_i < static_cast<scion::isize>(width_proj); ++v_i) {
            T.col(v_i) = evecs[v_i];
        }

        rtk::vecx<S> p = T.inverse() * starting_point;

        return std::vector<S>{p.data(), p.data() + p.size()};
    }

    /// \brief Computes the $x_i$.
    /// \param width_proj Number of control points.
    /// \param evecs Eigenvectors.
    /// \param starting_point
    /// \return $x_i$.
    template<std::floating_point S>
    auto compute_xis_force_positive(scion::usize const width_proj, std::span<rtk::vecx<S>> evecs, rtk::vecx<S> const& starting_point)->std::vector<S> {
        std::vector<S> xis = rtk::dif::compute_xis<S>(width_proj, evecs, starting_point);

        for (scion::usize index = 0; index < width_proj; ++index) {
            if (xis[index] < 0) {
                xis[index] = -xis[index];
                evecs[index] = -evecs[index];
            }
        }

        return std::move(xis);
    }

    /// \brief Computes $\alpha_i$.
    /// \param width_proj Number of control points.
    /// \param eigenvals Eigenvalues.
    /// \return $\alpha_i$.
    template<std::floating_point S>
    auto compute_ais([[maybe_unused]] scion::usize const width_proj, std::span<const S> eigenvals) -> std::vector<S> {
        S const lambda_1 = std::log(std::abs(eigenvals[1]));

        // WAIT(GCC13 deprecated)
        #ifdef __cpp_lib_ranges_to_container
            return eigenvals | std::views::transform([lambda_1](S const& eigenval) -> S { return std::log(std::abs(eigenval)) / lambda_1; }) | std::ranges::to<std::vector<S>>();
        #else
        {
            std::vector<S> ais;
            ais.reserve(width_proj);
            for (S const& eigenval : eigenvals) ais.push_back(std::log(std::abs(eigenval)) / lambda_1);
            return ais;
        }
        #endif
    }
}
