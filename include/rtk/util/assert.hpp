#pragma once

#include "scion/core/core.hpp"
#include "scion/util/util.hpp"

#include <concepts>

/// \file rtk/util/assert.hpp
/// \brief Provides debug assertion tooling.

namespace rtk
{
    /// \brief Default epsilon used for low precision epsilon-comparisons.
    /// \note This is only the default, every epsilon-comparing function should accept its epsilon as a template argument with this as default.
    /// \note If you need a better epsilon, consider using a dynamic epsilon.
    /// \note If you need perfect comparisons, use fixed point instead of floating point.
    /// \note Some functions might defer to another library's epsilon comparison (i.e. Eigen's).
    template <std::floating_point Sca = scion::f32>
    inline constexpr Sca EPS_LP = static_cast<Sca>(0.1);

    /// \brief Default epsilon used for medium precision epsilon-comparisons.
    /// \note This is only the default, every epsilon-comparing function should accept its epsilon as a template argument with this as default.
    /// \note If you need a better epsilon, consider using a dynamic epsilon.
    /// \note If you need perfect comparisons, use fixed point instead of floating point.
    /// \note Some functions might defer to another library's epsilon comparison (i.e. Eigen's).
    template <std::floating_point Sca = scion::f32>
    inline constexpr Sca EPS_MP = static_cast<Sca>(0.001);

    /// \brief Default epsilon used for high precision epsilon-comparisons.
    /// \note This is only the default, every epsilon-comparing function should accept its epsilon as a template argument with this as default.
    /// \note If you need a better epsilon, consider using a dynamic epsilon.
    /// \note If you need perfect comparisons, use fixed point instead of floating point.
    /// \note Some functions might defer to another library's epsilon comparison (i.e. Eigen's).
    template <std::floating_point Sca = scion::f32>
    inline constexpr Sca EPS_HP = static_cast<Sca>(0.00001);

    /// \brief Compares two floating point values with an epsilon threshold.
    /// \tparam Sca The floating point type.
    /// \tparam Eps The epsilon threshold.
    /// \tparam Thr The threshold under which comparison uses absolute delta.
    /// \param lhs The left operand.
    /// \param rhs The right operand.
    /// \returns 0.0 if the operands are equal, 1.0 if the left operand is greater, -1.0 if the right operand is greater.
    template <std::floating_point Sca = scion::f32, Sca Eps = EPS_HP<Sca>, Sca Thr = std::numeric_limits<Sca>::min()>
    requires (std::numeric_limits<Sca>::epsilon() <= Eps && Eps < static_cast<Sca>(1.0))
    inline constexpr auto eps_compare(Sca const lhs, Sca const rhs) -> Sca {
        if (lhs == rhs) return static_cast<Sca>(0.0);

        Sca const diff = std::abs(lhs - rhs);
        Sca const norm = std::min(std::abs(lhs + rhs), std::numeric_limits<Sca>::max());

        if (diff < std::max(Thr, Eps * norm)) {
            return static_cast<Sca>(0.0);
        } else {
            return lhs < rhs ? static_cast<Sca>(-1.0) : static_cast<Sca>(1.0);
        }
    }

    /// \brief Checks whether a value is valid.
    /// \param val Value to check.
    /// \return Whether the value is valid.
    /// \note A value is valid if it is neither `NaN` nor `Inf`.
    template <std::floating_point Sca = scion::f32>
    constexpr auto is_val_valid(Sca const & val) -> bool {
        return !std::isnan(val) && !std::isinf(val);
    }

    /// \brief Checks whether a vector is valid.
    /// \param vec Vector to check.
    /// \return Whether the vector is valid.
    /// \note A vector is valid if all of its coefficients are valid values.
    /// \see `rtk::is_val_valid`.
    template <std::floating_point Sca = scion::f32, scion::usize Dim = 0>
    constexpr auto is_vec_valid(rtk::vec<Sca, Dim> const& vec) -> bool {
        for (scion::isize idx = 0; idx < vec.size(); ++idx) if (!rtk::is_val_valid(vec(idx))) return false;

        return true;
    }

    /// \brief Checks whether a matrix is valid.
    /// \param mat Matrix to check.
    /// \return Whether the matrix is valid.
    /// \note A matrix is valid if all of its coefficients are valid values.
    /// \see `rtk::is_val_valid`.
    template <std::floating_point Sca = scion::f32, scion::usize Row = 0, scion::usize Col = 0>
    constexpr auto is_mat_valid(rtk::mat<Sca, Row, Col> const & mat) -> bool {
        for (scion::isize row = 0; row < mat.rows(); ++row) for (scion::isize col = 0; col < mat.cols(); ++col) if (!rtk::is_val_valid(mat(row, col))) return false;

        return true;
    }

    /// \brief Checks whether a value is a valid eigenvalue for a contraction.
    /// \param val Value to check.
    /// \return Whether the value is a valid eigenvalue for a contraction.
    /// \note A value is a valid eigenvalue for a contraction if it is a valid value and is at most 1.0.
    /// \see `rtk::is_val_valid`.
    template <std::floating_point Sca = scion::f32>
    inline constexpr auto is_val_eigenval_contraction(Sca const & v) -> bool {
        return std::abs(v) <= static_cast<Sca>(1.0);
    }

    /// \brief Checks whether a vector is a valid eigenvector.
    /// \param vec Vector to check.
    /// \return Whether the vector is a valid eigenvector.
    /// \note A vector is a valid eigenvector if it is a valid vector and is not the null-vector.
    /// \see `rtk::is_val_valid`.
    template <std::floating_point Sca = scion::f32, scion::usize Dim = 0>
    inline constexpr auto is_vec_eigenvec(vec<Sca, Dim> const & vec) -> bool {
        for (scion::isize index = 0; index < vec.size(); ++index) if (!rtk::is_val_valid(vec(index)) || vec.isZero()) return false;

        return true;
    }

    /// \brief Checks whether a vector describes a vector in barycentric coordinates.
    /// \param vec Vector to check.
    /// \return Whether the vector is describes a vector in barycentric coordinates.
    /// \note A vector describes a vector in barycentric coordinates if it is a valid vector and if its coefficients sum to 0.
    /// \see `rtk::is_vec_valid`.
    template <std::floating_point Sca = scion::f32, scion::usize Dim = 0, Sca Eps = rtk::EPS_LP<Sca>>
    requires (std::numeric_limits<Sca>::epsilon() <= Eps && Eps < static_cast<Sca>(1.0))
    inline constexpr auto is_vec_vec(rtk::vec<Sca, Dim> const & vec) -> bool {
        if (!rtk::is_vec_valid(vec)) return false;

        Sca acc = static_cast<Sca>(0.0);

        for (scion::isize index = 0; index < vec.size(); ++index) acc += vec(index);

        return -Eps <= acc && acc <= Eps;
    }

    /// \brief Checks whether a vector describes a point in barycentric coordinates.
    /// \param vec Vector to check.
    /// \return Whether the vector is describes a point in barycentric coordinates.
    /// \note A vector describes a point in barycentric coordinates if it is a valid vector and if its coefficients sum to 1.
    /// \see `rtk::is_vec_valid`.
    template <std::floating_point Sca = scion::f32, scion::usize Dim = 0, Sca Eps = rtk::EPS_LP<Sca>>
    requires (std::numeric_limits<Sca>::epsilon() <= Eps && Eps < static_cast<Sca>(1.0))
    inline constexpr auto is_vec_pnt(rtk::vec<Sca, Dim> const & vec) -> bool {
        if (!rtk::is_vec_valid(vec)) return false;

        Sca acc = static_cast<Sca>(0.0);

        for (scion::isize index = 0; index < vec.size(); ++index) acc += vec(index);

        return static_cast<Sca>(1.0) - Eps <= std::abs(acc) && std::abs(acc) <= static_cast<Sca>(1.0) + Eps;
    }

    /// \brief Checks whether a matrix is a semi-right-stochastic matrix.
    /// \param mat Matrix to check.
    /// \return Whether a matrix is a semi-right-stochastic matrix.
    /// \note A matrix is a semi-right-stochastic matrix if it is a valid matrix and is a stochastic matrix not restricted to non-negative coefficients.
    /// \see `rtk::is_mat_valid`.
    template <std::floating_point Sca = scion::f32, scion::usize Row = 0, scion::usize Col = 0, Sca Eps = rtk::EPS_LP<Sca>>
    requires (std::numeric_limits<Sca>::epsilon() <= Eps && Eps < static_cast<Sca>(1.0))
    inline constexpr auto is_mat_semi_right_stochastic(mat<Sca, Row, Col> const & mat) -> bool {
        if (!rtk::is_mat_valid(mat)) return false;

        for (scion::isize row = 0; row < mat.rows(); ++row) {
            Sca acc = static_cast<Sca>(0);

            for (scion::isize col = 0; col < mat.cols(); ++col) {
                if (!rtk::is_val_valid(mat(col, row))) return false;

                acc += mat(col, row);
            }

            if (std::abs(acc - static_cast<Sca>(1.0)) > Eps) return false;
        }

        return true;
    }
}
