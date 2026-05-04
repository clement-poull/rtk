#pragma once

#include "scion/core/core.hpp"

#include "rtk/core/types.hpp"

/// \file rtk/core/math.hpp
/// \brief provides some math functions requires by other modules

namespace rtk {
    template<typename S>
    [[nodiscard]] auto vector_product_tensor_flatten(rtk::vecx<S> const& lhs, rtk::vecx<S> const& rhs) -> rtk::vecx<S> {
        rtk::vecx<S> out = rtk::vecx<S>::Zero(lhs.size() * rhs.size());

        for (scion::usize i = 0; i < lhs.size(); ++i) {
            for (scion::usize j = 0; j < rhs.size(); ++j) {
                out(i * rhs.size() + j) = lhs(i) * rhs(j);
            }
        }

        return out;
    }

    template<typename S>
    [[nodiscard]] auto matrix_product_tensor_flatten(rtk::matx<S> const& lhs, rtk::matx<S> const& rhs) -> rtk::matx<S> {
        rtk::matx<S> out = rtk::matx<S>::Zero(lhs.size() * rhs.size(), lhs.size() * rhs.size());

        if (lhs.cols() != lhs.rows()) {
            throw std::logic_error("non square matrices");
        }

        if (rhs.cols() != rhs.rows()) {
            throw std::logic_error("non square matrices");
        }

        if (lhs.cols() != rhs.cols()) {
            throw std::logic_error("non square matrices");
        }

        scion::usize const rows = lhs.rows();
        scion::usize const cols = lhs.cols();

        for (scion::usize col_i = 0; col_i < cols; ++col_i) {
            for (scion::usize row_i = 0; row_i < rows; ++row_i) {
                for (scion::usize col_j = 0; col_j < cols; ++col_j) {
                    for (scion::usize row_j = 0; row_j < rows; ++row_j) {
                        out(col_i * cols + col_j, row_i * rows + row_j) = lhs(col_i, row_i) * rhs(col_j, row_j);
                    }
                }
            }
        }

        return out;
    }
}
