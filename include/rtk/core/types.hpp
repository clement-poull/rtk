#pragma once

#include "scion/core/core.hpp"

#include "Eigen/Dense"

#include <type_traits>

/// \file rtk/core/types.hpp
/// \brief Provides coherent types.

namespace rtk {
    /// \brief Base vector type (uses Eigen's Vector type).
    ///
    /// Types commonly used in computer graphics have specific alias.
    /// Aliases are also provided for dynamically-sized vectors with the aforementioned types.
    template<typename Sca, scion::usize Dim>
    using vec = std::conditional_t<Dim == 0, Eigen::VectorX<Sca>, Eigen::Vector<Sca, Dim>>;

    /// \brief Base vec type with dynamic coefficients.
    template<typename Sca>
    using vecx = rtk::vec<Sca, 0>;

    /// \brief Base vec with IEEE754 simple precision coefficients.
    template<scion::usize Dim>
    using vecf = rtk::vec<scion::f32, Dim>;
    /// \brief Base vec with IEEE754 double precision coefficients.
    template<scion::usize Dim>
    using vecd = rtk::vec<scion::f64, Dim>;
    /// \brief Base vec with unsigned integer coefficients.
    template<scion::usize Dim>
    using vecu = rtk::vec<scion::i32, Dim>;
    /// \brief Base vec with signed integer coefficients.
    template<scion::usize Dim>
    using veci = rtk::vec<scion::u32, Dim>;

    /// \brief Vec type with dynamic IEEE754 simple precision coefficients.
    using vecxf = rtk::vecf<0>;
    /// \brief Vec with dynamic IEEE754 double precision coefficients.
    using vecxd = rtk::vecd<0>;
    /// \brief Vec with dynamic unsigned integer coefficients.
    using vecxu = rtk::vecu<0>;
    /// \brief Vec with dynamic signed integer coefficients.
    using vecxi = rtk::veci<0>;

    /// \brief Vec type with 2 IEEE754 simple precision coefficients.
    using vec2f = rtk::vecf<2>;
    /// \brief Vec type with 2 IEEE754 double precision coefficients.
    using vec2d = rtk::vecd<2>;
    /// \brief Vec type with 2 unsigned integer coefficients.
    using vec2u = rtk::vecu<2>;
    /// \brief Vec type with 2 signed integer coefficients.
    using vec2i = rtk::veci<2>;

    /// \brief Vec type with 3 IEEE754 simple precision coefficients.
    using vec3f = rtk::vecf<3>;
    /// \brief Vec type with 3 IEEE754 double precision coefficients.
    using vec3d = rtk::vecd<3>;
    /// \brief Vec type with 3 unsigned integer coefficients.
    using vec3u = rtk::vecu<3>;
    /// \brief Vec type with 3 signed integer coefficients.
    using vec3i = rtk::veci<3>;

    /// \brief Vec type with 4 IEEE754 simple precision coefficients.
    using vec4f = rtk::vecf<4>;
    /// \brief Vec type with 4 IEEE754 double precision coefficients.
    using vec4d = rtk::vecd<4>;
    /// \brief Vec type with 4 unsigned integer coefficients.
    using vec4u = rtk::vecu<4>;
    /// \brief Vec type with 4 signed integer coefficients.
    using vec4i = rtk::veci<4>;

    /// \brief Base matrix type (uses Eigen's Matrix type).
    ///
    /// Types commonly used in computer graphics have specific alias:
    /// Sca = f32, f64, u32, i32
    /// Dim = 2x2, 2x3, 2x4, 3x2, 3x3, 3x4, 4x2, 4x3, 4x4
    ///
    /// Aliases are also provided for dynamically-sized mat and dynamically-sized square mat with the aforementioned types.
    template<typename Sca, scion::usize Row, scion::usize Col>
    using mat = std::conditional_t<Row == 0 || Col == 0, Eigen::MatrixX<Sca>, Eigen::Matrix<Sca, Row, Col>>;

    /// \brief Base matrix type with dynamic coefficients.
    template<typename Sca>
    using matx = rtk::mat<Sca, 0, 0>;

    /// \brief Base matrix type with IEEE754 simple precision coefficients.
    template<scion::usize Row, scion::usize Col>
    using matf = rtk::mat<scion::f32, Row, Col>;
    /// \brief Base matrix type with IEEE754 double precision coefficients.
    template<scion::usize Row, scion::usize Col>
    using matd = rtk::mat<scion::f64, Row, Col>;
    /// \brief Base matrix type with unsigned integer coefficients.
    template<scion::usize Row, scion::usize Col>
    using matu = rtk::mat<scion::u32, Row, Col>;
    /// \brief Base matrix type with signed integer coefficients.
    template<scion::usize Row, scion::usize Col>
    using mati = rtk::mat<scion::i32, Row, Col>;

    /// \brief Matrix type with dynamic IEEE754 simple precision coefficients.
    using matxf = rtk::mat<scion::f32, 0, 0>;
    /// \brief Matrix type with dynamic IEEE754 double precision coefficients.
    using matxd = rtk::mat<scion::f64, 0, 0>;
    /// \brief Matrix type with dynamic unsigned integer coefficients.
    using matxu = rtk::mat<scion::u32, 0, 0>;
    /// \brief Matrix type with dynamic signed integer coefficients.
    using matxi = rtk::mat<scion::i32, 0, 0>;

    /// \brief Matrix type with 2x2 IEEE754 simple precision coefficients.
    using mat2x2f = rtk::matf<2, 2>;
    /// \brief Matrix type with 2x2 IEEE754 double precision coefficients.
    using mat2x2d = rtk::matd<2, 2>;
    /// \brief Matrix type with 2x2 unsigned integer coefficients.
    using mat2x2u = rtk::matu<2, 2>;
    /// \brief Matrix type with 2x2 signed integer coefficients.
    using mat2x2i = rtk::mati<2, 2>;

    /// \brief Matrix type with 2x3 IEEE754 simple precision coefficients.
    using mat2x3f = rtk::matf<2, 3>;
    /// \brief Matrix type with 2x3 IEEE754 double precision coefficients.
    using mat2x3d = rtk::matd<2, 3>;
    /// \brief Matrix type with 2x3 unsigned integer coefficients.
    using mat2x3u = rtk::matu<2, 3>;
    /// \brief Matrix type with 2x3 signed integer coefficients.
    using mat2x3i = rtk::mati<2, 3>;

    /// \brief Matrix type with 2x4 IEEE754 simple precision coefficients.
    using mat2x4f = rtk::matf<2, 4>;
    /// \brief Matrix type with 2x4 IEEE754 double precision coefficients.
    using mat2x4d = rtk::matd<2, 4>;
    /// \brief Matrix type with 2x4 unsigned integer coefficients.
    using mat2x4u = rtk::matu<2, 4>;
    /// \brief Matrix type with 2x4 signed integer coefficients.
    using mat2x4i = rtk::mati<2, 4>;

    /// \brief Matrix type with 3x2 IEEE754 simple precision coefficients.
    using mat3x2f = rtk::matf<3, 2>;
    /// \brief Matrix type with 3x2 IEEE754 double precision coefficients.
    using mat3x2d = rtk::matd<3, 2>;
    /// \brief Matrix type with 3x2 unsigned integer coefficients.
    using mat3x2u = rtk::matu<3, 2>;
    /// \brief Matrix type with 3x2 signed integer coefficients.
    using mat3x2i = rtk::mati<3, 2>;

    /// \brief Matrix type with 3x3 IEEE754 simple precision coefficients.
    using mat3x3f = rtk::matf<3, 3>;
    /// \brief Matrix type with 3x3 IEEE754 double precision coefficients.
    using mat3x3d = rtk::matd<3, 3>;
    /// \brief Matrix type with 3x3 unsigned integer coefficients.
    using mat3x3u = rtk::matu<3, 3>;
    /// \brief Matrix type with 3x3 signed integer coefficients.
    using mat3x3i = rtk::mati<3, 3>;

    /// \brief Matrix type with 3x4 IEEE754 simple precision coefficients.
    using mat3x4f = rtk::matf<3, 4>;
    /// \brief Matrix type with 3x4 IEEE754 double precision coefficients.
    using mat3x4d = rtk::matd<3, 4>;
    /// \brief Matrix type with 3x4 unsigned integer coefficients.
    using mat3x4u = rtk::matu<3, 4>;
    /// \brief Matrix type with 3x4 signed integer coefficients.
    using mat3x4i = rtk::mati<3, 4>;

    /// \brief Matrix type with 4x2 IEEE754 simple precision coefficients.
    using mat4x2f = rtk::matf<4, 2>;
    /// \brief Matrix type with 4x2 IEEE754 double precision coefficients.
    using mat4x2d = rtk::matd<4, 2>;
    /// \brief Matrix type with 4x2 unsigned integer coefficients.
    using mat4x2u = rtk::matu<4, 2>;
    /// \brief Matrix type with 4x2 signed integer coefficients.
    using mat4x2i = rtk::mati<4, 2>;

    /// \brief Matrix type with 4x3 IEEE754 simple precision coefficients.
    using mat4x3f = rtk::matf<4, 3>;
    /// \brief Matrix type with 4x3 IEEE754 double precision coefficients.
    using mat4x3d = rtk::matd<4, 3>;
    /// \brief Matrix type with 4x3 unsigned integer coefficients.
    using mat4x3u = rtk::matu<4, 3>;
    /// \brief Matrix type with 4x3 signed integer coefficients.
    using mat4x3i = rtk::mati<4, 3>;

    /// \brief Matrix type with 4x4 IEEE754 simple precision coefficients.
    using mat4x4f = rtk::matf<4, 4>;
    /// \brief Matrix type with 4x4 IEEE754 double precision coefficients.
    using mat4x4d = rtk::matd<4, 4>;
    /// \brief Matrix type with 4x4 unsigned integer coefficients.
    using mat4x4u = rtk::matu<4, 4>;
    /// \brief Matrix type with 4x4 signed integer coefficients.
    using mat4x4i = rtk::mati<4, 4>;

    /// \brief Square matrix type with 2x2 IEEE754 simple precision coefficients.
    using mat2f = rtk::mat2x2f;
    /// \brief Square matrix type with 2x2 IEEE754 double precision coefficients.
    using mat2d = rtk::mat2x2d;
    /// \brief Square matrix type with 2x2 unsigned integer coefficients.
    using mat2u = rtk::mat2x2u;
    /// \brief Square matrix type with 2x2 signed integer coefficients.
    using mat2i = rtk::mat2x2i;

    /// \brief Square matrix type with 3x3 IEEE754 simple precision coefficients.
    using mat3f = rtk::mat3x3f;
    /// \brief Square matrix type with 3x3 IEEE754 double precision coefficients.
    using mat3d = rtk::mat3x3d;
    /// \brief Square matrix type with 3x3 unsigned integer coefficients.
    using mat3u = rtk::mat3x3u;
    /// \brief Square matrix type with 3x3 signed integer coefficients.
    using mat3i = rtk::mat3x3i;

    /// \brief Square matrix type with 4x4 IEEE754 simple precision coefficients.
    using mat4f = rtk::mat4x4f;
    /// \brief Square matrix type with 4x4 IEEE754 double precision coefficients.
    using mat4d = rtk::mat4x4d;
    /// \brief Square matrix type with 4x4 unsigned integer coefficients.
    using mat4u = rtk::mat4x4u;
    /// \brief Square matrix type with 4x4 signed integer coefficients.
    using mat4i = rtk::mat4x4i;
}

static_assert(sizeof(rtk::vec2f) == sizeof(scion::f32) * 2, "vec2f is not contiguous");
static_assert(sizeof(rtk::vec2d) == sizeof(scion::f64) * 2, "vec2d is not contiguous");
static_assert(sizeof(rtk::vec2i) == sizeof(scion::i32) * 2, "vec2i is not contiguous");
static_assert(sizeof(rtk::vec2u) == sizeof(scion::u32) * 2, "vec2u is not contiguous");
static_assert(sizeof(rtk::vec3f) == sizeof(scion::f32) * 3, "vec3f is not contiguous");
static_assert(sizeof(rtk::vec3d) == sizeof(scion::f64) * 3, "vec3d is not contiguous");
static_assert(sizeof(rtk::vec3i) == sizeof(scion::i32) * 3, "vec3i is not contiguous");
static_assert(sizeof(rtk::vec3u) == sizeof(scion::u32) * 3, "vec3u is not contiguous");
static_assert(sizeof(rtk::vec4f) == sizeof(scion::f32) * 4, "vec4f is not contiguous");
static_assert(sizeof(rtk::vec4d) == sizeof(scion::f64) * 4, "vec4d is not contiguous");
static_assert(sizeof(rtk::vec4i) == sizeof(scion::i32) * 4, "vec4i is not contiguous");
static_assert(sizeof(rtk::vec4u) == sizeof(scion::u32) * 4, "vec4u is not contiguous");
static_assert(sizeof(rtk::mat2x2f) == sizeof(scion::f32) * 2 * 2, "mat2x2f is not contiguous");
static_assert(sizeof(rtk::mat2x2d) == sizeof(scion::f64) * 2 * 2, "mat2x2d is not contiguous");
static_assert(sizeof(rtk::mat2x2i) == sizeof(scion::i32) * 2 * 2, "mat2x2i is not contiguous");
static_assert(sizeof(rtk::mat2x2u) == sizeof(scion::u32) * 2 * 2, "mat2x2u is not contiguous");
static_assert(sizeof(rtk::mat3x2f) == sizeof(scion::f32) * 3 * 2, "mat3x2f is not contiguous");
static_assert(sizeof(rtk::mat3x2d) == sizeof(scion::f64) * 3 * 2, "mat3x2d is not contiguous");
static_assert(sizeof(rtk::mat3x2i) == sizeof(scion::i32) * 3 * 2, "mat3x2i is not contiguous");
static_assert(sizeof(rtk::mat3x2u) == sizeof(scion::u32) * 3 * 2, "mat3x2u is not contiguous");
static_assert(sizeof(rtk::mat4x2f) == sizeof(scion::f32) * 4 * 2, "mat4x2f is not contiguous");
static_assert(sizeof(rtk::mat4x2d) == sizeof(scion::f64) * 4 * 2, "mat4x2d is not contiguous");
static_assert(sizeof(rtk::mat4x2i) == sizeof(scion::i32) * 4 * 2, "mat4x2i is not contiguous");
static_assert(sizeof(rtk::mat4x2u) == sizeof(scion::u32) * 4 * 2, "mat4x2u is not contiguous");
static_assert(sizeof(rtk::mat2x3f) == sizeof(scion::f32) * 2 * 3, "mat2x3f is not contiguous");
static_assert(sizeof(rtk::mat2x3d) == sizeof(scion::f64) * 2 * 3, "mat2x3d is not contiguous");
static_assert(sizeof(rtk::mat2x3i) == sizeof(scion::i32) * 2 * 3, "mat2x3i is not contiguous");
static_assert(sizeof(rtk::mat2x3u) == sizeof(scion::u32) * 2 * 3, "mat2x3u is not contiguous");
static_assert(sizeof(rtk::mat3x3f) == sizeof(scion::f32) * 3 * 3, "mat3x3f is not contiguous");
static_assert(sizeof(rtk::mat3x3d) == sizeof(scion::f64) * 3 * 3, "mat3x3d is not contiguous");
static_assert(sizeof(rtk::mat3x3i) == sizeof(scion::i32) * 3 * 3, "mat3x3i is not contiguous");
static_assert(sizeof(rtk::mat3x3u) == sizeof(scion::u32) * 3 * 3, "mat3x3u is not contiguous");
static_assert(sizeof(rtk::mat4x3f) == sizeof(scion::f32) * 4 * 3, "mat4x3f is not contiguous");
static_assert(sizeof(rtk::mat4x3d) == sizeof(scion::f64) * 4 * 3, "mat4x3d is not contiguous");
static_assert(sizeof(rtk::mat4x3i) == sizeof(scion::i32) * 4 * 3, "mat4x3i is not contiguous");
static_assert(sizeof(rtk::mat4x3u) == sizeof(scion::u32) * 4 * 3, "mat4x3u is not contiguous");
static_assert(sizeof(rtk::mat2x4f) == sizeof(scion::f32) * 2 * 4, "mat2x4f is not contiguous");
static_assert(sizeof(rtk::mat2x4d) == sizeof(scion::f64) * 2 * 4, "mat2x4d is not contiguous");
static_assert(sizeof(rtk::mat2x4i) == sizeof(scion::i32) * 2 * 4, "mat2x4i is not contiguous");
static_assert(sizeof(rtk::mat2x4u) == sizeof(scion::u32) * 2 * 4, "mat2x4u is not contiguous");
static_assert(sizeof(rtk::mat3x4f) == sizeof(scion::f32) * 3 * 4, "mat3x4f is not contiguous");
static_assert(sizeof(rtk::mat3x4d) == sizeof(scion::f64) * 3 * 4, "mat3x4d is not contiguous");
static_assert(sizeof(rtk::mat3x4i) == sizeof(scion::i32) * 3 * 4, "mat3x4i is not contiguous");
static_assert(sizeof(rtk::mat3x4u) == sizeof(scion::u32) * 3 * 4, "mat3x4u is not contiguous");
static_assert(sizeof(rtk::mat4x4f) == sizeof(scion::f32) * 4 * 4, "mat4x4f is not contiguous");
static_assert(sizeof(rtk::mat4x4d) == sizeof(scion::f64) * 4 * 4, "mat4x4d is not contiguous");
static_assert(sizeof(rtk::mat4x4i) == sizeof(scion::i32) * 4 * 4, "mat4x4i is not contiguous");
static_assert(sizeof(rtk::mat4x4u) == sizeof(scion::u32) * 4 * 4, "mat4x4u is not contiguous");

static_assert(std::is_standard_layout_v<rtk::vec2f>, "vec2f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec2d>, "vec2d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec2i>, "vec2i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec2u>, "vec2u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec3f>, "vec3f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec3d>, "vec3d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec3i>, "vec3i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec3u>, "vec3u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec4f>, "vec4f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec4d>, "vec4d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec4i>, "vec4i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::vec4u>, "vec4u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x2f>, "mat2x2f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x2d>, "mat2x2d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x2i>, "mat2x2i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x2u>, "mat2x2u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x2f>, "mat3x2f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x2d>, "mat3x2d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x2i>, "mat3x2i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x2u>, "mat3x2u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x2f>, "mat4x2f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x2d>, "mat4x2d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x2i>, "mat4x2i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x2u>, "mat4x2u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x3f>, "mat2x3f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x3d>, "mat2x3d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x3i>, "mat2x3i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x3u>, "mat2x3u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x3f>, "mat3x3f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x3d>, "mat3x3d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x3i>, "mat3x3i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x3u>, "mat3x3u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x3f>, "mat4x3f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x3d>, "mat4x3d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x3i>, "mat4x3i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x3u>, "mat4x3u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x4f>, "mat2x4f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x4d>, "mat2x4d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x4i>, "mat2x4i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat2x4u>, "mat2x4u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x4f>, "mat3x4f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x4d>, "mat3x4d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x4i>, "mat3x4i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat3x4u>, "mat3x4u is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x4f>, "mat4x4f is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x4d>, "mat4x4d is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x4i>, "mat4x4i is not standard layout");
static_assert(std::is_standard_layout_v<rtk::mat4x4u>, "mat4x4u is not standard layout");
