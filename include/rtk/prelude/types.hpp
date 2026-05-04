#pragma once

#include "scion/core/core.hpp"

#include "rtk/core/types.hpp"

/// \file rtk/prelude/types.hpp
/// \brief Re-exports the types in `rtk/core/types.hpp` to the root namespace.
/// \note Do not import this file in library header, only use it in source files or client code.

template<typename Sca, scion::usize Dim>
using vec = rtk::vec<Sca, Dim>;

template<scion::usize Dim>
using vecf = rtk::vecf<Dim>;
template<scion::usize Dim>
using vecd = rtk::vecd<Dim>;
template<scion::usize Dim>
using vecu = rtk::vecu<Dim>;
template<scion::usize Dim>
using veci = rtk::veci<Dim>;

using vecxf = rtk::vecxf;
using vecxd = rtk::vecxd;
using vecxu = rtk::vecxu;
using vecxi = rtk::vecxi;

using vec2f = rtk::vec2f;
using vec2d = rtk::vec2d;
using vec2u = rtk::vec2u;
using vec2i = rtk::vec2i;

using vec3f = rtk::vec3f;
using vec3d = rtk::vec3d;
using vec3u = rtk::vec3u;
using vec3i = rtk::vec3i;

using vec4f = rtk::vec4f;
using vec4d = rtk::vec4d;
using vec4u = rtk::vec4u;
using vec4i = rtk::vec4i;

template<typename Sca, scion::usize Row, scion::usize Col>
using mat = rtk::mat<Sca, Row, Col>;

template<scion::usize Row, scion::usize Col>
using matf = rtk::matf<Row, Col>;
template<scion::usize Row, scion::usize Col>
using matd = rtk::matd<Row, Col>;
template<scion::usize Row, scion::usize Col>
using matu = rtk::matu<Row, Col>;
template<scion::usize Row, scion::usize Col>
using mati = rtk::mati<Row, Col>;

using matxf = rtk::matxf;
using matxd = rtk::matxd;
using matxu = rtk::matxu;
using matxi = rtk::matxi;

using mat2x2f = rtk::mat2x2f;
using mat2x2d = rtk::mat2x2d;
using mat2x2u = rtk::mat2x2u;
using mat2x2i = rtk::mat2x2i;

using mat2x3f = rtk::mat2x3f;
using mat2x3d = rtk::mat2x3d;
using mat2x3u = rtk::mat2x3u;
using mat2x3i = rtk::mat2x3i;

using mat2x4f = rtk::mat2x4f;
using mat2x4d = rtk::mat2x4d;
using mat2x4u = rtk::mat2x4u;
using mat2x4i = rtk::mat2x4i;

using mat3x2f = rtk::mat3x2f;
using mat3x2d = rtk::mat3x2d;
using mat3x2u = rtk::mat3x2u;
using mat3x2i = rtk::mat3x2i;

using mat3x3f = rtk::mat3x3f;
using mat3x3d = rtk::mat3x3d;
using mat3x3u = rtk::mat3x3u;
using mat3x3i = rtk::mat3x3i;

using mat3x4f = rtk::mat3x4f;
using mat3x4d = rtk::mat3x4d;
using mat3x4u = rtk::mat3x4u;
using mat3x4i = rtk::mat3x4i;

using mat4x2f = rtk::mat4x2f;
using mat4x2d = rtk::mat4x2d;
using mat4x2u = rtk::mat4x2u;
using mat4x2i = rtk::mat4x2i;

using mat4x3f = rtk::mat4x3f;
using mat4x3d = rtk::mat4x3d;
using mat4x3u = rtk::mat4x3u;
using mat4x3i = rtk::mat4x3i;

using mat4x4f = rtk::mat4x4f;
using mat4x4d = rtk::mat4x4d;
using mat4x4u = rtk::mat4x4u;
using mat4x4i = rtk::mat4x4i;

using mat2f = rtk::mat2f;
using mat2d = rtk::mat2d;
using mat2u = rtk::mat2u;
using mat2i = rtk::mat2i;

using mat3f = rtk::mat3f;
using mat3d = rtk::mat3d;
using mat3u = rtk::mat3u;
using mat3i = rtk::mat3i;

using mat4f = rtk::mat4f;
using mat4d = rtk::mat4d;
using mat4u = rtk::mat4u;
using mat4i = rtk::mat4i;
