#pragma once

#include "scion/scion.hpp"

#include "rtk/core/core.hpp"
#include "rtk/util/util.hpp"

#include "rtk/constants.hpp"

/// \file rtk/rtk.hpp
/// \brief Main file of the library.

namespace rtk::version {
    inline constexpr scion::usize MAJOR = RTK_VERSION_MAJOR;
    inline constexpr scion::usize MINOR = RTK_VERSION_MINOR;
    inline constexpr scion::usize PATCH = RTK_VERSION_PATCH;
}
