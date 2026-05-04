#pragma once

#include "scion/util/log.hpp"

#include <concepts>
#include <format>
#include <utility>

/// \file rtk/util/log.hpp
/// \brief Provides leveled logging capabilities.

#ifdef RTK_NOLOGS
    #define RTK_NOLOG_ERROR
    #define RTK_NOLOG_INFO
    #define RTK_NOLOG_WARN
    #define RTK_NOLOG_DEBUG
    #define RTK_NOLOG_TRACE
#endif

#ifdef SCION_CONFIGURATION_RELEASE
    #define RTK_NOLOG_DEBUG
    #define RTK_NOLOG_TRACE
#endif

#ifdef SCION_CONFIGURATION_PROFILE
    #define RTK_NOLOG_TRACE
#endif

namespace rtk::log
{
    #ifndef SCION_NOLOG_ERROR
        SCION_IMPL_LOG_FUNCTION_ON(error, scion::log::ANSI_ERROR)
    #else
        SCION_IMPL_LOG_FUNCTION_NO(error, scion::log::ANSI_ERROR)
    #endif

    #ifndef SCION_NOLOG_WARN
        SCION_IMPL_LOG_FUNCTION_ON(warn, scion::log::ANSI_WARN)
    #else
        SCION_IMPL_LOG_FUNCTION_NO(warn, scion::log::ANSI_WARN)
    #endif

    #ifndef SCION_NOLOG_INFO
        SCION_IMPL_LOG_FUNCTION_ON(info, scion::log::ANSI_INFO)
    #else
        SCION_IMPL_LOG_FUNCTION_NO(info, scion::log::ANSI_INFO)
    #endif

    #ifndef SCION_NOLOG_DEBUG
        SCION_IMPL_LOG_FUNCTION_ON(debug, scion::log::ANSI_DEBUG)
    #else
        SCION_IMPL_LOG_FUNCTION_NO(debug, scion::log::ANSI_DEBUG)
    #endif

    #ifndef RTK_NOLOG_TRACE
        SCION_IMPL_LOG_FUNCTION_ON(trace, scion::log::ANSI_TRACE)
    #else
        SCION_IMPL_LOG_FUNCTION_NO(trace, scion::log::ANSI_TRACE)
    #endif

    struct error_t {
        template <typename... Args>
        inline void operator()(char const * prefix, std::format_string<Args...> fmt, Args &&... args) const {
            rtk::log::error(prefix, fmt, std::forward<Args>(args)...);
        }
    };

    struct warn_t {
        template <typename... Args>
        inline void operator()(char const * prefix, std::format_string<Args...> fmt, Args &&... args) const {
            rtk::log::warn(prefix, fmt, std::forward<Args>(args)...);
        }
    };

    struct info_t {
        template <typename... Args>
        inline void operator()(char const * prefix, std::format_string<Args...> fmt, Args &&... args) const {
            rtk::log::info(prefix, fmt, std::forward<Args>(args)...);
        }
    };

    struct debug_t {
        template <typename... Args>
        inline void operator()(char const * prefix, std::format_string<Args...> fmt, Args &&... args) const {
            rtk::log::debug(prefix, fmt, std::forward<Args>(args)...);
        }
    };

    struct trace_t {
        template <typename... Args>
        inline void operator()(char const * prefix, std::format_string<Args...> fmt, Args &&... args) const {
            rtk::log::trace(prefix, fmt, std::forward<Args>(args)...);
        }
    };

    struct silent_t {
        template <typename... Args>
        inline void operator()(char const *, std::format_string<Args...>, Args&&...) const {
            // DO NOTHING
        }
    };

    template <typename T>
    concept functor = std::same_as<T, error_t> || std::same_as<T, warn_t> || std::same_as<T, info_t> || std::same_as<T, debug_t> || std::same_as<T, trace_t> || std::same_as<T, silent_t>;
}
