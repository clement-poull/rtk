#pragma once

#include "scion/core/core.hpp"

#include "rtk/core/types.hpp"

#include <concepts>
#include <utility>

/// \file rtk/core/vertices.hpp
/// \brief provides a struct to represent a vertex with its attributes

namespace rtk
{
    /// \brief a vertex element
    /// \tparam LI local index
    /// \tparam T type of the vertex element
    template <scion::usize LI, typename T>
    class vertex_elem {
        T m_value;

    protected:
        /// \brief constructs a vertex element
        /// \param value value of the vertex element
        explicit constexpr vertex_elem(T const & value) : m_value{value} {
            // do nothing
        }

        /// \brief constructs a vertex element
        /// \param value value of the vertex element
        explicit constexpr vertex_elem(T && value) : m_value{std::move(value)} {
            // do nothing
        }

        /// \brief accessor for the vertex element
        /// \tparam QI queried index
        /// \return value of the vertex element
        template <scion::usize QI>
            requires (QI == LI)
        constexpr auto get() -> T & {
            return this->m_value;
        }

        /// \brief accessor for the vertex element
        /// \tparam QI queried index
        /// \return value of the vertex element
        template <scion::usize QI>
            requires (QI == LI)
        constexpr auto get() const -> T const & {
            return this->m_value;
        }
    };

    template <typename I, typename... T>
    class vertex_impl;

    /// \brief implementation details for a vertex
    /// \tparam I index sequence for the vertex implementation
    /// \tparam T type sequence for the vertex implementation
    template <scion::usize... I, typename... T>
    class vertex_impl<std::index_sequence<I...>, T...> : rtk::vertex_elem<I, T>... {
    public:
        /// \brief constructs all the vertex's elements
        /// \param args
        explicit constexpr vertex_impl(T const &... args) : rtk::vertex_elem<I, T>(args)... {}

        /// \brief constructs all the vertex's elements
        /// \param args
        explicit constexpr vertex_impl(T &&... args) : rtk::vertex_elem<I, T>(std::move(args))... {}

        using rtk::vertex_elem<I, T>::get...;
    };

    /// \brief represents a vertex with its data
    /// \tparam S floating point to use for the vertex
    /// \tparam WS widths of the vectors storing data of the vertex (position, normal, colour, ...)
    template <std::floating_point S, scion::usize... WS>
        requires (sizeof...(WS) >= 1)
    struct vertex : rtk::vertex_impl<std::make_index_sequence<sizeof...(WS)>, rtk::vec<S, WS>...> {
        using rtk::vertex_impl<std::make_index_sequence<sizeof...(WS)>, rtk::vec<S, WS>...>::vertex_impl;

        static constexpr scion::usize SIZE = sizeof(S) * (0 + ... + WS);

        static_assert(sizeof(typename rtk::vertex<S, WS...>::vertex_impl) == sizeof(S) * (0 + ... + WS), "rtk::vertex has incoherent size (probably due to padding)");
    };
}
