#pragma once

#include "rtk/core/vertices.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <expected>
#include <ranges>

/// \file rtk/util/exporter.hpp
/// \brief Provides data export tooling.

namespace rtk
{
    template <typename Range>
    concept range_indices = std::ranges::range<Range> && (std::is_same_v<std::ranges::range_value_t<Range>, scion::usize> || std::is_same_v<std::ranges::range_value_t<Range>, scion::u32>);

    template <typename Range, scion::usize N = 0>
    concept range_indices_tuple = std::ranges::range<Range> && scion::uniform_tuple<std::ranges::range_value_t<Range>, N> && (std::is_same_v<std::tuple_element_t<0, std::ranges::range_value_t<Range>>, scion::usize> || std::is_same_v<std::tuple_element_t<0, std::ranges::range_value_t<Range>>, scion::u32>);

    /// \brief writes an obj file
    /// \tparam VP f32 per vertex position
    /// \tparam VN f32 per vertex normal
    /// \tparam RVerts range type to use for the indices
    /// \param path
    /// \param verts
    /// \param verts_indices indices to use in verts (use std::views::iota(0UZ, vertices.size()))
    template <scion::usize VP = 3, scion::usize VN = 0, scion::usize STRIDE = VP + VN, rtk::range_indices RVerts>
        requires (STRIDE >= VP + VN && VP <= 4 && VN <= 4)
    auto export_obj_verts(std::string const & path, scion::f32 const * verts, RVerts verts_indices) -> std::expected<void, std::string> {
        std::ofstream file(path);
        if (!file.is_open()) {
            return std::unexpected<std::string>(std::format("could not open file {}", path));
        }

        for (std::ranges::range_value_t<RVerts> const index : verts_indices) {
            scion::f32 const * const vertex = verts + index * STRIDE;

            if constexpr (VP != 0) {
                scion::f32 const * const vertex_p = verts + index * STRIDE;

                file << "v";

                for (scion::usize i = 0; i < VP; ++i) {
                    file << " " << vertex_p[i];
                }

                if constexpr (VP == 1) {
                    file << " " << 0;
                }

                if constexpr (VP == 2) {
                    file << " " << 0;
                }

                file << "\n";
            }

            if constexpr (VN != 0) {
                scion::f32 const * const vertex_n = verts + index * STRIDE + VP;

                file << "vn";

                for (scion::usize i = 0; i < VN; ++i) {
                    file << " " << vertex_n[i];
                }

                if constexpr (VN == 1) {
                    file << " " << 0;
                }

                if constexpr (VN == 2) {
                    file << " " << 0;
                }

                file << "\n";
            }
        }

        file.close();
        if (file.fail()) {
            return std::unexpected<std::string>(std::format("could not close file {}", path));
        }

        return {};
    }

    template <scion::usize VP = 3, scion::usize VN = 0, scion::usize STRIDE = VP + VN, rtk::range_indices RVerts, rtk::range_indices_tuple<2> RLines>
        requires (STRIDE >= VP + VN && VP <= 4 && VN <= 4)
    auto export_obj_lines(std::string const & path, scion::f32 const * verts, RVerts verts_indices, RLines lines_indices) -> std::expected<void, std::string> {
        std::ofstream file(path);
        if (!file.is_open()) {
            return std::unexpected<std::string>(std::format("could not open file {}", path));
        }

        for (std::ranges::range_value_t<RVerts> const index : verts_indices) {
            scion::f32 const * const vertex = verts + index * STRIDE;

            if constexpr (VP != 0) {
                scion::f32 const * const vertex_p = verts + index * STRIDE;

                file << "v";

                for (scion::usize i = 0; i < VP; ++i) {
                    file << " " << vertex_p[i];
                }

                if constexpr (VP == 1) {
                    file << " " << 0;
                }

                if constexpr (VP == 2) {
                    file << " " << 0;
                }

                file << "\n";
            }

            if constexpr (VN != 0) {
                scion::f32 const * const vertex_n = verts + index * STRIDE + VP;

                file << "vn";

                for (scion::usize i = 0; i < VN; ++i) {
                    file << " " << vertex_n[i];
                }

                if constexpr (VN == 1) {
                    file << " " << 0;
                }

                if constexpr (VN == 2) {
                    file << " " << 0;
                }

                file << "\n";
            }
        }

        for (std::ranges::range_value_t<RLines> const index : lines_indices) {
            file << "l " << 1 + std::get<0>(index) << " " << 1 + std::get<1>(index) << "\n";
        }

        file.close();
        if (file.fail()) {
            return std::unexpected<std::string>(std::format("could not close file {}", path));
        }

        return {};
    }

    template <scion::usize VP = 3, scion::usize VN = 0, scion::usize STRIDE = VP + VN, rtk::range_indices RVerts>
        requires (STRIDE >= VP + VN && VP <= 4 && VN <= 4)
    auto export_obj_line_strip(std::string const & path, scion::f32 const * verts, RVerts verts_indices) -> std::expected<void, std::string> {
        std::ofstream file(path);
        if (!file.is_open()) {
            return std::unexpected<std::string>(std::format("could not open file {}", path));
        }

        std::stringstream ss_line;
        scion::usize count_vertices = 0;
        ss_line << "l";

        for (std::ranges::range_value_t<RVerts> const index : verts_indices) {
            scion::f32 const * const vertex = verts + index * STRIDE;

            if constexpr (VP != 0) {
                scion::f32 const * const vertex_p = verts + index * STRIDE;

                file << "v";

                for (scion::usize i = 0; i < VP; ++i) {
                    file << " " << vertex_p[i];
                }

                if constexpr (VP == 1) {
                    file << " " << 0;
                }

                if constexpr (VP == 2) {
                    file << " " << 0;
                }

                file << "\n";
            }

            if constexpr (VN != 0) {
                scion::f32 const * const vertex_n = verts + index * STRIDE + VP;

                file << "vn";

                for (scion::usize i = 0; i < VN; ++i) {
                    file << " " << vertex_n[i];
                }

                if constexpr (VN == 1) {
                    file << " " << 0;
                }

                if constexpr (VN == 2) {
                    file << " " << 0;
                }

                file << "\n";
            }

            ss_line << " " << ++count_vertices;
        }

        ss_line << "\n";

        if (count_vertices != 0) {
            file << ss_line.str();
        }

        file.close();
        if (file.fail()) {
            return std::unexpected<std::string>(std::format("could not close file {}", path));
        }

        return {};
    }

    template <scion::usize VP = 3, scion::usize VN = 0, scion::usize STRIDE = VP + VN, rtk::range_indices RVerts, rtk::range_indices_tuple RFaces>
        requires (STRIDE >= VP + VN && VP <= 4 && VN <= 4)
    auto export_obj_faces(std::string const & path, scion::f32 const * verts, RVerts verts_indices, RFaces faces_indices) -> std::expected<void, std::string> {
        std::ofstream file(path);
        if (!file.is_open()) {
            return std::unexpected<std::string>(std::format("could not open file {}", path));
        }

        for (std::ranges::range_value_t<RVerts> const index : verts_indices) {
            scion::f32 const * const vertex = verts + index * STRIDE;

            if constexpr (VP != 0) {
                scion::f32 const * const vertex_p = verts + index * STRIDE;

                file << "v";

                for (scion::usize i = 0; i < VP; ++i) {
                    file << " " << vertex_p[i];
                }

                if constexpr (VP == 1) {
                    file << " " << 0;
                }

                if constexpr (VP == 2) {
                    file << " " << 0;
                }

                file << "\n";
            }

            if constexpr (VN != 0) {
                scion::f32 const * const vertex_n = verts + index * STRIDE + VP;

                file << "vn";

                for (scion::usize i = 0; i < VN; ++i) {
                    file << " " << vertex_n[i];
                }

                if constexpr (VN == 1) {
                    file << " " << 0;
                }

                if constexpr (VN == 2) {
                    file << " " << 0;
                }

                file << "\n";
            }
        }

        for (std::ranges::range_value_t<RFaces> const index : faces_indices) {
            file << "f";
            std::apply([&file](auto&&... i) {
                if constexpr (VN == 0) {
                    ((file << " " << i + 1), ...);
                } else {
                    ((file << " " << i + 1 << "//" << i + 1), ...);
                }
            }, index);
            file << "\n";
        }

        file.close();
        if (file.fail()) {
            return std::unexpected<std::string>(std::format("could not close file {}", path));
        }

        return {};
    }
}
