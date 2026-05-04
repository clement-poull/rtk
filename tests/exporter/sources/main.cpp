#include "scion/scion.hpp"
#include "scion/prelude/prelude.hpp"

#include "rtk/rtk.hpp"
#include "rtk/prelude/prelude.hpp"

#include <array>
#include <cmath>
#include <numbers>

auto main() -> int {
    std::array<rtk::vertex<f32, 2, 3>, 4> square = {
        rtk::vertex<f32, 2, 3>{
            vec2f{1.0f, 1.0f}.normalized(),
            vec3f{1.0f, 0.0f, 0.0f},
        }, rtk::vertex<f32, 2, 3>{
            vec2f{-1.0f, 1.0f}.normalized(),
            vec3f{0.5f, 0.0f, 0.5f},
        }, rtk::vertex<f32, 2, 3>{
            vec2f{-1.0f, -1.0f}.normalized(),
            vec3f{1.0f, 0.0f, 0.0f},
        }, rtk::vertex<f32, 2, 3>{
            vec2f{1.0f, -1.0f}.normalized(),
            vec3f{0.5f, 0.0f, 0.5f},
        },
    };

    if (auto result = rtk::export_obj_verts<2, 0, 5>("test_verts.obj", reinterpret_cast<f32 const*>(square.data()), std::views::iota(0UZ, square.size())); !result.has_value()) {
        rtk::log::error("test", "{}", result.error());
    }

    if (auto result = rtk::export_obj_lines<2, 0, 5>("test_lines.obj", reinterpret_cast<f32 const*>(square.data()), std::views::iota(0UZ, square.size()), std::views::iota(0UZ, square.size() + 1) | std::views::transform([&square](auto const& index){ return index % square.size(); }) | std::views::adjacent<2>); !result.has_value()) {
        rtk::log::error("test", "{}", result.error());
    }

    if (auto result = rtk::export_obj_line_strip<2, 0, 5>("test_line_strip.obj", reinterpret_cast<f32 const*>(square.data()), std::views::iota(0UZ, square.size())); !result.has_value()) {
        rtk::log::error("test", "{}", result.error());
    }

    {
        std::vector<scion::usize> indices;
        indices.emplace_back(0);
        indices.emplace_back(1);
        indices.emplace_back(2);

        indices.emplace_back(0);
        indices.emplace_back(2);
        indices.emplace_back(3);

        if (auto result = rtk::export_obj_faces<2, 0, 5>("test_faces.obj", reinterpret_cast<f32 const*>(square.data()), std::views::iota(0UZ, square.size()), indices | std::views::chunk(3) | std::views::transform([](auto chunk) -> std::tuple<scion::usize, scion::usize, scion::usize> {
            auto it = chunk.begin();
            if constexpr (std::ranges::sized_range<decltype(chunk)>) {
                if (std::ranges::size(chunk) == 3) {
                    return std::tuple(*it, *std::next(it), *std::next(it, 2));
                }
            }

            return std::tuple(0, 0, 0);
        })); !result.has_value()) {
            rtk::log::error("test", "{}", result.error());
        }
    }

    return 0;
}
