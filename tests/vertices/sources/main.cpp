#include "scion/scion.hpp"
#include "scion/prelude/prelude.hpp"

#include "rtk/rtk.hpp"
#include "rtk/prelude/prelude.hpp"

#include <array>
#include <cmath>
#include <numbers>

template <typename Scalar, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
struct std::formatter<Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>, char> {
    std::formatter<std::string, char> string_formatter;

    constexpr auto parse(std::format_parse_context& ctx) {
        return string_formatter.parse(ctx);
    }

    auto format(const Eigen::Matrix<Scalar, Rows, Cols, Options, MaxRows, MaxCols>& m, std::format_context& ctx) const {
        std::ostringstream oss;
        oss << m;
        return string_formatter.format(oss.str(), ctx);
    }
};

auto main() -> int {
    // square inside the unit circle
    // vertices of the main diagonal are purple
    // vertices of the main anti-diagonal are red
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

    int exit_code = 0;

    auto check = [&exit_code](auto const& a, auto const& b) {
        if (a != b) {
            rtk::log::error("assert", "test failed: {} != {}", a, b);

            exit_code += 1;
        }
    };

    auto check_epsilon = [&exit_code](auto const& e, auto const& a, auto const& b) {
        if (std::abs(a - b) > e) {
            rtk::log::error("assert", "test failed: {} != {}", a, b);

            exit_code += 1;
        }
    };

    // check if the positions are right on all vertices
    check(square[0].get<0>(), vec2f{1.0f, 1.0f}.normalized());
    check(square[1].get<0>(), vec2f{-1.0f, 1.0f}.normalized());
    check(square[2].get<0>(), vec2f{-1.0f, -1.0f}.normalized());
    check(square[3].get<0>(), vec2f{1.0f, -1.0f}.normalized());

    // check if the two first vertices are orthogonal
    check(square[0].get<0>().dot(square[1].get<0>()), 0.0f);

    auto const* pointer = reinterpret_cast<f32 const*>(square.data());

    // check if the order of the values is okay for the first vertex
    check_epsilon(0.000001f, pointer[0], std::cos(1.0f * std::numbers::pi_v<f32> / 4.0f));
    check_epsilon(0.000001f, pointer[1], std::sin(1.0f * std::numbers::pi_v<f32> / 4.0f));
    check(pointer[2], 1.0f);
    check(pointer[3], 0.0f);
    check(pointer[4], 0.0f);

    // check if the order of the values is okay for the second vertex
    check_epsilon(0.000001f, pointer[5], std::cos(3.0f * std::numbers::pi_v<f32> / 4.0f));
    check_epsilon(0.000001f, pointer[6], std::sin(3.0f * std::numbers::pi_v<f32> / 4.0f));
    check(pointer[7], 0.5f);
    check(pointer[8], 0.0f);
    check(pointer[9], 0.5f);

    // check if the size of the struct is as expected
    check(sizeof(square), sizeof(f32) * 2 * 4 + sizeof(f32) * 3 * 4);

    if (exit_code == 0) {
        rtk::log::info("assert", "test succeeded");
    }

    return exit_code;
}
