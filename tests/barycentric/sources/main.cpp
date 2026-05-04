#include "scion/scion.hpp"
#include "scion/prelude/prelude.hpp"

#include "rtk/rtk.hpp"
#include "rtk/prelude/prelude.hpp"

#include <vector>

auto main() -> int {
    rtk::projection projection_r2 = rtk::projection::make<2>(2, 2).value(); {
        projection_r2.points<2>().value()[0] = vec2f{-1.0f, 0.0f};
        projection_r2.points<2>().value()[1] = vec2f{1.0f, 0.0f};
    }

    std::vector<vecxf> points_b2;
    points_b2.push_back(vecxf{{0.9f, 0.1f}});
    points_b2.push_back(vecxf{{0.5f, 0.5f}});
    points_b2.push_back(vecxf{{0.1f, 0.9f}});

    std::vector<rtk::vecf<2>> points_r2 = projection_r2.project<2>(points_b2).value();

    for (usize i = 0; i < points_r2.size(); ++i) {
        rtk::log::info("", "({}, {}) -> ({}, {})", points_b2[i][0], points_b2[i][1], points_r2[i][0], points_r2[i][1]);
    }

    return 0;
}
