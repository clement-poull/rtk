#ifndef RTK_MEASURES_D3_UTILS_INCLUDED
#define RTK_MEASURES_D3_UTILS_INCLUDED

#include "scion/types.hpp"

#include <optional>
#include <span>
#include <vector>

namespace rtk::measures::d3 {
    struct data_box_counting {
        std::vector<bool> boxes;
    };

    struct data {
        std::span<const scion::f32> vertices;

        scion::f32 meanline;

        scion::f32 x_min;
        scion::f32 x_max;
        scion::f32 y_min;
        scion::f32 y_max;

        std::optional<data_box_counting> box_counting;
    };

    auto compute_sa(rtk::measures::d3::data & data) -> scion::f32 {}

    auto compute_sq(rtk::measures::d3::data & data) -> scion::f32 {}

    auto compute_ssk(rtk::measures::d3::data & data) -> scion::f32 {}

    auto compute_sku(rtk::measures::d3::data & data) -> scion::f32 {}

    auto compute_fd(rtk::measures::d3::data & data) -> scion::f32 {}
}

#endif
