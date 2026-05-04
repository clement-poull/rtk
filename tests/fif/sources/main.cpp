#include "scion/scion.hpp"
#include "scion/prelude/prelude.hpp"

#include "rtk/rtk.hpp"
#include "rtk/model/model.hpp"
#include "rtk/prelude/prelude.hpp"

#include <iostream>

auto main() -> int {
    std::vector<vec2f> points;
    points.emplace_back(0.000f,  0.0f);
    points.emplace_back(0.333f,  0.2f);
    points.emplace_back(0.666f, -0.2f);
    points.emplace_back(0.999f,  0.0f);
    usize const size = points.size();
    rtk::fif model = rtk::fif::make(std::move(points)).value();

    // std::vector<f32> as = std::views::iota(0uz, size - 1uz) | std::views::transform([&model](auto const& idx) -> f32 { return model.a(idx); }) | std::ranges::to<std::vector>();
    // std::vector<f32> cs = std::views::iota(0uz, size - 1uz) | std::views::transform([&model](auto const& idx) -> f32 { return model.c(idx); }) | std::ranges::to<std::vector>();
    // std::vector<f32> ds = std::views::iota(0uz, size - 1uz) | std::views::transform([&model](auto const& idx) -> f32 { return model.d(idx); }) | std::ranges::to<std::vector>();
    // std::vector<f32> es = std::views::iota(0uz, size - 1uz) | std::views::transform([&model](auto const& idx) -> f32 { return model.e(idx); }) | std::ranges::to<std::vector>();
    // std::vector<f32> fs = std::views::iota(0uz, size - 1uz) | std::views::transform([&model](auto const& idx) -> f32 { return model.f(idx); }) | std::ranges::to<std::vector>();

    std::vector<f32> as;
    std::vector<f32> cs;
    std::vector<f32> ds;
    std::vector<f32> es;
    std::vector<f32> fs;

    for (auto idx = 0uz; idx < size - 1uz; ++idx) {
        as.push_back(model.a(idx));
        cs.push_back(model.c(idx));
        ds.push_back(model.d(idx));
        es.push_back(model.e(idx));
        fs.push_back(model.f(idx));
    }

    {
        rtk::log::info("test", "checking contractions");

        // check a valid model

        if (!rtk::fif::check_contractions(as, cs, ds)) {
            rtk::log::error("test", "failed at {}", __LINE__);
        }

        // check an invalid model

        cs[1] = 10.0f;

        if (rtk::fif::check_contractions(as, cs, ds)) {
            rtk::log::error("test", "failed at {}", __LINE__);
        }

        cs[1] = model.c(1);

        rtk::log::info("test", "checking pseudo tangents");

        if (!rtk::fif::compute_pseudo_tangent_r(as, cs, ds, 0)->normalized().isApprox(rtk::fif::compute_tangent(as[0], cs[0], ds[0])->normalized())) {
            rtk::log::error("test", "failed at {}", __LINE__);
        }

        if (!rtk::fif::compute_pseudo_tangent_l(as, cs, ds, size - 1)->normalized().isApprox(rtk::fif::compute_tangent(as[size - 2], cs[size - 2], ds[size - 2])->normalized())) {
            rtk::log::error("test", "failed at {}", __LINE__);
        }
    }

    return 0;
}
