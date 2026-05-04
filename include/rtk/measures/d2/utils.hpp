#ifndef RTK_MEASURES_D2_UTILS_INCLUDED
#define RTK_MEASURES_D2_UTILS_INCLUDED

#include "scion/core/types.hpp"

#include <optional>
#include <span>
#include <vector>

namespace rtk::measures::d2 {
    struct data_box_counting {
        static constexpr scion::usize MIN_ITERATION_REGRESSION = 4;
        static constexpr scion::usize MAX_ITERATION_REGRESSION = 8;

        std::optional<scion::usize> min_iteration_regression;
        std::optional<scion::usize> max_iteration_regression;

        std::vector<std::vector<bool>> boxes;

        scion::f32 value;
    };

    struct data {
        std::span<const rtk::vec2f> vertices;

        bool sorted;

        scion::f32 meanline;

        scion::f32 x_min;
        scion::f32 x_max;
        scion::f32 y_min;
        scion::f32 y_max;

        std::optional<data_box_counting> box_counting;
    };

    auto preprocess(std::span<const rtk::vec2f> vertices) -> std::optional<rtk::measures::d2::data> {
        if (vertices.size() <= 1) return std::nullopt;

        scion::f32 integral = 0.0f;
        scion::f32 length = 0.0f;

        scion::f32 x_min = vertices[0].x();
        scion::f32 x_max = vertices[0].x();
        scion::f32 y_min = vertices[0].y();
        scion::f32 y_max = vertices[0].y();

        bool sorted = true;

        for (scion::usize idx = 1; idx < vertices.size(); ++idx) {
            rtk::vec2f const& ante = vertices[idx - 1];
            rtk::vec2f const& post = vertices[idx];

            integral += ((post.y() + ante.y()) / 2.0f) * std::abs(post.x() - ante.x());
            length += std::abs(post.x() - ante.x());
            x_min = std::min(x_min, post.x());
            x_max = std::max(x_max, post.x());
            y_min = std::min(y_min, post.y());
            y_max = std::max(y_max, post.y());
            sorted &= ante.x() <= post.x();
        }

        return rtk::measures::d2::data{
            .vertices = vertices,

            .sorted = sorted,

            .meanline = (length > 0.0f) ? integral / length : ((y_max + y_min) / 2.0f),

            .x_min = x_min,
            .x_max = x_max,
            .y_min = y_min,
            .y_max = y_max,
        };
    }

    auto compute_ra(rtk::measures::d2::data & data) -> std::optional<scion::f32> {
        if (!data.sorted) return std::nullopt;
        if (data.x_max <= data.x_min) return std::nullopt;

        scion::f32 acc = 0.0f;

        for (scion::usize idx = 1; idx < data.vertices.size(); ++idx) {
            rtk::vec2f const& ante = data.vertices[idx - 1];
            rtk::vec2f const& post = data.vertices[idx];

            scion::f32 const xante = ante.x();
            scion::f32 const yante = ante.y() - data.meanline;
            scion::f32 const yante_abs = std::abs(yante);

            scion::f32 const xpost = post.x();
            scion::f32 const ypost = post.y() - data.meanline;
            scion::f32 const ypost_abs = std::abs(ypost);

            scion::f32 const dx = xpost - xante;

            acc += dx * (yante_abs + ypost_abs) / 2.0f;
        }

        return acc / (data.x_max - data.x_min);
    }

    auto compute_rq(rtk::measures::d2::data & data) -> std::optional<scion::f32> {
        if (!data.sorted) return std::nullopt;
        if (data.x_max <= data.x_min) return std::nullopt;

        scion::f32 acc_2 = 0.0f;

        for (scion::usize idx = 1; idx < data.vertices.size(); ++idx) {
            rtk::vec2f const& ante = data.vertices[idx - 1];
            rtk::vec2f const& post = data.vertices[idx];

            scion::f32 const xante = ante.x();
            scion::f32 const yante = ante.y() - data.meanline;
            scion::f32 const yante_2 = yante * yante;

            scion::f32 const xpost = post.x();
            scion::f32 const ypost = post.y() - data.meanline;
            scion::f32 const ypost_2 = ypost * ypost;

            scion::f32 const dx = xpost - xante;

            acc_2 += dx * (yante_2 + ypost_2) / 2.0f;
        }

        return std::sqrt(acc_2 / (data.x_max - data.x_min));
    }

    auto compute_rsk(rtk::measures::d2::data & data) -> std::optional<scion::f32> {
        if (!data.sorted) return std::nullopt;
        if (data.x_max <= data.x_min) return std::nullopt;

        scion::f32 acc_2 = 0.0f;
        scion::f32 acc_3 = 0.0f;

        for (scion::usize idx = 1; idx < data.vertices.size(); ++idx) {
            rtk::vec2f const& ante = data.vertices[idx - 1];
            rtk::vec2f const& post = data.vertices[idx];

            scion::f32 const xante = ante.x();
            scion::f32 const yante = ante.y() - data.meanline;
            scion::f32 const yante_2 = yante * yante;
            scion::f32 const yante_3 = yante * yante * yante;

            scion::f32 const xpost = post.x();
            scion::f32 const ypost = post.y() - data.meanline;
            scion::f32 const ypost_2 = ypost * ypost;
            scion::f32 const ypost_3 = ypost * ypost * ypost;

            scion::f32 const dx = xpost - xante;

            acc_2 += dx * (yante_2 + ypost_2) / 2.0f;
            acc_3 += dx * (yante_3 + ypost_3) / 2.0f;
        }

        scion::f32 const rq = std::sqrt(acc_2 /  (data.x_max - data.x_min));

        return (acc_3 / (data.x_max - data.x_min)) / (rq * rq * rq);
    }

    auto compute_rku(rtk::measures::d2::data & data) -> std::optional<scion::f32> {
        if (!data.sorted) return std::nullopt;
        if (data.x_max <= data.x_min) return std::nullopt;

        scion::f32 acc_2 = 0.0f;
        scion::f32 acc_4 = 0.0f;

        for (scion::usize idx = 1; idx < data.vertices.size(); ++idx) {
            rtk::vec2f const& ante = data.vertices[idx - 1];
            rtk::vec2f const& post = data.vertices[idx];

            scion::f32 const xante = ante.x();
            scion::f32 const yante = ante.y() - data.meanline;
            scion::f32 const yante_2 = yante * yante;
            scion::f32 const yante_4 = yante * yante * yante * yante;

            scion::f32 const xpost = post.x();
            scion::f32 const ypost = post.y() - data.meanline;
            scion::f32 const ypost_2 = ypost * ypost;
            scion::f32 const ypost_4 = ypost * ypost * ypost * ypost;

            scion::f32 const dx = xpost - xante;

            acc_2 += dx * (yante_2 + ypost_2) / 2.0f;
            acc_4 += dx * (yante_4 + ypost_4) / 2.0f;
        }

        scion::f32 const rq = std::sqrt(acc_2 /  (data.x_max - data.x_min));

        return (acc_4 / (data.x_max - data.x_min)) / (rq * rq * rq * rq);
    }

    auto compute_fd(rtk::measures::d2::data & data) -> std::optional<scion::f32> {
        data.box_counting = rtk::measures::d2::data_box_counting{};

        data.box_counting.value().value = 0.0f;
        data.box_counting.value().boxes.clear();
        data.box_counting.value().min_iteration_regression.reset();
        data.box_counting.value().max_iteration_regression.reset();

        if (data.x_max <= data.x_min || data.y_max <= data.y_min) return std::nullopt;

        scion::f32 const lx = data.x_max - data.x_min;
        scion::f32 const ly = data.y_max - data.y_min;
        scion::f32 const scale = std::max(lx, ly);

        std::vector<rtk::vec2f> vertices_normalised;
        vertices_normalised.reserve(data.vertices.size());

        for (rtk::vec2f const& v : data.vertices) {
            scion::f32 const nx = (v.x() - data.x_min) / scale;
            scion::f32 const ny = (v.y() - data.y_min) / scale;
            vertices_normalised.emplace_back(nx, ny);
        }

        std::vector<scion::f32> log_inv_eps;
        std::vector<scion::f32> log_N;

        scion::usize const max_scales = 8;
        scion::f32 eps = 0.5f;

        for (int iter = 0; iter < max_scales; ++iter) {
            if (eps < 1.0f / (vertices_normalised.size() * 4)) break;

            scion::usize grid_res = static_cast<scion::usize>(std::ceil(1.0f / eps));

            std::vector<bool> occupied(grid_res * grid_res, false);
            scion::usize box_count = 0;

            auto mark_cell = [&](scion::isize ix, scion::isize iy) {
                if (ix < 0 || iy < 0) return;
                if (ix >= static_cast<scion::isize>(grid_res)) return;
                if (iy >= static_cast<scion::isize>(grid_res)) return;

                scion::usize index = static_cast<scion::usize>(iy) * grid_res + static_cast<scion::usize>(ix);

                if (!occupied[index]) {
                    occupied[index] = true;
                    ++box_count;
                }
            };

            auto traverse_segment = [&](scion::f32 x1, scion::f32 y1, scion::f32 x2, scion::f32 y2) {
                auto clamp_to_unit = [](scion::f32 v) -> scion::f32 {
                    if (v < 0.0f) return 0.0f;
                    if (v > 1.0f) return 1.0f;

                    return v;
                };

                x1 = clamp_to_unit(x1);
                y1 = clamp_to_unit(y1);
                x2 = clamp_to_unit(x2);
                y2 = clamp_to_unit(y2);

                auto cell_of = [eps, grid_res](scion::f32 v) -> scion::isize {
                    scion::isize i = static_cast<scion::isize>(std::floor(v / eps));
                    if (i < 0) i = 0;
                    if (i >= static_cast<scion::isize>(grid_res)) i = static_cast<scion::isize>(grid_res) - 1;

                    return i;
                };

                scion::isize ix = cell_of(x1);
                scion::isize iy = cell_of(y1);
                scion::isize ix_end = cell_of(x2);
                scion::isize iy_end = cell_of(y2);

                mark_cell(ix, iy);

                scion::f32 dx = x2 - x1;
                scion::f32 dy = y2 - y1;

                scion::isize step_x = (dx > 0) ? 1 : (dx < 0 ? -1 : 0);
                scion::isize step_y = (dy > 0) ? 1 : (dy < 0 ? -1 : 0);

                scion::f32 t_delta_x = (dx != 0.0f) ? std::abs(eps / dx) : std::numeric_limits<scion::f32>::infinity();
                scion::f32 t_delta_y = (dy != 0.0f) ? std::abs(eps / dy) : std::numeric_limits<scion::f32>::infinity();

                scion::f32 next_boundary_x = (step_x > 0) ? (static_cast<scion::f32>(ix + 1) * eps) : (static_cast<scion::f32>(ix) * eps);
                scion::f32 next_boundary_y = (step_y > 0) ? (static_cast<scion::f32>(iy + 1) * eps) : (static_cast<scion::f32>(iy) * eps);

                scion::f32 t_max_x = (dx != 0.0f) ? ((next_boundary_x - x1) / dx) : std::numeric_limits<scion::f32>::infinity();
                scion::f32 t_max_y = (dy != 0.0f) ? ((next_boundary_y - y1) / dy) : std::numeric_limits<scion::f32>::infinity();

                if (t_max_x < 0.0f) t_max_x = 0.0f;
                if (t_max_y < 0.0f) t_max_y = 0.0f;

                while (ix != ix_end || iy != iy_end) {
                    if (t_max_x < t_max_y) {
                        ix += step_x;
                        t_max_x += t_delta_x;
                    } else if (t_max_y < t_max_x) {
                        iy += step_y;
                        t_max_y += t_delta_y;
                    } else {
                        ix += step_x;
                        iy += step_y;
                        t_max_x += t_delta_x;
                        t_max_y += t_delta_y;
                    }

                    mark_cell(ix, iy);
                }
            };

            for (scion::usize idx = 0; idx < vertices_normalised.size(); idx += 2) {
                rtk::vec2f const& ante = vertices_normalised[idx];
                rtk::vec2f const& post = vertices_normalised[idx + 1];

                traverse_segment(ante.x(), ante.y(), post.x(), post.y());
            }

            if (box_count > 1) {
                data.box_counting.value().boxes.push_back(std::move(occupied));

                log_inv_eps.push_back(std::log(1.0f / eps));
                log_N.push_back(std::log(static_cast<scion::f32>(box_count)));
            }

            eps *= 0.5f;
        }

        if (log_N.size() < 5) {
            data.box_counting.value().boxes.clear();
            return std::nullopt;
        }

        scion::f32 sum_x = 0.0f, sum_y = 0.0f, sum_xy = 0.0f, sum_x2 = 0.0f;
        const size_t n = log_N.size();

        for (size_t i = 0; i < n; ++i) {
            scion::f32 x = log_inv_eps[i];
            scion::f32 y = log_N[i];
            sum_x  += x;
            sum_y  += y;
            sum_xy += x * y;
            sum_x2 += x * x;
        }

        scion::f32 denom = static_cast<scion::f32>(n) * sum_x2 - sum_x * sum_x;
        if (std::abs(denom) < 1e-8f) {
            data.box_counting.value().boxes.clear();
            return std::nullopt;
        }

        scion::f32 slope = (static_cast<scion::f32>(n) * sum_xy - sum_x * sum_y) / denom;

        data.box_counting.value().value = slope;
        data.box_counting.value().min_iteration_regression = 0;
        data.box_counting.value().max_iteration_regression = n - 1;

        return slope;
    }
}

#endif
