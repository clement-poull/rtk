#pragma once

#include "scion/core/core.hpp"

#include "rtk/core/types.hpp"

#include "rtk/util/log.hpp"
#include "rtk/util/assert.hpp"

#include <utility>
#include <expected>

/// \file rtk/model/fif.hpp
/// \brief Provides a set of function to work with Fractal Interpolation Functions and a class to represent a Fractal Interpolation Function and enforcing invariants.
///
/// Most functions use debug assertions.
/// Such assertions are enabled in consteval contexts or in debug builds.

namespace rtk {
    enum fif_error {
        not_contractive,
    };

    /// \brief Represents a Fractal Interpolation Function.
    /// \tparam Scalar The scalar type used for the coefficients.
    template <typename Scalar = scion::f32>
    class fif_impl {
    private: // static private members
        static inline constexpr void check_llrr(rtk::vec<Scalar, 2> const& L, rtk::vec<Scalar, 2> const& l, rtk::vec<Scalar, 2> const& r, rtk::vec<Scalar, 2> const& R, std::source_location const location = std::source_location::current()) {
            scion::assert_debug(rtk::is_vec_valid(L), "invalid interpolation point", location);
            scion::assert_debug(rtk::is_vec_valid(l), "invalid interpolation point", location);
            scion::assert_debug(rtk::is_vec_valid(r), "invalid interpolation point", location);
            scion::assert_debug(rtk::is_vec_valid(R), "invalid interpolation point", location);

            scion::assert_debug(l.x() < r.x(), "unsorted interpolation points", location);
            scion::assert_debug(L == l || L.x() < l.x(), "unsorted interpolation points", location);
            scion::assert_debug(r == R || r.x() < R.x(), "unsorted interpolation points", location);
        }

        static inline constexpr void check_points(std::span<rtk::vec<Scalar, 2> const> const points, std::source_location const location = std::source_location::current()) {
            for (scion::usize idx = 0; idx < points.size(); ++idx) scion::assert_debug(rtk::is_vec_valid(points[idx]), "invalid interpolation point", location);

            for (scion::usize idx = 0; idx < points.size() - 1; ++idx) scion::assert_debug(points[idx].x() < points[idx + 1].x(), "unsorted interpolation points", location);
        }

        static inline constexpr void check_a(Scalar const a, std::source_location const location = std::source_location::current()) {
            scion::assert_debug(rtk::is_val_valid(a), "invalid transformation entry", location);

            scion::assert_debug(a >= 0 && a < 1, "invalid transformation entry", location);
        }

        static inline constexpr void check_c(Scalar const c, std::source_location const location = std::source_location::current()) {
            scion::assert_debug(rtk::is_val_valid(c), "invalid transformation entry", location);
        }

        static inline constexpr void check_d(Scalar const d, std::source_location const location = std::source_location::current()) {
            scion::assert_debug(rtk::is_val_valid(d), "invalid transformation entry", location);

            scion::assert_debug(std::abs(d) < 1, "invalid transformation entry", location);
        }

        static inline constexpr void check_e(Scalar const e, std::source_location const location = std::source_location::current()) {
            scion::assert_debug(rtk::is_val_valid(e), "invalid transformation entry", location);
        }

        static inline constexpr void check_f(Scalar const f, std::source_location const location = std::source_location::current()) {
            scion::assert_debug(rtk::is_val_valid(f), "invalid transformation entry", location);
        }

        static inline constexpr void check_acd(Scalar const a, Scalar const c, Scalar const d, std::source_location const location = std::source_location::current()) {
            rtk::fif_impl<Scalar>::check_a(a, location);
            rtk::fif_impl<Scalar>::check_c(c, location);
            rtk::fif_impl<Scalar>::check_d(d, location);
        }

        static inline constexpr void check_acdef(Scalar const a, Scalar const c, Scalar const d, Scalar const e, Scalar const f, std::source_location const location = std::source_location::current()) {
            rtk::fif_impl<Scalar>::check_a(a, location);
            rtk::fif_impl<Scalar>::check_c(c, location);
            rtk::fif_impl<Scalar>::check_d(d, location);
            rtk::fif_impl<Scalar>::check_e(e, location);
            rtk::fif_impl<Scalar>::check_f(f, location);
        }

    protected: // static protected members

    public: // static public members
        /// \brief Checks that the sizes of the views are consistent, and returns it if so.
        /// \param head First coefficients view.
        /// \param tail Rest of the coefficients views.
        /// \return Size of all views, if it is unique.
        template <typename Head, typename... Tail>
            requires (std::convertible_to<Head, std::span<typename std::remove_reference_t<Head>::element_type>> && (std::convertible_to<Tail, std::span<typename std::remove_reference_t<Tail>::element_type>> && ...))
        [[nodiscard]] static inline constexpr auto check_coefficients(Head && head, Tail &&... tail) -> std::optional<scion::usize> {
            if ((... && (tail.size() == head.size()))) {
                return head.size();
            }

            return std::nullopt;
        }

        /// \brief Computes the value of the `a` coefficient for the given points.
        /// \param L Left extremity of the curve.
        /// \param l Left extremity of the trans.
        /// \param r Right extremity of the trans.
        /// \param R Right extremity of the curve.
        /// \return `a` coefficient for the given points.
        [[nodiscard]] static inline constexpr auto compute_a(rtk::vec<Scalar, 2> const& L, rtk::vec<Scalar, 2> const& l, rtk::vec<Scalar, 2> const& r, rtk::vec<Scalar, 2> const& R) -> Scalar {
            rtk::fif_impl<Scalar>::check_llrr(L, l, r, R);

            Scalar const gdx = R.x() - L.x();
            Scalar const ldx = r.x() - l.x();

            return ldx / gdx;
        }

        /// \brief Computes the value of the `c` coefficient for the given points and vertical scaling factor.
        /// \param L Left extremity of the curve.
        /// \param l Left extremity of the trans.
        /// \param r Right extremity of the trans.
        /// \param R Right extremity of the curve.
        /// \param d Vertical scaling factor of the trans.
        /// \return `c` coefficient for the given points and vertical scaling factor.
        [[nodiscard]] static inline constexpr auto compute_c(rtk::vec<Scalar, 2> const& L, rtk::vec<Scalar, 2> const& l, rtk::vec<Scalar, 2> const& r, rtk::vec<Scalar, 2> const& R, Scalar const d) -> Scalar {
            rtk::fif_impl<Scalar>::check_llrr(L, l, r, R);
            rtk::fif_impl<Scalar>::check_d(d);

            Scalar const gdx = R.x() - L.x();
            Scalar const gdy = R.y() - L.y();

            Scalar const ldy = r.y() - l.y();

            return ldy / gdx - d * gdy / gdx;
        }

        /// \brief Computes the value of the `e` coefficient for the given points
        /// \param L Left extremity of the curve.
        /// \param l Left extremity of the trans.
        /// \param r Right extremity of the trans.
        /// \param R Right extremity of the curve.
        /// \return `e` coefficient for the given points.
        [[nodiscard]] static inline constexpr auto compute_e(rtk::vec<Scalar, 2> const& L, rtk::vec<Scalar, 2> const& l, rtk::vec<Scalar, 2> const& r, rtk::vec<Scalar, 2> const& R) -> Scalar {
            rtk::fif_impl<Scalar>::check_llrr(L, l, r, R);

            Scalar const gdx = R.x() - L.x();

            Scalar const Rxlx = R.x() * l.x();
            Scalar const Lxrx = L.x() * r.x();

            return (Rxlx - Lxrx) / gdx;
        }

        /// \brief Computes the value of the `f` coefficient for the given points and vertical scaling factor.
        /// \param L Left extremity of the curve.
        /// \param l Left extremity of the trans.
        /// \param r Right extremity of the trans.
        /// \param R Right extremity of the curve.
        /// \param d Vertical scaling factor of the trans.
        /// \return `f` coefficients for the given points and vertical scaling factor.
        [[nodiscard]] static inline constexpr auto compute_f(rtk::vec<Scalar, 2> const& L, rtk::vec<Scalar, 2> const& l, rtk::vec<Scalar, 2> const& r, rtk::vec<Scalar, 2> const& R, Scalar const d) -> Scalar {
            rtk::fif_impl<Scalar>::check_llrr(L, l, r, R);
            rtk::fif_impl<Scalar>::check_d(d);

            Scalar const gdx = R.x() - L.x();

            Scalar const RxLy = R.x() * L.y();
            Scalar const LxRy = L.x() * R.y();

            Scalar const Rxly = R.x() * l.y();
            Scalar const Lxry = L.x() * r.y();

            return (Rxly - Lxry) / gdx - d * ((RxLy - LxRy) / gdx);
        }

        /// \brief Computes the `a` coefficients in-place.
        /// \param points Interpolation points of the FIF.
        /// \param as `a` coefficients store.
        /// \throws std::runtime_error points.size() must be greater than 1.
        /// \throws std::runtime_error the sizes must respect points.size() - 1 == as.size().
        static inline constexpr void compute_as(std::span<const rtk::vec<Scalar, 2>> const points, std::span<Scalar> const as) {
            rtk::fif_impl<Scalar>::check_points(points);

            if (points.size() <= 1) throw std::runtime_error("degenerate fif");
            if (as.size() != points.size() - 1) throw std::runtime_error("inconsistent sizes");

            scion::usize const size = points.size() - 1;

            for (scion::usize idx = 0; idx < size; ++idx) as[idx] = compute_a(points.front(), points[idx], points[idx + 1], points.back());
        }

        /// \brief Computes the `e` coefficients in-place.
        /// \param points Interpolation points of the FIF.
        /// \param es `e` coefficients store.
        /// \throws std::runtime_error points.size() must be greater than 1.
        /// \throws std::runtime_error the sizes must respect points.size() - 1 == es.size().
        static inline constexpr void compute_es(std::span<const rtk::vec<Scalar, 2>> const points, std::span<Scalar> const es) {
            rtk::fif_impl<Scalar>::check_points(points);

            if (points.size() <= 1) throw std::runtime_error("degenerate fif");
            if (es.size() != points.size() - 1) throw std::runtime_error("inconsistent sizes");

            scion::usize const size = points.size() - 1;

            for (scion::usize idx = 0; idx < size; ++idx) es[idx] = compute_e(points.front(), points[idx], points[idx + 1], points.back());
        }

        /// \brief Computes the `c`coefficients in-place.
        /// \param points Interpolation points of the FIF.
        /// \param ds Vertical scaling factors view.
        /// \param cs `c` coefficients store.
        /// \throws std::runtime_error points.size() must be greater than 1.
        /// \throws std::runtime_error the sizes must respect points.size() - 1 == ds.size() == cs.size().
        static inline constexpr void compute_cs(std::span<const rtk::vec<Scalar, 2>> const points, std::span<const Scalar> const ds, std::span<Scalar> const cs) {
            rtk::fif_impl<Scalar>::check_points(points);

            if (points.size() <= 1) throw std::runtime_error("degenerate fif");
            if (ds.size() != points.size() - 1 || cs.size() != points.size() - 1) throw std::runtime_error("inconsistent sizes");

            scion::usize const size = points.size() - 1;

            for (scion::usize idx = 0; idx < size; ++idx) cs[idx] = compute_c(points.front(), points[idx], points[idx + 1], points.back(), ds[idx]);
        }

        /// \brief Computes the `f`coefficients in-place.
        /// \param points Interpolation points of the FIF.
        /// \param ds Vertical scaling factors view.
        /// \param fs `f` coefficients store.
        /// \throws std::runtime_error points.size() must be greater than 1.
        /// \throws std::runtime_error the sizes must respect points.size() - 1 == ds.size() == fs.size().
        static inline constexpr void compute_fs(std::span<const rtk::vec<Scalar, 2>> const points, std::span<const Scalar> const ds, std::span<Scalar> const fs) {
            rtk::fif_impl<Scalar>::check_points(points);

            if (points.size() <= 1) throw std::runtime_error("degenerate fif");
            if (ds.size() != points.size() - 1 || fs.size() != points.size() - 1) throw std::runtime_error("inconsistent sizes");

            scion::usize const size = points.size() - 1;

            for (scion::usize idx = 0; idx < size; ++idx) fs[idx] = compute_f(points.front(), points[idx], points[idx + 1], points.back(), ds[idx]);
        }

        /// \brief Computes the `a` and `e` coefficients in-place.
        /// \param points Interpolation points of the FIF.
        /// \param as `a` coefficients store.
        /// \param es `e` coefficients store.
        /// \throws std::runtime_error points.size() must be greater than 1.
        /// \throws std::runtime_error the sizes must respect points.size() - 1 == as.size() == es.size().
        static inline constexpr void compute_aes(std::span<const rtk::vec<Scalar, 2>> const points, std::span<Scalar> const as, std::span<Scalar> const es) {
            if (points.size() <= 1) throw std::runtime_error("degenerate fif");
            if (as.size() != points.size() - 1 || es.size() != points.size() - 1) throw std::runtime_error("inconsistent sizes");

            scion::usize const size = points.size() - 1;

            rtk::vec<Scalar, 2> const L = points.front();
            rtk::vec<Scalar, 2> const R = points.back();

            Scalar const gdx = R.x() - L.x();

            for (scion::usize idx = 0; idx < size; ++idx) {
                rtk::vec<Scalar, 2> const l = points[idx];
                rtk::vec<Scalar, 2> const r = points[idx + 1];

                Scalar const ldx = r.x() - l.x();

                Scalar const Rxlx = R.x() * l.x();
                Scalar const Lxrx = L.x() * r.x();

                as[idx] = ldx / gdx;
                es[idx] = (Rxlx - Lxrx) / gdx;
            }
        }

        /// \brief Computes the `c` and `f` coefficients in-place.
        /// \param points Interpolation points of the FIF.
        /// \param ds Vertical scaling factors view.
        /// \param cs `c` coefficients store.
        /// \param fs `f` coefficients store.
        /// \throws std::runtime_error points.size() must be greater than 1.
        /// \throws std::runtime_error the sizes must respect points.size() - 1 == ds.size() == cs.size() == fs.size().
        static inline constexpr void compute_cfs(std::span<const rtk::vec<Scalar, 2>> const points, std::span<const Scalar> const ds, std::span<Scalar> const cs, std::span<Scalar> const fs) {
            if (points.size() <= 1) throw std::runtime_error("degenerate fif");
            if (ds.size() != points.size() - 1 || cs.size() != points.size() - 1 || fs.size() != points.size() - 1) throw std::runtime_error("inconsistent sizes");

            scion::usize const size = points.size() - 1;

            rtk::vec<Scalar, 2> const L = points.front();
            rtk::vec<Scalar, 2> const R = points.back();

            Scalar const gdx = R.x() - L.x();
            Scalar const gdy = R.y() - L.y();

            Scalar const RxLy = R.x() * L.y();
            Scalar const LxRy = L.x() * R.y();

            for (scion::usize idx = 0; idx < size; ++idx) {
                rtk::vec<Scalar, 2> const l = points[idx];
                rtk::vec<Scalar, 2> const r = points[idx + 1];

                 Scalar const ldy = r.y() - l.y();

                 Scalar const Rxly = R.x() * l.y();
                 Scalar const Lxry = L.x() * r.y();

                cs[idx] = ldy / gdx - ds[idx] * gdy / gdx;
                fs[idx] = (Rxly - Lxry) / gdx - ds[idx] * ((RxLy - LxRy) / gdx);
            }
        }

        /// \brief Computes the coefficients in-place.
        /// \param points Interpolation points of the FIF.
        /// \param ds Vertical scaling factors view.
        /// \param as `a` coefficients store.
        /// \param cs `c` coefficients store.
        /// \param es `e` coefficients store.
        /// \param fs `f` coefficients store.
        /// \throws std::runtime_error points.size() must be greater than 1.
        /// \throws std::runtime_error the sizes must respect points.size() - 1 == ds.size() == as.size() == cs.size() == es.size() == fs.size().
        static inline constexpr void compute_acefs(std::span<const rtk::vec<Scalar, 2>> const points, std::span<const Scalar> const ds, std::span<Scalar> const as, std::span<Scalar> const cs, std::span<Scalar> const es, std::span<Scalar> const fs) {
            if (points.size() <= 1) throw std::runtime_error("degenerate fif");
            if (ds.size() != points.size() - 1 || as.size() != points.size() - 1 || es.size() != points.size() - 1 || cs.size() != points.size() - 1 || fs.size() != points.size() - 1) throw std::runtime_error("inconsistent sizes");

            scion::usize const size = points.size() - 1;

            rtk::vec<Scalar, 2> const L = points.front();
            rtk::vec<Scalar, 2> const R = points.back();

            Scalar const gdx = R.x() - L.x();
            Scalar const gdy = R.y() - L.y();

            Scalar const RxLy = R.x() * L.y();
            Scalar const LxRy = L.x() * R.y();

            for (scion::usize idx = 0; idx < size; ++idx) {
                rtk::vec<Scalar, 2> const l = points[idx];
                rtk::vec<Scalar, 2> const r = points[idx + 1];

                Scalar const ldx = r.x() - l.x();
                Scalar const ldy = r.y() - l.y();

                Scalar const Rxlx = R.x() * l.x();
                Scalar const Lxrx = L.x() * r.x();
                Scalar const Rxly = R.x() * l.y();
                Scalar const Lxry = L.x() * r.y();

                as[idx] = ldx / gdx;
                es[idx] = (Rxlx - Lxrx) / gdx;
                cs[idx] = ldy / gdx - ds[idx] * gdy / gdx;
                fs[idx] = (Rxly - Lxry) / gdx - ds[idx] * ((RxLy - LxRy) / gdx);
            }
        }

        /// \brief Computes the mu of the trans.
        /// \param a `a` coefficient.
        /// \param c `c` coefficient.
        /// \param d `d` coefficient.
        /// \return mu of the trans.
        [[nodiscard]] static inline constexpr auto compute_mu(Scalar const a, Scalar const c, Scalar const d) -> rtk::vec<Scalar, 2> {
            rtk::fif_impl<Scalar>::check_acd(a, c, d);

            if (a == d) {
                return rtk::vec<Scalar, 2>{0.0, 1.0};
            }

            if (c == static_cast<Scalar>(0.0)) {
                return rtk::vec<Scalar, 2>{1.0, 0.0};
            }

            return rtk::vec<Scalar, 2>{a - d, c};
        }

        /// \brief Computes the nu of the trans.
        /// \param a `a` coefficient
        /// \param c `c` coefficient
        /// \param d `d` coefficient
        /// \return nu of the trans
        [[nodiscard]] static inline constexpr auto compute_nu(Scalar const a, Scalar const c, Scalar const d) -> rtk::vec<Scalar, 2> {
            rtk::fif_impl<Scalar>::check_acd(a, c, d);

            return rtk::vec<Scalar, 2>{0.0, 1.0};
        }

        /// \brief Computes the tangent along which points converge when trans is applied iteratively.
        /// \param a `a` coefficient.
        /// \param c `c` coefficient.
        /// \param d `d` coefficient.
        /// \return Tangent along which points converge when trans is applied iteratively.
        [[nodiscard]] static constexpr auto compute_tangent(Scalar const a, Scalar const c, Scalar const d) -> std::optional<rtk::vec<Scalar, 2>> {
            rtk::fif_impl<Scalar>::check_acd(a, c, d);

            // collinear case: no tangent
            if (a == d && c == static_cast<Scalar>(0.0)) return std::nullopt;

            // alternating case: no tangent
            if (a == -d) return std::nullopt;

            // mu-dominant case
            if (a > std::abs(d)) return rtk::fif_impl<Scalar>::compute_mu(a, c, d);

            // nu-dominant case
            if (a <= std::abs(d)) return rtk::fif_impl<Scalar>::compute_nu(a, c, d);

            throw std::logic_error("unexpected case"); // should be unreachable
        }

        /// \brief Computes the asymptotic direction along which the starting point will converge when trans is applied iteratively.
        /// \param a `a` coefficient.
        /// \param c `c` coefficient.
        /// \param d `d` coefficient.
        /// \return Asymptotic direction along which points converge when trans is applied iteratively.
        ///
        /// If the starting point is already at the fixed point, `std::nullopt` is returned.
        /// If the tangent exists, the asymptotic direction might still be different if the starting point is in already in an eigenspace.
        [[nodiscard]] static constexpr auto compute_asymptotic_direction(Scalar const a, Scalar const c, Scalar const d, Scalar const e, Scalar const f, rtk::vec<Scalar, 2> const starting_point) -> std::optional<std::variant<rtk::vec<Scalar, 2>, std::pair<rtk::vec<Scalar, 2>, rtk::vec<Scalar, 2>>>> {
            rtk::fif_impl<Scalar>::check_acdef(a, c, d, e, f);
            scion::assert_debug(rtk::is_vec_valid(starting_point), "invalid starting point");

            Scalar const dx = starting_point.x() - e;
            Scalar const dy = starting_point.y() - f;

            rtk::vec<Scalar, 2> const asym_e = rtk::vec<Scalar,2>{dx, dy};
            rtk::vec<Scalar, 2> const asym_o = rtk::vec<Scalar,2>{a * dx, c * dx + d * dy};

            // degenerate case: starting point is already at the fixed point, no asymptotic direction
            if (dx == static_cast<Scalar>(0.0) && dy == static_cast<Scalar>(0.0)) return std::nullopt;

            // collinear case: no dominant direction, the iterates converge on a straight line to the fixed point
            if (a == d && c == static_cast<Scalar>(0.0)) return asym_e;

            // alternating case: no dominant direction, the iterates converge on alternating straight lines to the fixed point depending on their parity
            if (a == -d) return std::pair{asym_e, asym_o};

            // starting point is already on the nu-eigenspace
            if (dx == static_cast<Scalar>(0.0)) return rtk::fif_impl<Scalar>::compute_nu(a, c, d);

            // starting point is already on the mu-eigenspace
            if (c * dx + (d - a) * dy == static_cast<Scalar>(0.0)) return rtk::fif_impl<Scalar>::compute_mu(a, c, d);

            // mu-dominant case
            if (a > std::abs(d)) return rtk::fif_impl<Scalar>::compute_mu(a, c, d);

            // nu-dominant case
            if (a <= std::abs(d)) return rtk::fif_impl<Scalar>::compute_nu(a, c, d);

            throw std::logic_error("unexpected case"); // should be unreachable
        }

        /// \brief Computes the left pseudo-tangent of the interpolation point at the specified index.
        /// \param as `a` coefficients view.
        /// \param cs `c` coefficients view.
        /// \param ds `d` coefficients view.
        /// \param index Index of the interpolation point.
        /// \return Left pseudo-tangent of the interpolation point at the specified index, if it exists.
        /// \throws std::runtime_error the sizes must respect 1 <= as.size() == cs.size() == ds.size().
        /// \throws std::runtime_error `index` must be in [1; min(as.size(), cs.size(), ds.size())].
        static constexpr auto compute_pseudo_tangent_l(std::span<Scalar const> const as, std::span<Scalar const> const cs, std::span<Scalar const> const ds, scion::usize const index) -> std::optional<rtk::vec<Scalar, 2>> {
            if (ds.size() != as.size() || ds.size() != cs.size()) throw std::runtime_error("inconsistent sizes");
            if (ds.size() == 0) throw std::runtime_error("degenerate fif");
            if (index == 0 || index > ds.size()) throw std::runtime_error("index out of bounds");

            scion::usize const idx_e = ds.size() - 1;
            scion::usize const idx_c = index - 1;

            return rtk::fif_impl<Scalar>::compute_tangent(as[idx_e], cs[idx_e], ds[idx_e]).transform([idx_c, &as, &cs, &ds](rtk::vec<Scalar, 2> const& tangent) -> rtk::vec<Scalar, 2> {
                return rtk::vec<Scalar, 2>{as[idx_c] * tangent.x(), cs[idx_c] * tangent.x() + ds[idx_c] * tangent.y() };
            });
        }

        /// \brief Computes the right pseudo-tangent of the interpolation point at the specified index.
        /// \param as `a` coefficients view.
        /// \param cs `c` coefficients view.
        /// \param ds `d` coefficients view.
        /// \param index Index of the interpolation point.
        /// \return Right pseudo-tangent of the interpolation point at the specified index, if it exists.
        /// \throws std::runtime_error the sizes must respect 1 <= as.size() == cs.size() == ds.size().
        /// \throws std::runtime_error `index` must be in [0; min(as.size(), cs.size(), ds.size()) - 1].
        static constexpr auto compute_pseudo_tangent_r(std::span<Scalar const> const as, std::span<Scalar const> const cs, std::span<Scalar const> const ds, scion::usize const index) -> std::optional<rtk::vec<Scalar, 2>> {
            if (ds.size() != as.size() || ds.size() != cs.size()) throw std::runtime_error("inconsistent sizes");
            if (ds.size() == 0) throw std::runtime_error("degenerate fif");
            if (index >= ds.size()) throw std::runtime_error("index out of bounds");

            scion::usize const idx_e = 0;
            scion::usize const idx_c = index;

            return rtk::fif_impl<Scalar>::compute_tangent(as[idx_e], cs[idx_e], ds[idx_e]).transform([idx_c, &as, &cs, &ds](rtk::vec<Scalar, 2> const& tangent) -> rtk::vec<Scalar, 2> {
                return rtk::vec<Scalar, 2>{as[idx_c] * tangent.x(), cs[idx_c] * tangent.x() + ds[idx_c] * tangent.y() };
            });
        };

        /// \brief Computes the left eigensector at the interpolation point at the specified index.
        /// \param as `a` coefficients view.
        /// \param cs `c` coefficients view.
        /// \param ds `d` coefficients view.
        /// \param index Index of the interpolation point.
        /// \return Left eigensector at the interpolation point at the specified index, if it exists.
        /// \throws std::runtime_error the sizes must respect 1 <= as.size() == cs.size() == ds.size().
        /// \throws std::runtime_error `index` must be in [1; min(as.size(), cs.size(), ds.size())].
        static constexpr auto compute_eigensector_l(std::span<const rtk::vec<Scalar, 2>> const points, std::span<const Scalar> const as, std::span<const Scalar> const cs, std::span<const Scalar> const ds, scion::usize const index) -> std::expected<std::pair<Scalar, Scalar>, std::monostate> {
            if (ds.size() != points.size() - 1 || as.size() != points.size() - 1 || cs.size() != points.size() - 1) throw std::runtime_error("inconsistent sizes");
            if (points.size() <= 1) throw std::runtime_error("degenerate fif");
            if (index == 0 || index > ds.size()) throw std::runtime_error("index out of bounds");

            scion::usize const idx_e = ds.size() - 1;
            scion::usize const idx_c = index - 1;

            Scalar const ae = as[idx_e];
            Scalar const ce = cs[idx_e];
            Scalar const de = ds[idx_e];

            Scalar const ac = as[idx_c];
            Scalar const cc = cs[idx_c];
            Scalar const dc = ds[idx_c];

            if (idx_c == idx_e) {
                if (ce > 0.0f) return std::make_pair(std::atan(ce / 2 * ae), std::numbers::pi_v<Scalar> / 2.0f);
                if (ce < 0.0f) return std::make_pair(-std::numbers::pi_v<Scalar> / 2.0f, std::atan(ce / 2 * ae));

                return std::unexpected(std::monostate{});
            } else {
                if (std::abs(de) >= ae) return std::unexpected(std::monostate{});

                return std::make_pair(std::atan((cc * (ae - de) - std::abs(ce)) / (ac * (ae - de))), std::atan((cc * (ae - de) + std::abs(ce)) / (ac * (ae - de))));
            }
        }

        /// \brief Computes the right eigensector at the interpolation point at the specified index.
        /// \param as `a` coefficients view.
        /// \param cs `c` coefficients view.
        /// \param ds `d` coefficients view.
        /// \param index index of the interpolation point.
        /// \return right eigensector at the interpolation point.
        /// \throws std::runtime_error the sizes must respect 1 <= as.size() == cs.size() == ds.size().
        /// \throws std::runtime_error `index` must be in [0; min(as.size(), cs.size(), ds.size()) - 1].
        static constexpr auto compute_eigensector_r(std::span<const rtk::vec<Scalar, 2>> const points, std::span<const Scalar> const as, std::span<const Scalar> const cs, std::span<const Scalar> const ds, scion::usize const index) -> std::expected<std::pair<Scalar, Scalar>, std::monostate> {
            if (ds.size() != points.size() - 1 || as.size() != points.size() - 1 || cs.size() != points.size() - 1) throw std::runtime_error("inconsistent sizes");
            if (points.size() <= 1) throw std::runtime_error("degenerate fif");
            if (index >= ds.size()) throw std::runtime_error("index out of bounds");

            scion::usize const idx_e = 0;
            scion::usize const idx_c = index;

            Scalar const ae = as[idx_e];
            Scalar const ce = cs[idx_e];
            Scalar const de = ds[idx_e];

            Scalar const ac = as[idx_c];
            Scalar const cc = cs[idx_c];
            Scalar const dc = ds[idx_c];

            if (idx_c == idx_e) {
                if (ce > 0.0f) return std::make_pair(std::atan(ce / 2 * ae), std::numbers::pi_v<Scalar> / 2.0f);
                if (ce < 0.0f) return std::make_pair(-std::numbers::pi_v<Scalar> / 2.0f, std::atan(ce / 2 * ae));

                return std::unexpected(std::monostate{});
            } else {
                if (std::abs(de) >= ae) return std::unexpected(std::monostate{});

                return std::make_pair(std::atan((cc * (ae - de) - std::abs(ce)) / (ac * (ae - de))), std::atan((cc * (ae - de) + std::abs(ce)) / (ac * (ae - de))));
            }
        }

        /// \brief Checks if the trans is a contraction.
        /// \param a `a` coefficient
        /// \param c `c` coefficient
        /// \param d `d` coefficient
        /// \return Checks if the trans is a contraction.
        static inline constexpr auto check_contraction(Scalar const a, Scalar const c, Scalar const d) -> bool {
            rtk::fif_impl<Scalar>::check_acd(a, c, d);

            Scalar const trace = a * a + c * c + d * d;
            Scalar const det = (a * d) * (a * d);
            Scalar const disc = trace * trace - 4 * det;

            Scalar const lambda_max = (trace + std::sqrt(disc)) / 2;

            return lambda_max < static_cast<Scalar>(1.0);
        }

        /// \brief checks if the trans at the specified index is a contraction mapping
        /// \param as `a` coefficients view
        /// \param cs `c` coefficients view
        /// \param ds `d` coefficients view
        /// \throws std::runtime_error the sizes must respect 1 <= as.size() == cs.size() == ds.size()
        static inline constexpr auto check_contractions(std::span<const Scalar> const as, std::span<const Scalar> const cs, std::span<const Scalar> const ds) -> bool {
            for (scion::usize idx = 0; idx < as.size(); ++idx) if (!rtk::fif_impl<Scalar>::check_contraction(as[idx], cs[idx], ds[idx])) return false;

            return true;
        }

        /// \brief makes a fif with the specified interpolation points
        /// \tparam Ptr std::unique_ptr or std::shared_ptr if you want a dynamically allocated fif
        /// \param points interpolation points to use
        /// \return newly created fif, or error
        template <template<typename> typename Ptr = std::type_identity>
            requires (std::is_same_v<Ptr<rtk::fif_impl<Scalar>>, std::type_identity<rtk::fif_impl<Scalar>>> || std::is_same_v<Ptr<rtk::fif_impl<Scalar>>, std::unique_ptr<rtk::fif_impl<Scalar>>> || std::is_same_v<Ptr<rtk::fif_impl<Scalar>>, std::shared_ptr<rtk::fif_impl<Scalar>>>)
        static auto make(std::vector<rtk::vec<Scalar, 2>> && points) -> std::expected<std::conditional_t<std::is_same_v<Ptr<rtk::fif_impl<Scalar>>, std::type_identity<rtk::fif_impl<Scalar>>>, rtk::fif_impl<Scalar>, Ptr<rtk::fif_impl<Scalar>>>, std::invalid_argument> {
            std::vector<Scalar> as(points.size() - 1);
            std::vector<Scalar> es(points.size() - 1);

            std::vector<Scalar> cs(points.size() - 1);
            std::vector<Scalar> fs(points.size() - 1);

            rtk::fif_impl<Scalar>::compute_aes(points, as, es);
            std::vector<Scalar> ds = as;
            rtk::fif_impl<Scalar>::compute_cfs(points, ds, cs, fs);

            if constexpr (std::is_same_v<Ptr<rtk::fif_impl<Scalar>>, std::type_identity<rtk::fif_impl<Scalar>>>) {
                return rtk::fif_impl<Scalar>{
                    {
                        std::move(points),
                        std::move(as),
                        std::move(cs),
                        std::move(ds),
                        std::move(es),
                        std::move(fs),
                        std::vector<rtk::vec<Scalar, 2>>{},
                        true,
                    }
                };
            } else {
                return Ptr<rtk::fif_impl<Scalar>>(new rtk::fif_impl<Scalar>{
                    {
                        std::move(points),
                        std::move(as),
                        std::move(cs),
                        std::move(ds),
                        std::move(es),
                        std::move(fs),
                        std::vector<rtk::vec<Scalar, 2>>{},
                        true,
                    }
                });
            }
        }

    private: // static private functions

    protected: // static protected functions

    public: // static public functions

    private: // private members

    protected: // protected members
        /// \brief member data
        struct M {
            std::vector<rtk::vec<Scalar, 2>> points;

            std::vector<Scalar> as;
            std::vector<Scalar> cs;
            std::vector<Scalar> ds;
            std::vector<Scalar> es;
            std::vector<Scalar> fs;

            std::vector<rtk::vec<Scalar, 2>> attractor{};

            scion::usize iterations = 4;
        } m;

    public: // public members

    private: // private methods

    protected: // protected methods
        explicit fif_impl(rtk::fif_impl<Scalar>::M m) noexcept : m{std::move(m)} {
            rtk::log::trace("init", "ctor for rtk::fif@{}", static_cast<void *>(this));
        }

    public: // public methods
        /// \brief copy constructor
        /// \param other other model to copy
        fif_impl(rtk::fif_impl<Scalar> const & other): m{other.m} {
            rtk::log::trace("copy", "ctor for rtk::fif@{} (from rtk::fif@{})", static_cast<void *>(this), static_cast<void const *>(&other));
        }

        /// \brief move constructor
        /// \param other other model to move
        fif_impl(rtk::fif_impl<Scalar> && other) noexcept: m{std::move(other.m)} {
            rtk::log::trace("move", "ctor for rtk::fif@{} (from rtk::fif@{})", static_cast<void *>(this), static_cast<void *>(&other));
        }

        /// \brief copy assignment operator
        /// \param other other model to copy
        /// \return instance to which it was copied
        rtk::fif_impl<Scalar> & operator=(rtk::fif_impl<Scalar> const & other) {
            this->m = other.m;

            rtk::log::trace("copy", "otor for rtk::fif@{} (from rtk::fif@{})", static_cast<void *>(this), static_cast<void const *>(&other));

            return *this;
        }

        /// \brief move assignment operator
        /// \param other other model to move
        /// \return instance to which it was moved
        rtk::fif_impl<Scalar> & operator=(rtk::fif_impl<Scalar> && other) noexcept {
            this->m = std::move(other.m);

            rtk::log::trace("move", "otor for rtk::fif@{} (from rtk::fif@{})", static_cast<void *>(this), static_cast<void *>(&other));

            return *this;
        }

        /// \brief destructor for a fif
        ~fif_impl() noexcept {
            rtk::log::trace("drop", "dtor for rtk::fif@{}", static_cast<void *>(this));
        }

        /// \brief moves an interpolation point
        /// \param index index of the point to move
        /// \param point new coordinates of the point
        /// \return new index of the point
        auto move_point(scion::usize const index, rtk::vec<Scalar, 2> const& point) -> std::expected<scion::usize, fif_error> {
            if (index >= m.points.size()) throw std::runtime_error("point index out of bounds");

            // If there's only one point, update it directly and exit early.
            if (m.points.size() == 1) {
                m.points[0] = point;

                return 0;
            }

            bool const extremity_old = index == 0 || index == m.points.size() - 1;
            bool const extremity_new = index == 0 || index == m.points.size() - 1;

            bool const extremity_change = extremity_old || extremity_new;

            m.points.erase(m.points.begin() + static_cast<scion::isize>(index));
            auto const it = std::ranges::lower_bound(m.points, point, [](rtk::vec<Scalar, 2> const& a, rtk::vec<Scalar, 2> const& b) -> bool {
                return a.x() != b.x() ? a.x() < b.x() : a.y() < b.y();
            });
            m.points.insert(it, point);

            rtk::fif_impl<Scalar>::compute_acefs(m.points, m.ds, m.as, m.cs, m.es, m.fs);

            if (!rtk::fif_impl<Scalar>::check_contractions(m.as, m.cs, m.ds)) {
                for (scion::usize idx = 0; idx < m.points.size() - 1; ++idx) if (!rtk::fif_impl<Scalar>::check_contraction(m.as[idx], m.cs[idx], m.ds[idx])) rtk::log::error("error", "fif has a transformation that is not a contraction mapping at index {})\n", idx);

                return rtk::fif_error::not_contractive;
            }

            return std::distance(m.points.begin(), it);
        }

        void compute_attractor_recurse(scion::usize const iterations, rtk::vec<Scalar, 2> const& pa, rtk::vec<Scalar, 2> const& pb) {
            if (iterations == 0) {
                m.attractor.emplace_back(pa);
                m.attractor.emplace_back(pb);

                return;
            }

            for (scion::usize i = 0; i < m.points.size() - 1; ++i) {
                Scalar const a = m.as[i];
                Scalar const e = m.es[i];
                Scalar const d = m.ds[i];
                Scalar const c = m.cs[i];
                Scalar const f = m.fs[i];

                this->compute_attractor_recurse(iterations - 1, {a * pa.x() + e, c * pa.x() + d * pa.y() + f}, {a * pb.x() + e, c * pb.x() + d * pb.y() + f});
            }
        }

        void compute_attractor(scion::usize const iterations) {
            if (m.points.size() < 2) {
                return;
            }

            m.attractor.clear();

            compute_attractor_recurse(iterations, m.points.front(), m.points.back());
        }

        auto attractor() -> std::span<const rtk::vec<Scalar, 2>> {
            return m.attractor;
        }

        /// \brief returns the number of interpolation points
        /// \return number of interpolation points
        [[nodiscard]] constexpr inline auto size_points() const noexcept -> scion::usize {
            return m.points.size();
        }

        /// \brief returns the number of transformations
        /// \return number of transformations
        [[nodiscard]] constexpr inline auto size_trans() const noexcept -> scion::usize {
            return m.as.size();
        }

        /// \brief returns the interpolation point at the specified index
        /// \param index index of the interpolation point
        /// \return interpolation point
        /// \throws std::runtime_error index must be [0; this->size_points() - 1]
        /// \note if you want a nothrow alternative, use the subset operator on the result of this->points()
        [[nodiscard]] constexpr inline auto point(scion::usize const index) const -> rtk::vec<Scalar, 2> const& {
            if (index >= m.points.size()) throw std::runtime_error("index out of bounds");

            return m.points[index];
        }

        /// \brief returns the interpolation points
        /// \return interpolation points
        [[nodiscard]] constexpr inline auto points() const noexcept -> std::span<const rtk::vec<Scalar, 2>> {
            return m.points;
        }

        /// \brief returns the `a` coefficient at the specified index
        /// \param index index of the coefficient
        /// \return `a` coefficient
        /// \throws std::runtime_error index must be [0; this->size_trans() - 1]
        /// \note if you want a nothrow alternative, use the subset operator on the result of this->as()
        [[nodiscard]] constexpr inline auto a(scion::usize const index) const -> Scalar {
            if (index >= m.as.size()) throw std::runtime_error("index out of bounds");

            return m.as[index];
        }

        /// \brief returns the `c` coefficient at the specified index
        /// \param index index of the coefficient
        /// \return `c` coefficient
        /// \throws std::runtime_error index must be [0; this->size_trans() - 1]
        /// \note if you want a nothrow alternative, use the subset operator on the result of this->cs()
        [[nodiscard]] constexpr inline auto c(scion::usize const index) const -> Scalar {
            if (index >= m.cs.size()) throw std::runtime_error("index out of bounds");

            return m.cs[index];
        }

        /// \brief returns the `d` coefficient at the specified index
        /// \param index index of the coefficient
        /// \return `d` coefficient
        /// \throws std::runtime_error index must be [0; this->size_trans() - 1]
        /// \note if you want a nothrow alternative, use the subset operator on the result of this->ds()
        [[nodiscard]] constexpr inline auto d(scion::usize const index) const -> Scalar {
            if (index >= m.ds.size()) throw std::runtime_error("index out of bounds");

            return m.ds[index];
        }

        /// \brief returns the `e` coefficient at the specified index
        /// \param index index of the coefficient
        /// \return `e` coefficient
        /// \throws std::runtime_error index must be [0; this->size_trans() - 1]
        /// \note if you want a nothrow alternative, use the subset operator on the result of this->es()
        [[nodiscard]] constexpr inline auto e(scion::usize const index) const -> Scalar {
            if (index >= m.es.size()) throw std::runtime_error("index out of bounds");

            return m.es[index];
        }

        /// \brief returns the `f` coefficient at the specified index
        /// \param index index of the coefficient
        /// \return `f` coefficient
        /// \throws std::runtime_error index must be [0; this->size_trans() - 1]
        /// \note if you want a nothrow alternative, use the subset operator on the result of this->fs()
        [[nodiscard]] constexpr inline auto f(scion::usize const index) const -> Scalar {
            if (index >= m.fs.size()) throw std::runtime_error("index out of bounds");

            return m.fs[index];
        }

        /// \brief Returns the `a` coefficients.
        /// \return `a` coefficients.
        [[nodiscard]] constexpr inline auto as() const noexcept -> std::span<const Scalar> {
            return m.as;
        }

        /// \brief Returns the `c` coefficients.
        /// \return `c` coefficients.
        [[nodiscard]] constexpr inline auto cs() const noexcept -> std::span<const Scalar> {
            return m.cs;
        }

        /// \brief Returns the `d` coefficients.
        /// \return `d` coefficients.
        [[nodiscard]] constexpr inline auto ds() const noexcept -> std::span<const Scalar> {
            return m.ds;
        }

        /// \brief Returns the `e` coefficients.
        /// \return `e` coefficients.
        [[nodiscard]] constexpr inline auto es() const noexcept -> std::span<const Scalar> {
            return m.es;
        }

        /// \brief Returns the `f` coefficients.
        /// \return `f` coefficients.
        [[nodiscard]] constexpr inline auto fs() const noexcept -> std::span<const Scalar> {
            return m.fs;
        }

        /// \see rtk::fif_impl::compute_mu
        [[nodiscard]] auto mu(scion::usize const index) const -> rtk::vec<Scalar, 2> {
            if (index >= m.as.size()) throw std::runtime_error("index out of bounds");

            return rtk::fif_impl<Scalar>::compute_mu(m.as[index], m.cs[index], m.ds[index]);
        }

        /// \see rtk::fif_impl::compute_nu
        [[nodiscard]] auto nu([[maybe_unused]] scion::usize const index) const -> rtk::vec<Scalar, 2> {
            if (index >= m.as.size()) throw std::runtime_error("index out of bounds");

            return rtk::fif_impl<Scalar>::compute_nu(m.as[index], m.cs[index], m.ds[index]);
        }

        /// \see rtk::fif_impl::compute_tangent
        [[nodiscard]] auto tangent(scion::usize const index) const -> std::optional<rtk::vec<Scalar, 2>> {
            if (index >= m.as.size()) throw std::runtime_error("index out of bounds");

            return rtk::fif_impl<Scalar>::compute_tangent(m.as[index], m.cs[index], m.ds[index]);
        }

        /// \see rtk::fif_impl::compute_asymptotic_direction
        [[nodiscard]] auto asymptotic_direction(scion::usize const index, rtk::vec<Scalar, 2> const starting_point) -> std::optional<std::variant<rtk::vec<Scalar, 2>, std::pair<rtk::vec<Scalar, 2>, rtk::vec<Scalar, 2>>>> {
            if (index >= m.as.size()) throw std::runtime_error("index out of bounds");

            return rtk::fif_impl<Scalar>::compute_asymptotic_direction(m.as[index], m.cs[index], m.ds[index], m.es[index], m.fs[index], starting_point);
        }

        /// \see rtk::fif_impl::compute_pseudo_tangent_r
        [[nodiscard]] auto pseudo_tangent_r(scion::usize const index) const -> std::optional<rtk::vec<Scalar, 2>> {
            return rtk::fif_impl<Scalar>::compute_pseudo_tangent_r(m.as, m.cs, m.ds, index);
        }

        /// \see rtk::fif_impl::compute_pseudo_tangent_l
        [[nodiscard]] auto pseudo_tangent_l(scion::usize const index) const -> std::optional<rtk::vec<Scalar, 2>> {
            return rtk::fif_impl<Scalar>::compute_pseudo_tangent_l(m.as, m.cs, m.ds, index);
        }

        /// \brief sets the vertical scaling factor of the trans at the specified index to the given value
        /// \param index index of the trans
        /// \param d new vertical scaling factor
        /// \throws std::runtime_error d must be strictly in [-1; 1]
        void set_di(scion::usize const index, Scalar const d) {
            if (std::norm(d) >= static_cast<Scalar>(1.0)) {
                throw std::runtime_error("not a contraction");
            }

            m.ds[index] = d;

            rtk::fif_impl<Scalar>::compute_acefs(m.points, m.ds, m.as, m.cs, m.es, m.fs);

            if (!rtk::fif_impl<Scalar>::check_contractions(m.as, m.cs, m.ds)) {
                for (scion::usize idx = 0; idx < m.points.size() - 1; ++idx) if (!rtk::fif_impl<Scalar>::check_contraction(m.as[idx], m.cs[idx], m.ds[idx])) rtk::log::error("error", "fif has a transformation that is not a contraction mapping at index {})\n", idx);

                return;
            }
        }

        /// \brief sets the vertical scaling factor of the trans at the specified index so that the angle of the tangent is the specified angle, if possible, and returns it
        /// \param index index of the trans
        /// \param angle desired angle
        /// \return value given to the vertical scaling factor
        auto set_angle(scion::usize const index, Scalar const angle) -> std::optional<Scalar> {
            Scalar const dx = m.points[index + 1].x() - m.points[index].x();
            Scalar const dX = m.points.back().x() - m.points.front().x();
            Scalar const dy = m.points[index + 1].y() - m.points[index].y();
            Scalar const dY = m.points.back().y() - m.points.front().y();

            Scalar const dp = (dx + dy * std::tan(angle)) / (dX + dY * std::tan(angle));
            Scalar const dm = (dx - dy * std::tan(angle)) / (dX - dY * std::tan(angle));

            if (std::norm(dp) <= m.as[index]) {
                m.ds[index] = dp;
            } else if (std::norm(dm) <= m.as[index]) {
                m.ds[index] = dm;
            } else {
                return std::nullopt;
            }

            rtk::fif_impl<Scalar>::compute_acefs(m.points, m.ds, m.as, m.cs, m.es, m.fs);

            if (!rtk::fif_impl<Scalar>::check_contractions(m.as, m.cs, m.ds)) {
                for (scion::usize idx = 0; idx < m.points.size() - 1; ++idx) if (!rtk::fif_impl<Scalar>::check_contraction(m.as[idx], m.cs[idx], m.ds[idx])) rtk::log::error("error", "fif has a transformation that is not a contraction mapping at index {})\n", idx);

                return std::nullopt;
            }

            return m.ds[index];
        }

        /// \see rtk::fif_impl::compute_eigensector_r
        auto eigensector_r(scion::usize const index) const -> std::expected<std::pair<Scalar, Scalar>, std::monostate> {
            return rtk::fif_impl<Scalar>::compute_eigensector_r(m.points, m.as, m.cs, m.ds, index);
        }

        /// \see rtk::fif_impl::compute_eigensector_l
        auto eigensector_l(scion::usize const index) const -> std::expected<std::pair<Scalar, Scalar>, std::monostate> {
            return rtk::fif_impl<Scalar>::compute_eigensector_l(m.points, m.as, m.cs, m.ds, index);
        }
    };

    using fif = rtk::fif_impl<>;
}
