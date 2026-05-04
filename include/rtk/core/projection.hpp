#pragma once

#include "scion/core/core.hpp"
#include "scion/util/util.hpp"

#include "rtk/core/types.hpp"

#include "rtk/util/log.hpp"

#include <expected>
#include <memory>
#include <ranges>
#include <variant>
#include <vector>

/// \file rtk/core/projection.hpp
/// \brief provides a class representing a projection for free-form deformations

namespace rtk
{
    /// \class rtk::projection
    /// \brief represents a projection for free-form deformations
    ///
    /// \invariant the points of the set are always valid (no NaN or infinite values)
    /// \invariant the points of the set are always of the same width
    ///
    /// If client code obtain a mutable reference or pointer to internal data of the class, the client becomes responsible for keeping up the invariants.
    ///
    /// Many function or methods provide a \c Dim template value parameter.
    /// Using it allows usage of 2D and 3D points (rtk::vec2f or rtk::vec3f).
    /// By default, points are dynamically-sized.
    class projection final {
    private: // static private members

    protected: // static protected members

    public: // static public members
        /// \brief makes a new projection with canonical name, the given dimension and width, and zero-initializes all points
        /// \tparam Dim dimension of the points (2 or 3 for fixed, 0 for dynamic)
        /// \param dimension dimension of the projection
        /// \param width width of the projection
        /// \return newly created projection on success, error description on failure
        /// \errors dimension must be at least 2
        /// \errors dimension must be equal to the template parameter Dim if it is non-zero
        /// \errors width must be at least 2
        template <scion::usize Dim = 0, template<typename> typename Ptr = std::type_identity>
            requires ((Dim == 0 || Dim == 2 || Dim == 3) && (std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>> || std::is_same_v<Ptr<rtk::projection>, std::unique_ptr<rtk::projection>> || std::is_same_v<Ptr<rtk::projection>, std::shared_ptr<rtk::projection>>))
        [[nodiscard]] static auto make(scion::usize const dimension, scion::usize const width) -> std::expected<std::conditional_t<std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>>, rtk::projection, Ptr<rtk::projection>>, std::invalid_argument> {
            if (dimension < 2) {
                return std::unexpected{std::invalid_argument{"dimension must be at least 2"}};
            }

            if (width < 2) {
                return std::unexpected{std::invalid_argument{"width must be at least 2"}};
            }

            if constexpr (Dim != 0) {
                if (dimension != Dim) {
                    return std::unexpected{std::invalid_argument{"dimension is inconsistent"}};
                }
            }

            if constexpr (std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>>) {
                return rtk::projection{
                    {
                        std::format("r{}w{}", dimension, width),
                        true,
                        std::vector<rtk::vecf<Dim>>{width, rtk::vecf<Dim>::Zero(static_cast<scion::isize>(dimension))},
                        0,
                    }
                };
            } else {
                return Ptr<rtk::projection>(new rtk::projection{
                    {
                        std::format("r{}w{}", dimension, width),
                        true,
                        std::vector<rtk::vecf<Dim>>{width, rtk::vecf<Dim>::Zero(static_cast<scion::isize>(dimension))},
                        0,
                    }
                });
            }
        }

        /// \brief makes a new projection with canonical name and the given points
        /// \tparam Dim dimension of the points (2 or 3 for fixed, 0 for dynamic)
        /// \param points points to use for the projection
        /// \return newly created projection on success, error description on failure
        /// \errors dimension must be at least 2
        /// \errors dimension must be equal to the template parameter Dim if it is non-zero
        /// \errors width must be at least 2
        /// \errors points dimension inconsistent
        template <scion::usize Dim = 0, template<typename> typename Ptr = std::type_identity>
            requires ((Dim == 0 || Dim == 2 || Dim == 3) && (std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>> || std::is_same_v<Ptr<rtk::projection>, std::unique_ptr<rtk::projection>> || std::is_same_v<Ptr<rtk::projection>, std::shared_ptr<rtk::projection>>))
        [[nodiscard]] static auto make(std::vector<vecf<Dim>> && points) -> std::expected<std::conditional_t<std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>>, rtk::projection, Ptr<rtk::projection>>, std::invalid_argument> {
            if (points.size() < 2) {
                return std::unexpected{std::invalid_argument{"width must be at least 2"}};
            }

            if (points[0].size() < 2) {
                return std::unexpected{std::invalid_argument{"dimension must be at least 2"}};
            }

            if (scion::usize const width = points[0].size(); std::ranges::any_of(points, [width](rtk::vecxf const & point) {
                return point.size() != width;
            })) {
                return std::unexpected{std::invalid_argument{"points dimension inconsistent"}};
            }

            if constexpr (Dim != 0) {
                if (points[0].size() != Dim) {
                    return std::unexpected{std::invalid_argument{"dimension is inconsistent"}};
                }
            }

            if constexpr (std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>>) {
                return rtk::projection{
                    {
                        std::format("r{}w{}", points[0].size(), points.size()),
                        true,
                        std::move(points),
                        0,
                    }
                };
            } else {
                return Ptr<rtk::projection>(new rtk::projection{
                    {
                        std::format("r{}w{}", points[0].size(), points.size()),
                        true,
                        std::move(points),
                        0,
                    }
                });
            }
        }

        /// \brief makes a new projection with the given name, dimension and width, and zero-initializes all points
        /// \tparam Dim dimension of the points (2 or 3 for fixed, 0 for dynamic)
        /// \param name name of the projection
        /// \param dimension dimension of the projection
        /// \param width width of the projection
        /// \return newly created projection on success, error description on failure
        /// \errors dimension must be at least 2
        /// \errors dimension must be equal to the template parameter Dim if it is non-zero
        /// \errors width must be at least 2
        /// \errors name is not a non-empty string of alphanumeric characters, underscores or dashes
        template <scion::usize Dim = 0, template<typename> typename Ptr = std::type_identity>
            requires ((Dim == 0 || Dim == 2 || Dim == 3) && (std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>> || std::is_same_v<Ptr<rtk::projection>, std::unique_ptr<rtk::projection>> || std::is_same_v<Ptr<rtk::projection>, std::shared_ptr<rtk::projection>>))
        [[nodiscard]] static auto make(std::string && name, scion::usize const dimension, scion::usize const width) -> std::expected<std::conditional_t<std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>>, rtk::projection, Ptr<rtk::projection>>, std::invalid_argument> {
            if (dimension < 2) {
                return std::unexpected{std::invalid_argument{"dimension must be at least 2"}};
            }

            if (width < 2) {
                return std::unexpected{std::invalid_argument{"width must be at least 2"}};
            }

            if (name.empty()) {
                return std::unexpected{std::invalid_argument{"name is empty"}};
            }

            if (std::ranges::any_of(name, [](char const c) -> bool {
                return !(std::isalnum(c) || c == '_' || c == '-');
            })) {
                return std::unexpected{std::invalid_argument{"name is not alphanumeric, dash or underscore"}};
            }

            if constexpr (Dim != 0) {
                if (dimension != Dim) {
                    return std::unexpected{std::invalid_argument{"dimension is inconsistent"}};
                }
            }

            if constexpr (std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>>) {
                return rtk::projection{
                    {
                        std::move(name),
                        false,
                        std::vector<rtk::vecf<Dim>>{width, rtk::vecf<Dim>::Zero(static_cast<scion::isize>(dimension))},
                        0,
                    }
                };
            } else {
                return Ptr<rtk::projection>(new rtk::projection{
                    {
                        std::move(name),
                        false,
                        std::vector<rtk::vecf<Dim>>{width, rtk::vecf<Dim>::Zero(static_cast<scion::isize>(dimension))},
                        0,
                    }
                });
            }
        }

        /// \brief makes a new projection with the given name and points
        /// \tparam Dim dimension of the points (2 or 3 for fixed, 0 for dynamic)
        /// \param name name of the projection
        /// \param points points to use for the projection
        /// \return newly created projection on success, error description on failure
        /// \errors dimension must be at least 2
        /// \errors dimension must be equal to the template parameter Dim if it is non-zero
        /// \errors width must be at least 2
        /// \errors points dimension inconsistent
        /// \errors name is not a non-empty string of alphanumeric characters, underscores or dashes
        template <scion::usize Dim = 0, template<typename> typename Ptr = std::type_identity>
            requires ((Dim == 0 || Dim == 2 || Dim == 3) && (std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>> || std::is_same_v<Ptr<rtk::projection>, std::unique_ptr<rtk::projection>> || std::is_same_v<Ptr<rtk::projection>, std::shared_ptr<rtk::projection>>))
        [[nodiscard]] static auto make(std::string && name, std::vector<rtk::vecf<Dim>> && points) -> std::expected<std::conditional_t<std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>>, rtk::projection, Ptr<rtk::projection>>, std::invalid_argument> {
            if (points.size() < 2) {
                return std::unexpected{std::invalid_argument{"width must be at least 2"}};
            }

            if (points[0].size() < 2) {
                return std::unexpected{std::invalid_argument{"dimension must be at least 2"}};
            }

            if (scion::usize const width = points[0].size(); std::ranges::any_of(points, [width](rtk::vecxf const & point) {
                return point.size() != width;
            })) {
                return std::unexpected{std::invalid_argument{"points dimension inconsistent"}};
            }

            if (name.empty()) {
                return std::unexpected{std::invalid_argument{"name is empty"}};
            }

            if (std::ranges::any_of(name, [](char const c) -> bool {
                return !(std::isalnum(c) || c == '_' || c == '-');
            })) {
                return std::unexpected{std::invalid_argument{"name is not alphanumeric, dash or underscore"}};
            }

            if constexpr (Dim != 0) {
                if (points[0].size() != Dim) {
                    return std::unexpected{std::invalid_argument{"dimension is inconsistent"}};
                }
            }

            if constexpr (std::is_same_v<Ptr<rtk::projection>, std::type_identity<rtk::projection>>) {
                return rtk::projection{
                    {
                        std::move(name),
                        false,
                        std::move(points),
                        0,
                    }
                };
            } else {
                return Ptr<rtk::projection>(new rtk::projection{
                    {
                        std::move(name),
                        false,
                        std::move(points),
                        0,
                    }
                });
            }
        }

    private: // static private functions

    protected: // static protected functions

    public: // static public functions

    private: // private members
        struct M {
            /// \brief stores the name of the projection
            std::string name;

            /// \brief Stores \c true if the name is canonical
            bool name_canonical;

            /// \brief stores the points used for the projection
            std::variant<std::vector<rtk::vecf<2>>, std::vector<rtk::vecf<3>>, std::vector<rtk::vecxf>> points;

            /// \brief is incremented every time the projection is modified
            ///
            /// Client code is able to poll this to check if an update is required on their part.
            scion::usize counter;
        } m;

    protected: // protected members

    public: // public members

    private: // private methods
        /// \brief Private constructor for a projection, specifically for R2.
        /// \param m member data
        explicit projection(rtk::projection::M && m) noexcept : m{std::move(m)} {
            rtk::log::trace("init", "ctor for rtk::projection@{}", static_cast<void *>(this));
        }

    protected: // protected methods

    public: // public methods
        /// \brief copy constructor for a projection
        /// \param other instance to copy
        projection(rtk::projection const & other): m{other.m} {
            rtk::log::trace("copy", "ctor for rtk::projection@{} (from rtk::projection@{})", static_cast<void *>(this), static_cast<void const *>(&other));
        }

        /// \brief move constructor for a projection
        /// \param other instance to move
        projection(rtk::projection && other) noexcept: m{std::move(other.m)} {
            rtk::log::trace("move", "ctor for rtk::projection@{} (from rtk::projection@{})", static_cast<void *>(this), static_cast<void *>(&other));
        }

        /// \brief copy-assignment-operator for a projection
        /// \param other instance to copy
        /// \return instance to which it was copied
        rtk::projection & operator=(rtk::projection const & other) {
            this->m = other.m;

            rtk::log::trace("copy", "otor for rtk::projection@{} (from rtk::projection@{})", static_cast<void *>(this), static_cast<void const *>(&other));

            return *this;
        }

        /// \brief move-assignment-operator for a projection
        /// \param other instance to move
        /// \return instance to which it was moved
        rtk::projection & operator=(rtk::projection && other) noexcept {
            this->m = std::move(other.m);

            rtk::log::trace("move", "otor for rtk::projection@{} (from rtk::projection@{})", static_cast<void *>(this), static_cast<void *>(&other));

            return *this;
        }

        /// \brief destructor for a projection
        ~projection() {
            rtk::log::trace("drop", "dtor for rtk::projection@{}", static_cast<void *>(this));
        }

        /// \brief returns a reference to the points of the projection
        /// \tparam Dim dimension of the points (2 or 3 for fixed, 0 for dynamic)
        /// \return points of the projection, if they are of the queried dimension
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        [[nodiscard]] auto points() const noexcept -> std::optional<std::span<rtk::vecf<Dim> const>> {
            if (!std::holds_alternative<std::vector<rtk::vecf<Dim>>>(m.points)) {
                return std::nullopt;
            }

            return std::span{std::get<std::vector<rtk::vecf<Dim>>>(m.points)};
        }

        /// \brief returns a reference to the points of the projection
        /// \tparam Dim dimension of the points (2 or 3 for fixed, 0 for dynamic)
        /// \return points of the projection, if they are of the queried dimension
        /// \note client code is responsible for restoring the invariants that are violated when the returned view is mutated
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        [[nodiscard]] auto points() noexcept -> std::optional<std::span<rtk::vecf<Dim>>> {
            if (!std::holds_alternative<std::vector<rtk::vecf<Dim>>>(m.points)) {
                return std::nullopt;
            }

            return std::span{std::get<std::vector<rtk::vecf<Dim>>>(m.points)};
        }

        /// \brief returns the width of the projection
        /// \return width of the projection
        /// \throws std::logic_error on DEBUG or PROFILE if there are points with different dimensions or if the projection has less than 2 points
        /// \note this value is not cached, and is especially slow with dynamically sized points: they must all be checked for consistency in DEBUG or PROFILE
        [[nodiscard]] auto dimension() const -> scion::usize {
            if (std::holds_alternative<std::vector<rtk::vec2f>>(m.points)) {
                if constexpr (scion::ENABLE_ASSERTIONS) {
                    scion::assert_debug(std::get<std::vector<rtk::vec2f>>(m.points).size() >= 2);
                }

                return 2;
            }

            if (std::holds_alternative<std::vector<rtk::vec3f>>(m.points)) {
                if constexpr (scion::ENABLE_ASSERTIONS) {
                    scion::assert_debug(std::get<std::vector<rtk::vec3f>>(m.points).size() >= 2);
                }

                return 3;
            }

            if (std::holds_alternative<std::vector<rtk::vecxf>>(m.points)) {
                std::vector<rtk::vecxf> const & points = std::get<std::vector<rtk::vecxf>>(m.points);

                if constexpr (scion::ENABLE_ASSERTIONS) {
                    scion::assert_debug(points.size() >= 2);
                }

                scion::isize const dimension = points[0].size();

                if constexpr (scion::ENABLE_ASSERTIONS) {
                    scion::assert_debug(std::ranges::none_of(points, [dimension](rtk::vecxf const & point) -> bool {
                        return point.size() != dimension;
                    }));
                }

                return dimension;
            }

            std::unreachable();
        }

        /// \brief returns the width of the projection
        /// \return width of the projection
        /// \throws std::logic_error on DEBUG or PROFILE if the projection has less than 2 points
        /// \note this value is not cached
        [[nodiscard]] auto width() const -> scion::usize {
            if (std::holds_alternative<std::vector<rtk::vec2f>>(m.points)) {
                if constexpr (scion::ENABLE_ASSERTIONS) {
                    scion::assert_debug(std::get<std::vector<rtk::vec2f>>(m.points).size() >= 2);
                }

                return std::get<std::vector<rtk::vec2f>>(m.points).size();
            }

            if (std::holds_alternative<std::vector<rtk::vec3f>>(m.points)) {
                if constexpr (scion::ENABLE_ASSERTIONS) {
                    scion::assert_debug(std::get<std::vector<rtk::vec3f>>(m.points).size() >= 2);
                }

                return std::get<std::vector<rtk::vec3f>>(m.points).size();
            }

            if (std::holds_alternative<std::vector<rtk::vecxf>>(m.points)) {
                if constexpr (scion::ENABLE_ASSERTIONS) {
                    scion::assert_debug(std::get<std::vector<rtk::vecxf>>(m.points).size() >= 2);
                }

                return std::get<std::vector<rtk::vecxf>>(m.points).size();
            }

            std::unreachable();
        }

        /// \brief returns the width of the projection of the specified dimension
        /// \tparam Dim dimension of the points (2 or 3 for fixed, 0 for dynamic)
        /// \return width of the projection
        /// \throws std::bad_variant_access if the projection is not of the queried dimension
        /// \throws std::logic_error on DEBUG or PROFILE if the projection has less than 2 points
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        [[nodiscard]] auto width() const -> scion::usize {
            if constexpr (scion::ENABLE_ASSERTIONS) {
                scion::assert_debug(std::get<std::vector<rtk::vecf<Dim>>>(m.points).size() >= 2);
            }

            return std::get<std::vector<rtk::vecf<Dim>>>(m.points).size();
        }

        /// \brief projects the given points
        /// \tparam Dim dimension of the points (2 or 3 for fixed, 0 for dynamic)
        /// \param points points to project
        /// \return projected points
        /// \errors the points must be projectable by the projection
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        [[nodiscard]] auto project(std::span<rtk::vecxf> const & points) const -> std::expected<std::vector<rtk::vecf<Dim>>, std::domain_error> {
            std::vector<rtk::vecf<Dim>> out;

            scion::isize const width = static_cast<scion::isize>(this->width());
            scion::isize const dimension = static_cast<scion::isize>(this->dimension());

            if constexpr (Dim == 0) {
                if (std::holds_alternative<std::vector<rtk::vec2f>>(m.points)) {
                    auto const & projection = std::get<std::vector<rtk::vec2f>>(m.points);

                    for (rtk::vecxf const & point : points) {
                        if (point.size() != width) {
                            return std::unexpected<std::domain_error>{"points is not projectable by current projection"};
                        }

                        rtk::vecxf p = rtk::vecxf::Zero(width);
                        for (scion::isize i = 0; i < width; ++i) {
                            p += point(i) * projection[i];
                        }
                        out.emplace_back(p);
                    }

                    return out;
                }

                if (std::holds_alternative<std::vector<rtk::vec3f>>(m.points)) {
                    auto const & projection = std::get<std::vector<rtk::vec3f>>(m.points);

                    for (rtk::vecxf const & point : points) {
                        if (point.size() != width) {
                            return std::unexpected<std::domain_error>{"points is not projectable by current projection"};
                        }

                        rtk::vecxf p = rtk::vecxf::Zero(width);
                        for (scion::isize i = 0; i < width; ++i) {
                            p += point(i) * projection[i];
                        }
                        out.emplace_back(p);
                    }

                    return out;
                }
            }

            if (!std::holds_alternative<std::vector<rtk::vecf<Dim>>>(m.points)) {
                return std::unexpected<std::domain_error>{"points are not projectable by current projection"};
            }

            auto const & projection = std::get<std::vector<rtk::vecf<Dim>>>(m.points);

            for (rtk::vecxf const & point : points) {
                if (point.size() != width) {
                    return std::unexpected<std::domain_error>{"points is not projectable by current projection"};
                }

                rtk::vecf<Dim> p = rtk::vecf<Dim>::Zero(dimension);
                for (scion::isize i = 0; i < width; ++i) {
                    p += point(i) * projection[i];
                }
                out.emplace_back(p);
            }

            return out;
        }

        /// \brief projects the given point
        /// \tparam Dim dimension of the points (2 or 3 for fixed, 0 for dynamic)
        /// \param point point to project
        /// \return projected point
        /// \errors the points must be projectable by the projection
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        [[nodiscard]] auto project(rtk::vecxf const & point) const -> std::expected<rtk::vecf<Dim>, std::domain_error> {
            scion::isize const width = static_cast<scion::isize>(this->width());
            scion::isize const dimension = static_cast<scion::isize>(this->dimension());

            if (point.size() != width) {
                return std::unexpected<std::domain_error>{"points is not projectable by current projection"};
            }

            if constexpr (Dim == 0) {
                if (std::holds_alternative<std::vector<rtk::vec2f>>(m.points)) {
                    std::span<rtk::vec2f const> const projection = std::get<std::vector<rtk::vec2f>>(m.points);

                    rtk::vecxf p = rtk::vecxf::Zero(width);
                    for (scion::isize i = 0; i < width; ++i) {
                        p += point(i) * projection[i];
                    }

                    return p;
                }

                if (std::holds_alternative<std::vector<rtk::vec3f>>(m.points)) {
                    std::span<rtk::vec3f const> const projection = std::get<std::vector<rtk::vec3f>>(m.points);

                    rtk::vecxf p = rtk::vecxf::Zero(width);
                    for (scion::isize i = 0; i < width; ++i) {
                        p += point(i) * projection[i];
                    }

                    return p;
                }
            }

            if (!std::holds_alternative<std::vector<rtk::vecf<Dim>>>(m.points)) {
                return std::unexpected<std::domain_error>{"points is not projectable by current projection"};
            }

            std::span<rtk::vecf<Dim> const> const projection = std::span{std::get<std::vector<rtk::vecf<Dim>>>(m.points)};

            rtk::vecf<Dim> p = rtk::vecf<Dim>::Zero(dimension);
            for (scion::isize i = 0; i < width; ++i) {
                p += point(i) * projection[i];
            }

            return p;
        }

        /// \brief creates a new point
        /// \param point point to create
        /// \return nothing on success, error description on failure
        /// \errors new point is not compatible with the projection
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        [[nodiscard]] auto create_point(rtk::vecf<Dim> && point) -> std::expected<std::monostate, std::invalid_argument> {
            scion::usize const dimension = this->dimension();
            scion::usize const width = this->width();

            if (!std::holds_alternative<std::vector<rtk::vecf<Dim>>>(m.points)) {
                return std::unexpected<std::invalid_argument>{"new point is not compatible with the projection"};
            }

            std::get<std::vector<rtk::vecf<Dim>>>(m.points).emplace_back(std::move(point));

            if (m.name_canonical) {
                m.name = std::format("r{}w{}", dimension, width + 1);
            }

            this->change();

            return std::monostate{};
        }

        /// \brief creates a new point at the specified index
        /// \param point point to create
        /// \param index index at which the point will be inserted
        /// \return nothing on success, error description on failure
        /// \errors the points must be compatible with the projection
        /// \errors the index must be a valid index (from 0 to width included)
        template <scion::usize Dim = 0>
            requires (Dim == 0 || Dim == 2 || Dim == 3)
        [[nodiscard]] auto create_point(vecf<Dim> const point, scion::usize const index) -> std::expected<std::monostate, std::variant<std::out_of_range, std::invalid_argument>> {
            scion::usize const dimension = this->dimension();
            scion::usize const width = this->width();

            if (index > width) {
                return std::unexpected<std::out_of_range>{"index is out of range"};
            }

            if (!std::holds_alternative<std::vector<rtk::vecf<Dim>>>(m.points)) {
                return std::unexpected<std::invalid_argument>{"new point is not compatible with the projection"};
            }

            std::get<std::vector<rtk::vecf<Dim>>>(m.points).insert(std::get<std::vector<rtk::vecf<Dim>>>(m.points).begin() + index, point);

            if (m.name_canonical) {
                m.name = std::format("r{}w{}", dimension, width + 1);
            }

            this->change();

            return std::monostate{};
        }

        /// \brief Deletes the point at the specified index.
        /// \param index Index of the point to delete.
        /// \return Nothing on success, error description on failure.
        /// \errors The index must be in range.
        /// \errors The projection must have a width of at least 2.
        [[nodiscard]] auto delete_point(scion::usize const index) -> std::expected<std::monostate, std::variant<std::out_of_range, std::invalid_argument>> {
            scion::usize const dimension = this->dimension();
            scion::usize const width = this->width();

            if (width <= 2) {
                return std::unexpected<std::out_of_range>("projection must have a width greater than 1");
            }

            if (width <= index) {
                return std::unexpected<std::out_of_range>("index is out of range");
            }

            if (std::holds_alternative<std::vector<vecf<2>>>(m.points)) {
                std::get<std::vector<vecf<2>>>(m.points).erase(std::get<std::vector<vecf<2>>>(m.points).begin() + index);
            }

            if (std::holds_alternative<std::vector<vecf<3>>>(m.points)) {
                std::get<std::vector<vecf<3>>>(m.points).erase(std::get<std::vector<vecf<3>>>(m.points).begin() + index);
            }

            if (std::holds_alternative<std::vector<vecf<0>>>(m.points)) {
                std::get<std::vector<vecf<0>>>(m.points).erase(std::get<std::vector<vecf<0>>>(m.points).begin() + index);
            }

            if (m.name_canonical) {
                m.name = std::format("r{}w{}", dimension, width - 1);
            }

            this->change();

            return std::monostate{};
        }

        /// \brief Returns the name of the projection.
        /// \return Name of the projection. Lifetime valid as long as the name doesn't change.
        [[nodiscard]] auto name() const noexcept -> std::string const & {
            return m.name;
        }

        /// \brief Increments the update counter.
        void change() noexcept {
            m.counter += 1;
        }

        /// \brief Returns the current value of the update counter.
        /// \return Current value of the update counter.
        [[nodiscard]] auto counter() const noexcept -> scion::usize {
            return m.counter;
        }
    };
}
