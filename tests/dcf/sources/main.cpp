#include "scion/prelude/prelude.hpp"

#include "rtk/dif/dcf.hpp"

#include <Eigen/Core>

#include <span>
#include <vector>
#include <iostream>

struct A {

    f32 x;
    f32 y;

    explicit A() {}
    explicit A(f32 x, f32 y) {}

    A(A const& a) {}
    A(A && A) {}

    A& operator=(A const& a) { return *this;}
    A& operator=(A && A) { return *this;}

public:

    friend A operator+(A const lhs, const A& rhs) {
        return A {
                lhs.x + rhs.x,
                lhs.y + rhs.y,
        };
    }
};

namespace rtk {
    template<typename S>
    auto vector_product_tensor_flatten(A const& lhs, A rhs) -> A {
        return lhs + rhs;
    }
}

auto main() -> int {
    constexpr usize width_proj = 3;

    std::vector<f32> eigenvalues;
    eigenvalues.emplace_back(1.0);
    eigenvalues.emplace_back(0.5);
    eigenvalues.emplace_back(0.25);

    std::vector<Eigen::VectorX<f32>> eigenvectors;
    eigenvectors.emplace_back(width_proj);
    eigenvectors.back() <<  1.0,  0.0,  0.0;
    eigenvectors.emplace_back(width_proj);
    eigenvectors.back() << -1.0,  1.0,  0.0;
    eigenvectors.emplace_back(width_proj);
    eigenvectors.back() <<  1.0, -2.0,  1.0;

    Eigen::MatrixX<scion::f32> L = Eigen::MatrixX<scion::f32>::Zero(static_cast<i32>(width_proj), static_cast<i32>(width_proj));
    Eigen::MatrixX<scion::f32> V = Eigen::MatrixX<scion::f32>::Zero(static_cast<i32>(width_proj), static_cast<i32>(width_proj));

    for (isize d = 0; d < width_proj; ++d) {
        L(d, d) = eigenvalues[d];

        for (isize e = 0; e < width_proj; ++e) {
            V(d, e) = eigenvectors[e][d];
        }
    }

    Eigen::MatrixX<f32> transformation = V * L * V.inverse();

    std::cout << transformation << std::endl;

    return 0;
}
