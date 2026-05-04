// #include "rtk/model/pifs.hpp"
//
// #include "Eigen/Eigen"
//
// #include <ranges>
//
// #include "rtk/log.hpp"
//
// void rtk::pifs::update_transformations() {
//     rtk::log::trace("model", "updating transformations of {}", static_cast<void *>(this));
//
//     for (scion::usize i_tran = 0; i_tran < this->m_width_tran; ++i_tran) {
//         Eigen::MatrixX<scion::f32> L;
//         L.setZero(static_cast<Eigen::Index>(this->m_width_proj), static_cast<Eigen::Index>(this->m_width_proj));
//
//         Eigen::MatrixX<scion::f32> V;
//         V.setZero(static_cast<Eigen::Index>(this->m_width_proj), static_cast<Eigen::Index>(this->m_width_proj));
//
//         for (scion::isize d = 0; d < this->m_width_proj; ++d) {
//             L(d, d) = this->m_eigenvalues[i_tran][d];
//
//             for (scion::isize e = 0; e < this->m_width_proj; ++e) {
//                 V(d, e) = this->m_eigenvectors[i_tran][e][d];
//             }
//         }
//
//         this->m_transformations[i_tran] = V * L * V.inverse();
//     }
//
//     this->m_update_processed.transformations = true;
//     this->request_update_attractor_bn();
// }
//
// void rtk::pifs::update_eigenproperties() {
//     rtk::log::trace("model", "updating eigenproperties of {}", static_cast<void *>(this));
//
//     for (scion::usize i_tran = 0; i_tran < this->m_width_tran; ++i_tran) {
//         this->m_eigenvalues[i_tran].clear();
//         this->m_eigenvectors[i_tran].clear();
//
//         Eigen::EigenSolver<Eigen::MatrixX<scion::f32>> solver{this->m_transformations[i_tran]};
//         Eigen::VectorX<scion::f32> solver_output_eigenvalues = solver.eigenvalues().real();
//         Eigen::MatrixX<scion::f32> solver_output_eigenvectors = solver.eigenvectors().real();
//
//         for(scion::usize i = 0; i < this->m_width_proj; ++i) {
//             this->m_eigenvalues[i_tran].push_back(solver_output_eigenvalues[i]);
//             this->m_eigenvectors[i_tran].push_back(solver_output_eigenvectors.col(i));
//         }
//
//         std::ranges::sort(std::ranges::zip_view(this->m_eigenvalues[i_tran], this->m_eigenvectors[i_tran]), std::less<>{}, [](auto const& a) {
//             return -std::abs(std::get<0>(a));
//         });
//     }
//
//     this->m_update_processed.eigenproperties = true;
//     this->request_update_attractor_bn();
// }
//
// void update_attractor_bn_recurse(std::vector<Eigen::VectorX<scion::f32>> & attractor, std::vector<Eigen::MatrixX<scion::f32>> const& transformations, std::pair<Eigen::VectorX<scion::f32> const&, Eigen::VectorX<scion::f32> const&> primitive, Eigen::MatrixX<scion::f32> const& matrix, scion::usize iteration) {
//     if (iteration == 0) {
//         attractor.emplace_back(matrix * std::get<0>(primitive));
//         attractor.emplace_back(matrix * std::get<1>(primitive));
//     } else {
//         for (Eigen::MatrixX<scion::f32> const& transformation : transformations) {
//             update_attractor_bn_recurse(attractor, transformations, primitive, matrix * transformation, iteration - 1);
//         }
//     }
// }
//
// void rtk::pifs::update_attractor_bn() {
//     rtk::log::trace("model", "updating attractor bn of {}", static_cast<void *>(this));
//
//     if (this->m_transformations.empty()) {
//         throw std::logic_error("transformations are missing");
//     }
//
//     this->m_attractor_bn.clear();
//
//     Eigen::MatrixX<scion::f32> matrix(this->m_width_proj, this->m_width_proj);
//     matrix.setIdentity();
//
//     update_attractor_bn_recurse(this->m_attractor_bn, this->m_transformations, std::make_pair(this->m_eigenvectors.front().front(), this->m_eigenvectors.back().front()), matrix, this->m_data_model.iterations);
//
//     this->m_update_processed.attractor_bn = true;
//     this->request_update_attractor_rn();
// }
//
// void rtk::pifs::update_attractor_rn() {
//     rtk::log::trace("model", "updating transformations of {}", static_cast<void *>(this));
//
//     if (this->m_attractor_bn.empty()) {
//         throw std::logic_error("modeling space attractor is missing");
//     }
//
//     this->m_attractor_rn.emplace<std::monostate>();
//
//     if (std::holds_alternative<std::weak_ptr<rtk::projection<2>>>()) {
//         if (auto result = this->m_projection_d2->get()->project(this->m_attractor_bn); result.has_value()) {
//             this->m_attractor_rn.emplace<std::vector<Eigen::Vector<scion::f32, 2>>>(result.value());
//         } else {
//             dcf::log::error("error", "could not compute projected attractor: {}", result.error().what());
//         }
//     }
//
//     if (std::holds_alternative<std::weak_ptr<rtk::projection<3>>>(this->m_projection)) {
//         if (std::get<std::weak_ptr<rtk::projection<3>>>(this->m_projection).expired()) {
//             this->unset_projection();
//
//             return;
//         }
//
//         const auto ptr = std::get<std::weak_ptr<rtk::projection<3>>>(this->m_projection).lock();
//
//         if (auto result = ptr->project(this->m_attractor_bn); result.has_value()) {
//             this->m_attractor_rn.emplace<std::vector<Eigen::Vector<f32, 3>>>(result.value());
//         } else {
//             dcf::log::error("error", "could not compute projected attractor: {}", result.error().what());
//         }
//
//         return;
//     }
//
//     if (std::holds_alternative<std::weak_ptr<rtk::projection<0>>>(this->m_projection)) {
//         if (std::get<std::weak_ptr<rtk::projection<0>>>(this->m_projection).expired()) {
//             this->unset_projection();
//
//             return;
//         }
//
//         const auto ptr = std::get<std::weak_ptr<rtk::projection<0>>>(this->m_projection).lock();
//
//         throw std::logic_error("computing an arbitrary dimensional projected attractor is not supported");
//
//         return;
//     }
//
//     this->m_update_processed.attractor_rn = true;
// }
