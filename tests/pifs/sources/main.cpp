#include "scion/scion.hpp"
#include "scion/prelude/prelude.hpp"

#include "rtk/rtk.hpp"
#include "rtk/model/model.hpp"
#include "rtk/prelude/prelude.hpp"

#include <iostream>

auto main() -> int {
    rtk::pifs model = rtk::pifs<>::make(2, 3).value();

    {
        rtk::log::info("test", "populating model");
        model.eval_mut(0, 0) = 1.00f;
        model.eval_mut(0, 1) = 0.50f;
        model.eval_mut(0, 2) = 0.25f;
        model.evec_mut(0, 0) = (vecxf(3) <<  1.0f,  0.0f,  0.0f).finished();
        model.evec_mut(0, 1) = (vecxf(3) << -1.0f,  1.0f,  0.0f).finished();
        model.evec_mut(0, 2) = (vecxf(3) <<  1.0f, -2.0f,  1.0f).finished();

        model.eval_mut(1, 0) = 1.00f;
        model.eval_mut(1, 1) = 0.50f;
        model.eval_mut(1, 2) = 0.25f;
        model.evec_mut(1, 0) = (vecxf(3) <<  0.0f,  0.0f,  1.0f).finished();
        model.evec_mut(1, 1) = (vecxf(3) <<  0.0f, -1.0f,  1.0f).finished();
        model.evec_mut(1, 2) = (vecxf(3) <<  1.0f, -2.0f,  1.0f).finished();

        model.request_update_transformations();

        static_cast<void>(model.validate<rtk::log::error_t>());
        model.update();

        std::cout << model.tran(0) << std::endl;
        std::cout << model.tran(1) << std::endl;
        std::cout << model.evec(0, 0) << std::endl;
        std::cout << model.evec(0, 1) << std::endl;
        std::cout << model.evec(0, 2) << std::endl;
        std::cout << model.evec(1, 0) << std::endl;
        std::cout << model.evec(1, 1) << std::endl;
        std::cout << model.evec(1, 2) << std::endl;
        std::cout << model.eval(0, 1) << std::endl;
        std::cout << model.eval(0, 2) << std::endl;
        std::cout << model.eval(1, 0) << std::endl;
        std::cout << model.eval(1, 1) << std::endl;
        std::cout << model.eval(1, 2) << std::endl;
    }

    {
        rtk::log::info("test", "populating model eigen");

        std::vector<std::vector<f32>> eigenvals;
        eigenvals.emplace_back(model.evals(0).begin(), model.evals(0).end());
        eigenvals.emplace_back(model.evals(1).begin(), model.evals(1).end());
        std::vector<std::vector<vecxf>> eigenvecs;
        eigenvecs.emplace_back(model.evecs(0).begin(), model.evecs(0).end());
        eigenvecs.emplace_back(model.evecs(1).begin(), model.evecs(1).end());

        rtk::pifs model_eigen = rtk::pifs<>::make(std::move(eigenvals), std::move(eigenvecs)).value();

        static_cast<void>(model_eigen.validate<rtk::log::error_t>());

        std::cout << model_eigen.tran(0) << std::endl;
        std::cout << model_eigen.tran(1) << std::endl;
        std::cout << model_eigen.evec(0, 0) << std::endl;
        std::cout << model_eigen.evec(0, 1) << std::endl;
        std::cout << model_eigen.evec(0, 2) << std::endl;
        std::cout << model_eigen.evec(1, 0) << std::endl;
        std::cout << model_eigen.evec(1, 1) << std::endl;
        std::cout << model_eigen.evec(1, 2) << std::endl;
        std::cout << model_eigen.eval(0, 1) << std::endl;
        std::cout << model_eigen.eval(0, 2) << std::endl;
        std::cout << model_eigen.eval(1, 0) << std::endl;
        std::cout << model_eigen.eval(1, 1) << std::endl;
        std::cout << model_eigen.eval(1, 2) << std::endl;
    }

    {
        rtk::log::info("test", "populating model trans");

        std::vector<matxf> trans;
        trans.emplace_back(model.tran(0));
        trans.emplace_back(model.tran(1));

        rtk::pifs model_trans = rtk::pifs<>::make(std::move(trans)).value();
        static_cast<void>(model_trans.validate<rtk::log::error_t>());

        std::cout << model_trans.tran(0) << std::endl;
        std::cout << model_trans.tran(1) << std::endl;
        std::cout << model_trans.evec(0, 0) << std::endl;
        std::cout << model_trans.evec(0, 1) << std::endl;
        std::cout << model_trans.evec(0, 2) << std::endl;
        std::cout << model_trans.evec(1, 0) << std::endl;
        std::cout << model_trans.evec(1, 1) << std::endl;
        std::cout << model_trans.evec(1, 2) << std::endl;
        std::cout << model_trans.eval(0, 0) << std::endl;
        std::cout << model_trans.eval(0, 1) << std::endl;
        std::cout << model_trans.eval(0, 2) << std::endl;
        std::cout << model_trans.eval(1, 0) << std::endl;
        std::cout << model_trans.eval(1, 1) << std::endl;
        std::cout << model_trans.eval(1, 2) << std::endl;
    }

    return 0;
}
