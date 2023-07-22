#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "sycomore/epg/operators.h"

#include "../type_casters.h"

template<typename Container, typename Shape>
xt::xtensor_fixed<typename Container::value_type, Shape>
as_xtensor_fixed(Container const & array, Shape &&)
{
    xt::xtensor_fixed<typename Container::value_type, Shape> result;
    std::copy(array.begin(), array.begin()+result.size(), result.begin());
    return result;
}

void wrap_epg_operators(pybind11::module & m)
{
    using namespace pybind11;
    using namespace pybind11::literals;
    using namespace sycomore;
    using namespace sycomore::epg::operators;

    auto operators = m.def_submodule("operators");
    operators.def(
        "pulse_single_pool",
        [](Real angle, Real phase) {
            return as_xtensor_fixed(
                pulse_single_pool(angle, phase), xt::xshape<3, 3>{});
        },
        "angle"_a, "phase"_a);
    
    operators.def(
        "pulse_exchange",
        [](Real angle_a, Real phase_a, Real angle_b, Real phase_b) {
            return as_xtensor_fixed(
                pulse_exchange(angle_a, phase_a, angle_b, phase_b),
                xt::xshape<2, 3, 3>());
        },
        "angle_a"_a, "phase_a"_a, "angle_b"_a, "phase_b"_a);
    
    operators.def(
        "pulse_magnetization_transfer",
        [](Real angle_a, Real phase_a, Real saturation) {
            auto const value = pulse_magnetization_transfer(
                angle_a, phase_a, saturation);
            return std::make_pair(
                as_xtensor_fixed(value, xt::xshape<3, 3>()), value.back());
        },
        "angle_a"_a, "phase_a"_a, "saturation"_a);
    
    operators.def(
        "relaxation_single_pool", &relaxation_single_pool,
        "R1"_a, "R2"_a, "duration"_a);
    
    operators.def(
        "relaxation_exchange",
        [](
                Real R1_a, Real R2_a, Real R1_b, Real R2_b,
                Real k_a, Real k_b, Real delta_b, Real M0_a, Real M0_b,
                Real duration) {
            auto const result = relaxation_exchange(
                R1_a, R2_a, R1_b, R2_b, k_a, k_b, delta_b, M0_a, M0_b, duration);
            
            xt::xtensor_fixed<Complex, xt::xshape<4, 4>> Xi_T;
            Xi_T.fill(0.);
            auto const & Xi_T_linear = std::get<0>(result);
            Xi_T.unchecked(0, 0) = Xi_T_linear[0];
            Xi_T.unchecked(0, 2) = Xi_T_linear[1];
            Xi_T.unchecked(1, 1) = Xi_T_linear[2];
            Xi_T.unchecked(1, 3) = Xi_T_linear[3];
            Xi_T.unchecked(2, 0) = Xi_T_linear[4];
            Xi_T.unchecked(2, 2) = Xi_T_linear[5];
            Xi_T.unchecked(3, 1) = Xi_T_linear[6];
            Xi_T.unchecked(3, 3) = Xi_T_linear[7];
            
            return std::make_tuple(
                Xi_T, as_xtensor_fixed(std::get<1>(result), xt::xshape<2, 2>()),
                std::get<2>(result));
        },
        "R1_a"_a, "R2_a"_a, "R1_b"_a, "R2_b"_a, "k_a"_a, "k_b"_a, "delta_b"_a,
        "M0_a"_a, "M0_b"_a, "duration"_a);
    
    operators.def(
        "relaxation_magnetization_transfer",
        [](
                Real R1_a, Real R2_a, Real R1_b, Real k_a, Real k_b, 
                Real M0_a, Real M0_b, Real duration) {
            auto const result = relaxation_magnetization_transfer(
                R1_a, R2_a, R1_b, k_a, k_b, M0_a, M0_b, duration);
            
            return std::make_tuple(
                std::get<0>(result),
                as_xtensor_fixed(std::get<1>(result), xt::xshape<2, 2>()),
                std::get<2>(result));
        },
        "R1_a"_a, "R2_a"_a, "R1_b"_a, "k_a"_a, "k_b"_a, "M0_a"_a, "M0_b"_a,
        "duration"_a);
}
