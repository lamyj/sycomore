#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pycontainer.hpp>

#include "sycomore/sycomore.h"

#include "type_casters.h"

void wrap_Dimensions(pybind11::module &);
void wrap_Quantity(pybind11::module &);
void wrap_Vector2Q(pybind11::module &);
void wrap_Vector3Q(pybind11::module &);
void wrap_Vector4Q(pybind11::module &);
void wrap_Matrix2x2Q(pybind11::module &);
void wrap_Matrix3x3Q(pybind11::module &);
void wrap_TensorQ1(pybind11::module &);
void wrap_TensorQ2(pybind11::module &);
void wrap_TensorQ3(pybind11::module &);
void wrap_TensorQ4(pybind11::module &);
void wrap_ArrayQ(pybind11::module &);
void wrap_units(pybind11::module &);

void wrap_Pulse(pybind11::module &);
void wrap_HardPulseApproximation(pybind11::module &);
void wrap_Species(pybind11::module &);
void wrap_TimeInterval(pybind11::module &);

void wrap_epg(pybind11::module &);
void wrap_isochromat(pybind11::module &);

PYBIND11_MODULE(_sycomore, _sycomore)
{
    xt::import_numpy();
    
    wrap_Dimensions(_sycomore);
    wrap_Quantity(_sycomore);
    wrap_Vector2Q(_sycomore);
    wrap_Vector3Q(_sycomore);
    wrap_Vector4Q(_sycomore);
    wrap_Matrix2x2Q(_sycomore);
    wrap_Matrix3x3Q(_sycomore);
    wrap_TensorQ1(_sycomore);
    wrap_TensorQ2(_sycomore);
    wrap_TensorQ3(_sycomore);
    wrap_TensorQ4(_sycomore);
    wrap_ArrayQ(_sycomore);
    wrap_units(_sycomore);

    wrap_Pulse(_sycomore);
    wrap_HardPulseApproximation(_sycomore);
    wrap_Species(_sycomore);
    wrap_TimeInterval(_sycomore);

    wrap_epg(_sycomore);
    wrap_isochromat(_sycomore);

    using namespace pybind11;
    using namespace sycomore;

    _sycomore.attr("gamma") = sycomore::gamma;
    _sycomore.attr("gamma_bar") = sycomore::gamma_bar;

    _sycomore.def(
        "linspace", 
        overload_cast<Quantity, Quantity, std::size_t>(linspace<Quantity>),
        "Generate evenly-spaced samples");
    _sycomore.def(
        "linspace",
        overload_cast<Quantity, std::size_t>(linspace<Quantity>),
        "Generate evenly-spaced samples");

    _sycomore.def(
        "linspace",
        overload_cast<ArrayQ, ArrayQ, std::size_t>(sycomore::linspace<ArrayQ>),
        "Generate evenly-spaced samples");
    _sycomore.def(
        "linspace", [](ArrayQ span, std::size_t size){
            return sycomore::linspace(
                xt::eval(-span/2.), xt::eval(+span/2.), size);
        },
        "Generate evenly-spaced samples");
    
    _sycomore.def(
        "round", sycomore::round<Quantity>, "Round to x nearest multiple of r");
}
