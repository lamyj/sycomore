#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#define FORCE_IMPORT_ARRAY
#include <xtensor-python/pycontainer.hpp>

#include "sycomore/sycomore.h"

#include "type_casters.h"

void wrap_Dimensions(pybind11::module &);
void wrap_Quantity(pybind11::module &);
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
        overload_cast<Quantity, Quantity, std::size_t>(linspace<Quantity>));
    _sycomore.def(
        "linspace",
        overload_cast<Quantity, std::size_t>(linspace<Quantity>));

    _sycomore.def(
        "linspace",
        overload_cast<ArrayQ, ArrayQ, std::size_t>(sycomore::linspace<ArrayQ>));
    _sycomore.def(
        "linspace", [](ArrayQ span, std::size_t size){
            return sycomore::linspace(
                xt::eval(-span/2.), xt::eval(+span/2.), size);
        });
}
