#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "sycomore/sycomore.h"

void wrap_Dimensions(pybind11::module &);
void wrap_Quantity(pybind11::module &);
void wrap_units(pybind11::module &);

void wrap_Array(pybind11::module &);

void wrap_magnetization(pybind11::module &);

void wrap_GridScanner(pybind11::module &);
void wrap_Grid(pybind11::module &);

void wrap_Pulse(pybind11::module &);
void wrap_HardPulseApproximation(pybind11::module &);
void wrap_Species(pybind11::module &);
void wrap_TimeInterval(pybind11::module &);

void wrap_como(pybind11::module &);
void wrap_epg(pybind11::module &);

PYBIND11_MODULE(_sycomore, _sycomore)
{
    wrap_Dimensions(_sycomore);
    wrap_Quantity(_sycomore);
    wrap_units(_sycomore);

    wrap_Array(_sycomore);

    wrap_magnetization(_sycomore);

    wrap_GridScanner(_sycomore);
    wrap_Grid(_sycomore);

    wrap_Pulse(_sycomore);
    wrap_HardPulseApproximation(_sycomore);
    wrap_Species(_sycomore);
    wrap_TimeInterval(_sycomore);

    wrap_como(_sycomore);
    wrap_epg(_sycomore);

    using namespace pybind11;
    using namespace sycomore;

    _sycomore.attr("gamma") = sycomore::gamma;

    _sycomore.def(
        "linspace",
        static_cast<std::vector<Quantity> (*) (Quantity, Quantity, std::size_t)>(
            linspace<Quantity>));
    _sycomore.def(
        "linspace",
        static_cast<std::vector<Quantity> (*) (Quantity, std::size_t)>(
            linspace<Quantity>));

    _sycomore.def(
        "linspace",
        static_cast<std::vector<Point> (*) (Point, Point, std::size_t)>(
            linspace<Point>));
    _sycomore.def(
        "linspace",
        static_cast<std::vector<Point> (*) (Point, std::size_t)>(
            linspace<Point>));
}
