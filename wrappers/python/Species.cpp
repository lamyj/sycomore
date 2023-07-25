#include <algorithm>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

#include "type_casters.h"

namespace
{

void set_D(sycomore::Species & species, pybind11::object const & value)
{
    if(pybind11::isinstance<sycomore::Quantity>(value))
    {
        species.set_D(value.cast<sycomore::Quantity>());
    }
    else
    {
        species.set_D(value.cast<sycomore::Matrix3x3Q>());
    }
}

sycomore::Species constructor(
    sycomore::Quantity const & R1, sycomore::Quantity const & R2,
    pybind11::object D, sycomore::Quantity const & delta_omega)
{
    sycomore::Species species(R1, R2, {0, sycomore::Diffusion}, delta_omega);
    set_D(species, D);
    return species;
}

}

void wrap_Species(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::units;

    class_<Species>(
            m, "Species", "Species described by its NMR parameters")
        .def(
            init(&constructor),
            arg("R1"), arg("R2"),
            arg("D")=0*units::m*units::m/s, arg("delta_omega")=0*rad/s)
        .def_property(
            "R1", &Species::R1, &Species::set_R1,
            "Longitudinal relaxation rate.")
        .def_property(
            "T1", &Species::T1, &Species::set_T1,
            "Longitudinal relaxation time.")
        .def_property(
            "R2", &Species::R2, &Species::set_R2, 
            "Transversal relaxation rate.")
        .def_property(
            "T2", &Species::T2, &Species::set_T2,
            "Transversal relaxation time.")
        .def_property("D", &Species::D, set_D, "Diffusion tensor.")
        .def_property(
            "delta_omega", &Species::delta_omega, &Species::set_delta_omega,
            "Frequency offset.")
        .def(pickle(
            [](Species const & s) {
                auto const & D = s.D();
                return make_tuple(
                    s.R1(), s.R2(),
                    D.unchecked(0,0), D.unchecked(0,1), D.unchecked(0, 2),
                    D.unchecked(1,0), D.unchecked(1,1), D.unchecked(1,2),
                    D.unchecked(2,0), D.unchecked(2,1), D.unchecked(2,2),
                    s.delta_omega());
            },
            [](tuple t) {
                if(t.size() != 12)
                {
                    throw std::runtime_error("Invalid state!");
                }
                return Species(
                    t[0].cast<Quantity>(), t[1].cast<Quantity>(),
                    {
                        {t[2].cast<Quantity>(), t[3].cast<Quantity>(), t[4].cast<Quantity>()},
                        {t[5].cast<Quantity>(), t[6].cast<Quantity>(), t[7].cast<Quantity>()},
                        {t[8].cast<Quantity>(), t[9].cast<Quantity>(), t[10].cast<Quantity>()},
                    },
                    t[11].cast<Quantity>());
            }
        ));
}
