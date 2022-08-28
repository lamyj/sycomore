#include <algorithm>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/Species.h"
#include "sycomore/units.h"

namespace
{

void set_D(sycomore::Species & species, pybind11::object const & value)
{
    if(pybind11::isinstance<pybind11::array>(value))
    {
        auto numpy = pybind11::module::import("numpy");
        auto raveled = numpy.attr("ravel")(
                value.cast<pybind11::array>()
            ).attr("tolist")();
        sycomore::Array<sycomore::Quantity> array(pybind11::len(raveled));
        std::transform(
            raveled.begin(), raveled.end(), array.begin(),
            [](pybind11::handle const & x) {
                return x.cast<sycomore::Quantity>(); });
        species.set_D(array);
    }
    else if(pybind11::isinstance<pybind11::sequence>(value))
    {
        sycomore::Array<sycomore::Quantity> array(pybind11::len(value));
        std::transform(
            value.begin(), value.end(), array.begin(),
            [](pybind11::handle const & x) {
                return x.cast<sycomore::Quantity>(); });
        species.set_D(array);
    }
    else
    {
        species.set_D(value.cast<sycomore::Quantity>());
    }
}

sycomore::Species constructor(
    sycomore::Quantity const & R1, sycomore::Quantity const & R2,
    pybind11::object D, sycomore::Quantity const & R2_prime,
    sycomore::Quantity const & delta_omega, sycomore::Real w)
{
    sycomore::Species species(
        R1, R2, {0, sycomore::Diffusion}, R2_prime, delta_omega, w);
    set_D(species, D);
    return species;
}

}

void wrap_Species(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::units;

    class_<Species>(m, "Species")
        .def(
            init(&constructor),
            arg("R1"), arg("R2"),
            arg("D")=0*units::m*units::m/s, arg("R2_prime")=0_Hz,
            arg("delta_omega")=0*rad/s, arg("w")=1)
        .def_property(
            "R1", &Species::get_R1, &Species::set_R1,
            "Longitudinal relaxation rate.")
        .def_property_readonly(
            "T1", &Species::get_T1, "Longitudinal relaxation time.")
        .def_property(
            "R2", &Species::get_R2, &Species::set_R2, 
            "Transversal relaxation rate.")
        .def_property_readonly(
            "T2", &Species::get_T2, "Transversal relaxation time.")
        .def_property(
            "D", &Species::get_D, set_D, "Diffusion tensor.")
        .def_property(
            "R2_prime", &Species::get_R2_prime, &Species::set_R2_prime,
            "The part of the apparent transversal relaxation rate R2* "
            "attributed to the magnetic field inhomogeneity")
        .def_property_readonly(
            "T2_prime", &Species::get_T2_prime,
            "The part of the apparent transversal relaxation time T2* "
            "attributed to the magnetic field inhomogeneity")
        .def_property(
            "delta_omega", &Species::get_delta_omega, &Species::set_delta_omega,
            "Frequency offset.")
        .def_readwrite("w", &Species::w, "Relative weight.")
        .def(pickle(
            [](Species const & s) {
                auto const D = s.get_D();
                return make_tuple(
                    s.get_R1(), s.get_R2(),
                    D[0], D[1], D[2], D[3], D[4], D[5], D[6], D[7], D[8],
                    s.get_R2_prime(),
                    s.get_delta_omega(), s.w);
            },
            [](tuple t) {
                if(t.size() != 14)
                {
                    throw std::runtime_error("Invalid state!");
                }
                return Species(
                    t[0].cast<Quantity>(), t[1].cast<Quantity>(),
                    {
                        t[2].cast<Quantity>(), t[3].cast<Quantity>(), t[4].cast<Quantity>(),
                        t[5].cast<Quantity>(), t[6].cast<Quantity>(), t[7].cast<Quantity>(),
                        t[8].cast<Quantity>(), t[9].cast<Quantity>(), t[10].cast<Quantity>(),
                    },
                    t[11].cast<Quantity>(),
                    t[12].cast<Quantity>(), t[13].cast<Real>());
            }
        ));
}
