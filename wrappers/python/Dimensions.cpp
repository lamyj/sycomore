#include <sstream>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>

#include "sycomore/Dimensions.h"

void wrap_Dimensions(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;

    class_<Dimensions>(m, "Dimensions")
        .def(init<double,double,double,double,double,double,double>())
        .def_readwrite("length", &Dimensions::length)
        .def_readwrite("mass", &Dimensions::mass)
        .def_readwrite("time", &Dimensions::time)
        .def_readwrite("electric_current", &Dimensions::electric_current)
        .def_readwrite(
            "thermodynamic_temperature", &Dimensions::thermodynamic_temperature)
        .def_readwrite("amount_of_substance", &Dimensions::amount_of_substance)
        .def_readwrite("luminous_intensity", &Dimensions::luminous_intensity)
        .def(self == self)
        .def(self != self)
        .def(self *= self)
        .def(self /= self)
        .def(self * self)
        .def(self / self)
        .def(
            "__pow__",
            [](Dimensions const & d, int i) { return std::pow(d, i); },
            is_operator())
        .def(
            "__repr__",
            [](Dimensions const & d) {
                std::ostringstream s;
                s << d;
                return s.str();
            });

    m.attr("Length") = Length;
    m.attr("Time") = Time;
    m.attr("Mass") = Mass;
    m.attr("ElectricCurrent") = ElectricCurrent;
    m.attr("ThermodynamicTemperature") = ThermodynamicTemperature;
    m.attr("AmountOfSubstance") = AmountOfSubstance;
    m.attr("LuminousIntensity") = LuminousIntensity;

    m.attr("Surface") = Surface;
    m.attr("Volume") = Volume;

    m.attr("Velocity") = Velocity;
    m.attr("Acceleration") = Acceleration;

    m.attr("Angle") = Angle;
    m.attr("SolidAngle") = SolidAngle;

    m.attr("Frequency") = Frequency;
    m.attr("Force") = Force;
    m.attr("Pressure") = Pressure;
    m.attr("Energy") = Energy;
    m.attr("Power") = Power;
    m.attr("ElectricCharge") = ElectricCharge;
    m.attr("Voltage") = Voltage;
    m.attr("Capacitance") = Capacitance;
    m.attr("Resistance") = Resistance;
    m.attr("ElectricalConductance") = ElectricalConductance;
    m.attr("MagneticFlux") = MagneticFlux;
    m.attr("MagneticFluxDensity") = MagneticFluxDensity;
    m.attr("Inductance") = Inductance;
    m.attr("LuminousFlux") = LuminousFlux;
    m.attr("Illuminance") = Illuminance;
    m.attr("Radioactivity") = Radioactivity;
    m.attr("AbsorbedDose") = AbsorbedDose;
    m.attr("EquivalentDose") = EquivalentDose;
    m.attr("CatalyticActivity") = CatalyticActivity;

    m.attr("AngularFrequency") = AngularFrequency;
}
