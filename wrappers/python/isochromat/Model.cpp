#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <xtensor-python/pytensor.hpp>

#include "sycomore/isochromat/Model.h"

#include "../Quantity.h"

namespace
{

sycomore::isochromat::Model
constructor(
    pybind11::object T1, pybind11::object T2, pybind11::object M0,
    pybind11::object positions, pybind11::object delta_omega)
{
    using sycomore::Quantity;
    using sycomore::TensorR;
    using sycomore::TensorQ;
    using sycomore::wrappers::as_quantity;
    
    try
    {
        auto const T1_ = T1.cast<Quantity const &>();
        if(delta_omega.is(pybind11::none()))
        {
            return {
                T1_, T2.cast<Quantity const &>(),
                M0.cast<TensorR<1>>(), as_quantity<TensorQ<2>>(positions)};
        }
        else
        {
            return {
                T1_, T2.cast<Quantity const &>(),
                M0.cast<TensorR<1>>(), as_quantity<TensorQ<2>>(positions),
                delta_omega.cast<Quantity const &>()};
        }
    }
    catch(pybind11::cast_error &)
    {
        if(delta_omega.is(pybind11::none()))
        {
            return {
                as_quantity<TensorQ<1>>(T1), as_quantity<TensorQ<1>>(T2),
                M0.cast<TensorR<2>>(), as_quantity<TensorQ<2>>(positions)};
        }
        else
        {
            return {
                as_quantity<TensorQ<1>>(T1), as_quantity<TensorQ<1>>(T2),
                M0.cast<TensorR<2>>(), as_quantity<TensorQ<2>>(positions),
                as_quantity<TensorQ<1>>(delta_omega)
            };
        }
    }
}

sycomore::isochromat::Operator
build_pulse(
    sycomore::isochromat::Model const & model,
    pybind11::object angle, pybind11::object phase)
{
    using sycomore::Quantity;
    using sycomore::TensorQ;
    using sycomore::wrappers::as_quantity;
    
    try
    {
        if(phase.is(pybind11::none()))
        {
            return model.build_pulse(
                angle.cast<Quantity const &>());
        }
        else
        {
            return model.build_pulse(
                angle.cast<Quantity const &>(),
                phase.cast<Quantity const &>());
        }
    }
    catch(pybind11::cast_error &)
    {
        if(phase.is(pybind11::none()))
        {
            return model.build_pulse(as_quantity<TensorQ<1>>(angle));
        }
        else
        {
            return model.build_pulse(
                as_quantity<TensorQ<1>>(angle), as_quantity<TensorQ<1>>(phase));
        }
    }
}

sycomore::isochromat::Operator
build_phase_accumulation(
    sycomore::isochromat::Model const & model, pybind11::object angle)
{
    using sycomore::Quantity;
    using sycomore::TensorQ;
    using sycomore::wrappers::as_quantity;
    
    try
    {
        return model.build_phase_accumulation(angle.cast<Quantity const &>());
    }
    catch(pybind11::cast_error &)
    {
        return model.build_phase_accumulation(as_quantity<TensorQ<1>>(angle));
    }
}

sycomore::isochromat::Operator
build_time_interval(
    sycomore::isochromat::Model const & model,
    sycomore::Quantity const & duration, pybind11::object delta_omega,
    pybind11::object gradient)
{
    using sycomore::Quantity;
    using sycomore::TensorQ;
    using sycomore::wrappers::as_quantity;
    
    if(delta_omega.is(pybind11::none()) && gradient.is(pybind11::none()))
    {
        return model.build_time_interval(duration);
    }
    else if(gradient.is(pybind11::none()))
    {
        try
        {
            return model.build_time_interval(
                duration, delta_omega.cast<Quantity>());
        }
        catch(pybind11::cast_error &)
        {
            return model.build_time_interval(
                duration, as_quantity<TensorQ<1>>(delta_omega));
        }
    }
    else if(delta_omega.is(pybind11::none()))
    {
        using namespace sycomore::units;
        
        try
        {
            auto const gradient_ = as_quantity<TensorQ<1>>(gradient);
            Quantity const delta_omega_ = 0*rad/s;
            return model.build_time_interval(duration, delta_omega_, gradient_);
        }
        catch(pybind11::cast_error &)
        {
            auto const gradient_ = as_quantity<TensorQ<2>>(gradient);
            TensorQ<1> delta_omega_(TensorQ<1>::shape_type{gradient_.shape()[0]});
            delta_omega_.fill(0*rad/s);
            return model.build_time_interval(duration, delta_omega_, gradient_);
        }
    }
    else
    {
        try
        {
            auto const delta_omega_ = delta_omega.cast<Quantity>();
            return model.build_time_interval(
                duration,
                delta_omega_,
                as_quantity<TensorQ<1>>(gradient));
        }
        catch(pybind11::cast_error &)
        {
            return model.build_time_interval(
                duration,
                as_quantity<TensorQ<1>>(delta_omega),
                as_quantity<TensorQ<2>>(gradient));
        }
    }
}

}

void wrap_isochromat_Model(pybind11::module & m)
{
    using namespace pybind11;
    using namespace pybind11::literals;
    using namespace sycomore;
    using namespace sycomore::isochromat;

    class_<Model>(m, "Model")
        .def(
            init(&constructor),
            "T1"_a, "T2"_a, "M0"_a, "positions"_a, "delta_omega"_a=none(),
            "Create a model")
        .def(
            "build_pulse", &build_pulse,
            "angle"_a, "phase"_a=none(),
            "Create an RF pulse operator")
        .def(
            "build_time_interval", &build_time_interval,
            "duration"_a, "delta_omega"_a=none(), "gradient"_a=none(),
            "Create a time interval operator")
        .def(
            "build_time_interval",
            overload_cast<TimeInterval const &>(
                &Model::build_time_interval, const_),
            "time_interval"_a,
            "Create a time interval operator")
        .def(
            "build_relaxation", &Model::build_relaxation, "duration"_a,
            "Create a relaxation operator")
        .def(
            "build_phase_accumulation", &build_phase_accumulation,
            "angle"_a,
            "Create a phase accumulation operator")
        .def(
            "apply", &Model::apply, "operator"_a,
            "Apply an operator to the magnetization")
        .def_property_readonly("T1", &Model::T1, "T1 field")
        .def_property_readonly("R1", &Model::R1, "R1 field")
        .def_property_readonly("T2", &Model::T2, "T2 field")
        .def_property_readonly("R2", &Model::R2, "R2 field")
        .def_property_readonly("M0", &Model::M0, "M0 field")
        .def_property_readonly(
            "delta_omega", &Model::delta_omega, "Off-resonance field")
        .def_property_readonly(
            "magnetization", &Model::magnetization, "Magnetization field")
        .def_property_readonly(
            "positions", &Model::positions, "Positions of the isochromats");
}
