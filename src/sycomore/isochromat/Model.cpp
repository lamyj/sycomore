#include "Model.h"

#include <array>
#include <cmath>
#include <stdexcept>

#include <xtensor/xbuilder.hpp>
#include <xtensor/xmath.hpp>
#include <xtensor/xmanipulation.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

#include "sycomore/Quantity.h"
#include "sycomore/sycomore.h"
#include "sycomore/units.h"
#include "sycomore/isochromat/Operator.h"

namespace sycomore
{

namespace isochromat
{

Model
::Model(
    Quantity const & T1, Quantity const & T2, TensorR<1> const & M0,
    TensorQ<2> const & positions)
: Model(
    xt::repeat(TensorQ<1>{T1}, positions.shape()[0], 0),
    xt::repeat(TensorQ<1>{T2}, positions.shape()[0], 0),
    xt::eval(xt::repeat(xt::atleast_2d(M0), positions.shape()[0], 0)),
    positions)
{
    // Nothing else
}

Model
::Model(
    TensorQ<1> const & T1, TensorQ<1> const & T2, TensorR<2> const & M0,
    TensorQ<2> const & positions)
: _T1(T1.shape()), _T2(T2.shape()), _M0(xt::view(M0, xt::all(), 2)),
    _magnetization(TensorR<2>::shape_type{M0.shape()[0], 4}),
    _positions(positions.shape())
{
    if(
        T1.size() != T2.size() || T1.size() != M0.shape()[0]
        || T1.size() != positions.shape()[0])
    {
        throw std::runtime_error("Size mismatch");
    }
    this->_T1 = convert_to(T1, units::s);
    this->_T2 = convert_to(T2, units::s);
    
    xt::view(this->_magnetization, xt::all(), xt::range(0, 3)) = M0;
    xt::view(this->_magnetization, xt::all(), 3) = 1;
    
    this->_positions = convert_to(positions, units::m);
}

Operator
Model
::build_pulse(Quantity const & angle, Quantity const & phase) const
{
    return this->build_pulse(
        TensorQ<1>{angle}, TensorQ<1>{phase});
}

Operator
Model
::build_pulse(TensorQ<1> const & angle, TensorQ<1> const & phase) const
{
    if(
        angle.size() != phase.size()
        || (angle.size() != 1 && angle.size() != this->_positions.shape()[0]))
    {
        throw std::runtime_error("Size mismatch");
    }
    
    TensorR<1> const cos_angle = xt::cos(angle), cos_phase = xt::cos(phase);
    TensorR<1> const sin_angle = xt::sin(angle), sin_phase = xt::sin(phase);
    
    Operator::Array op = xt::zeros<Operator::Array::value_type>(
        Operator::Array::shape_type{angle.size(), 4, 4});
    for(std::size_t i=0; i<angle.size(); ++i)
    {
        auto const ca=cos_angle.unchecked(i), cp=cos_phase.unchecked(i);
        auto const sa=sin_angle.unchecked(i), sp=sin_phase.unchecked(i);
        auto const cp2=std::pow(cp, 2), sp2=std::pow(sp, 2);
        xt::view(op, i) = TensorR<2>{
            {sp2*ca - sp2 + 1, (1-ca)*sp*cp,  sa*sp, 0},
            {    (1-ca)*sp*cp,   sp2+ca*cp2, -sa*cp, 0},
            {          -sa*sp,        sa*cp,     ca, 0},
            {               0,            0,      0, 1}
        };
    }
    return {op};
}

Operator
Model
::build_time_interval(
    Quantity const & duration, Quantity const & delta_omega,
    TensorQ<1> const & gradient) const
{
    return this->build_time_interval(
        duration, TensorQ<1>{delta_omega}, xt::atleast_2d(gradient));
}

Operator
Model
::build_time_interval(
    Quantity const & duration, TensorQ<1> const & delta_omega,
    TensorQ<2> const & gradient) const
{
    auto const duration_s = duration.convert_to(units::s);
    auto const delta_omega_Hz = convert_to(delta_omega, units::Hz);
    auto const gradient_T_per_m = convert_to(gradient, units::T/units::m);
    
    auto const angular_frequency =
        2*M_PI * (
            // Field-related dephasing
            delta_omega_Hz
            // FIXME: delta_omega of species 
            // // Species-related dephasing, e.g. chemical shift or susceptibility
            // + species.delta_omega
        )
        // Gradient-related dephasing
        + gamma.magnitude*xt::sum(gradient_T_per_m * this->_positions, {1});
    
    auto op = this->build_phase_accumulation(
        units::rad * duration_s * angular_frequency);
    
    op.preMultiply(this->build_relaxation(duration));
    
    return op;
}

Operator
Model
::build_relaxation(Quantity const & duration) const
{
    Operator::Array op = xt::zeros<Operator::Array::value_type>(
        Operator::Array::shape_type{this->_positions.shape()[0], 4, 4});
    auto const duration_s = duration.convert_to(units::s);
    auto const E1 = xt::exp(-duration_s/this->_T1);
    auto const E2 = xt::exp(-duration_s/this->_T2);
    xt::view(op, xt::all(), 0, 0) = E2;
    xt::view(op, xt::all(), 1, 1) = E2;
    xt::view(op, xt::all(), 2, 2) = E1;
    xt::view(op, xt::all(), 3, 3) = 1;
    xt::view(op, xt::all(), 2, 3) = this->_M0*(1-E1);
    
    return {op};
}

Operator
Model
::build_phase_accumulation(Quantity const & angle) const
{
    return this->build_phase_accumulation(TensorQ<1>{angle});
}

Operator
Model
::build_phase_accumulation(TensorQ<1> const & angle) const
{
    TensorR<1> const cos_angle = xt::cos(angle);
    TensorR<1> const sin_angle = xt::sin(angle);
    
    Operator::Array op = xt::zeros<Operator::Array::value_type>(
        Operator::Array::shape_type{angle.size(), 4, 4});
    for(std::size_t i=0; i<angle.size(); ++i)
    {
        auto const ca=cos_angle.unchecked(i), sa=sin_angle.unchecked(i);
        xt::view(op, i) = TensorR<2>{
            {ca, -sa, 0, 0},
            {sa,  ca, 0, 0},
            { 0,   0, 1, 0},
            { 0,   0, 0, 1}
        };
    }
    return {op};
}

void
Model
::apply(Operator const & operator_)
{
    auto const & array = operator_.array();
    for(std::size_t n=0; n<this->_magnetization.shape()[0]; ++n)
    {
        auto l = xt::view(array, n);
        auto r = xt::view(this->_magnetization, n);
        auto temp = xt::zeros_like(r);
        for(std::size_t i=0; i<4; ++i)
        {
            temp.unchecked(i) += 
                l.unchecked(i,0)*r.unchecked(0)
                + l.unchecked(i,1)*r.unchecked(1)
                + l.unchecked(i, 2)*r.unchecked(2)
                + l.unchecked(i,3)*r.unchecked(3);
        }
        r = temp;
    }
}

TensorQ<1>
Model
::T1() const
{
    return this->_T1*units::s;
}

TensorQ<1>
Model
::T2() const
{
    return this->_T2*units::s;
}

TensorR<1> const &
Model
::M0() const
{
    return this->_M0;
}

TensorR<2>
Model
::magnetization() const
{
    return xt::eval(
        xt::view(this->_magnetization, xt::all(), xt::range(0, 3))
        / xt::expand_dims(xt::view(this->_magnetization, xt::all(), 3), 1));
}

TensorQ<2>
Model
::positions() const
{
    return this->_positions*units::m;
}

}

}
