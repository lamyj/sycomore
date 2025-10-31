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

#pragma warning(push)
#pragma warning(disable: 4267)
Model
::Model(
    Quantity const & T1, Quantity const & T2, TensorR<1> const & M0,
    TensorQ<2> const & positions, Quantity const & delta_omega)
: Model(
    TensorQ<1>(
        xt::repeat(TensorR<1>{T1.magnitude}, positions.shape()[0], 0),
        T1.dimensions),
    TensorQ<1>(
        xt::repeat(TensorR<1>{T2.magnitude}, positions.shape()[0], 0),
        T2.dimensions),
    xt::eval(xt::repeat(xt::atleast_2d(M0), positions.shape()[0], 0)),
    positions,
    TensorQ<1>(
        xt::repeat(TensorR<1>{delta_omega.magnitude}, positions.shape()[0], 0),
        delta_omega.dimensions))
{
    // Nothing else
}
#pragma warning(pop)

Model
::Model(
    TensorQ<1> const & T1, TensorQ<1> const & T2, TensorR<2> const & M0,
    TensorQ<2> const & positions, TensorQ<1> const & delta_omega)
: _T1(T1.shape()), _T2(T2.shape()), _M0(xt::view(M0, xt::all(), 2UL)),
    _delta_omega(delta_omega.size() == 0 ? T1.shape() : delta_omega.shape()),
    _magnetization(TensorR<2>::shape_type{M0.shape()[0], 4}),
    _positions(positions.shape())
{
    if(
        T1.size() != T2.size() || T1.size() != M0.shape()[0]
        || T1.size() != positions.shape()[0]
        || (delta_omega.size() != 0 && T1.size() != delta_omega.size()))
    {
        throw std::runtime_error("Size mismatch");
    }
    this->_T1 = T1.convert_to(units::s);
    this->_T2 = T2.convert_to(units::s);
    
    xt::view(this->_magnetization, xt::all(), xt::range(0, 3)) = M0;
    xt::view(this->_magnetization, xt::all(), 3UL) = 1;
    
    this->_delta_omega = 
        delta_omega.size() == 0
        ? xt::repeat(TensorR<1>{0.}, T1.size(), 0)
        : delta_omega.convert_to(units::Hz);
    
    this->_positions = positions.convert_to(units::m);
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
    TensorR<1> phase_dummy;
    if(phase.size() == 0)
    {
        phase_dummy = xt::repeat(TensorR<1>{0}, angle.size(), 0);
    }
    
    auto const & phase_ = (phase.size()>0?phase.magnitude:phase_dummy);
    
    if(
        angle.size() != phase_.size()
        || (angle.size() != 1 && angle.size() != this->_positions.shape()[0]))
    {
        throw std::runtime_error("Size mismatch");
    }
    
    auto const & angle_ = angle.magnitude;
    
    TensorR<1> const cos_angle = xt::cos(angle_), cos_phase = xt::cos(phase_);
    TensorR<1> const sin_angle = xt::sin(angle_), sin_phase = xt::sin(phase_);
    
    Operator::Array op = xt::zeros<Operator::Array::value_type>(
        Operator::Array::shape_type{angle.size(), 4, 4});
    for(std::size_t i=0; i<angle.size(); ++i)
    {
        // Using the Rodrigues' formula in matrix notation, the RF pulse is a
        // rotation of angle α around an axis (cos ϕ, sin ϕ, 0)
        //           0       0    sin ϕ
        // K =       0       0   -cos ϕ
        //      -sin ϕ   cos ϕ        0
        // And its square
        //           -sin² ϕ   sin ϕ⋅cos ϕ   0
        // K² =  sin ϕ⋅cos ϕ       -cos² ϕ   0
        //                 0             0  -1
        // The rotation matrix is then R = I + sin α K + (1-cos α) K²
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
        duration, TensorQ<1>{delta_omega},
        gradient.size() != 0 ?
        TensorQ<2>(xt::atleast_2d(gradient.magnitude), gradient.dimensions)
        : TensorQ<2>{});
}

Operator
Model
::build_time_interval(
    Quantity const & duration, TensorQ<1> const & delta_omega,
    TensorQ<2> const & gradient) const
{
    auto const & duration_ = duration.magnitude;
    auto const & delta_omega_ = delta_omega.magnitude;
    auto const & gradient_ = gradient.magnitude;
    
    TensorR<1> angular_frequency = 
        2*M_PI * (
            // Field-related dephasing
            delta_omega_
            // Species-related dephasing, e.g. chemical shift or susceptibility
            + this->_delta_omega);
    if(gradient.size() > 0)
    {
        angular_frequency += gamma.magnitude * xt::sum(
            gradient_ * this->_positions, {1});
    }
    
    auto op = this->build_relaxation(duration);
    
    op.pre_multiply(
        this->build_phase_accumulation(
            {duration_ * angular_frequency, Angle}));
    
    return op;
}

Operator
Model
::build_relaxation(Quantity const & duration) const
{
    Operator::Array op = xt::zeros<Operator::Array::value_type>(
        Operator::Array::shape_type{this->_positions.shape()[0], 4, 4});
    auto const & duration_ = duration.magnitude;
    auto const E1 = xt::exp(-duration_/this->_T1);
    auto const E2 = xt::exp(-duration_/this->_T2);
    xt::view(op, xt::all(), 0UL, 0UL) = E2;
    xt::view(op, xt::all(), 1UL, 1UL) = E2;
    xt::view(op, xt::all(), 2UL, 2UL) = E1;
    xt::view(op, xt::all(), 3UL, 3UL) = 1;
    xt::view(op, xt::all(), 2UL, 3UL) = this->_M0*(1-E1);
    
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
    auto const angle_ = angle.magnitude;
    TensorR<1> const cos_angle = xt::cos(angle_);
    TensorR<1> const sin_angle = xt::sin(angle_);
    
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
        auto temp = xt::empty_like(r);
        for(std::size_t i=0; i<4; ++i)
        {
            temp.unchecked(i) = 
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

TensorR<1> const &
Model
::delta_omega() const
{
    return this->_delta_omega;
}

TensorR<2>
Model
::magnetization() const
{
    return xt::eval(
        xt::view(this->_magnetization, xt::all(), xt::range(0, 3))
        / xt::expand_dims(xt::view(this->_magnetization, xt::all(), 3UL), 1));
}

TensorQ<2>
Model
::positions() const
{
    return this->_positions*units::m;
}

}

}
