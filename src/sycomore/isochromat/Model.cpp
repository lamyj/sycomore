#include "Model.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <xtensor/xmath.hpp>
#include <xtensor/xmanipulation.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

namespace sycomore
{

namespace isochromat
{

Model
::Model(
    Real T1, Real T2,
    xt::xtensor<Real, 1> const & M0, xt::xtensor<Real, 2> const & positions)
: _T1{std::array<std::size_t, 1>{1}}, _T2{std::array<std::size_t, 1>{1}},
    _M0{std::array<std::size_t, 1>{1}},
    _magnetization{std::array<std::size_t, 2>{positions.shape()[0], 4}},
    _positions(positions)
{
    this->_T1 = {T1};
    this->_T2 = {T2};
    this->_M0 = {M0[2]};
    xt::col(this->_magnetization, 0).fill(M0[0]);
    xt::col(this->_magnetization, 1).fill(M0[1]);
    xt::col(this->_magnetization, 2).fill(M0[2]);
    xt::col(this->_magnetization, 3).fill(1.);
}

Model
::Model(
    xt::xtensor<Real, 1> const & T1, xt::xtensor<Real, 1> const & T2,
    xt::xtensor<Real, 2> const & M0, xt::xtensor<Real, 2> const & positions)
: _T1(T1), _T2(T2), _M0(xt::view(M0, xt::all(), 2)), _magnetization(M0),
    _positions(positions)
{
    if(
        T1.size() != T2.size() || T1.size() != M0.shape()[0]
        || T1.size() != positions.shape()[0])
    {
        throw std::runtime_error("Size mismatch");
    }
}

Operator
Model
::build_pulse(Real angle, Real phase) const
{
    return this->build_pulse(
        xt::xtensor<Real, 1>{angle}, xt::xtensor<Real, 1>{phase});
}

Operator
Model
::build_pulse(
    xt::xtensor<Real, 1> const & angle, xt::xtensor<Real, 1> const & phase) const
{
    if(
        angle.size() != phase.size()
        || angle.size() != 1 && angle.size() != this->_positions.shape()[0])
    {
        throw std::runtime_error("Size mismatch");
    }
    
    Operator::Array op = xt::zeros<Operator::value_type>(
        Operator::shape_type{angle.size(), 4, 4});
    for(std::size_t i=0; i<angle.size(); ++i)
    {
        auto const ca=std::cos(angle[i]), cp=std::cos(phase[i]);
        auto const sa=std::sin(angle[i]), sp=std::sin(phase[i]);
        auto const cp2=std::pow(cp, 2), sp2=std::pow(sp, 2);
        xt::view(op, i) = xt::xtensor<Real, 2>{
            {sp2*ca - sp2 + 1, (1-ca)*sp*cp,  sa*sp, 0},
            {    (1-ca)*sp*cp,   sp2+ca*cp2, -sa*cp, 0},
            {          -sa*sp,        sa*cp,     ca, 0},
            {               0,            0,      0, 1}
        };
    }
    return op;
}

Operator
Model
::build_time_interval(
    Real duration, Real delta_omega, xt::xtensor<Real, 1> const & gradient) const
{
    return this->build_time_interval(
        duration, xt::xtensor<Real, 1>{delta_omega}, xt::atleast_2d(gradient));
}

Operator
Model
::build_time_interval(
    Real duration, xt::xtensor<Real, 1> const & delta_omega,
    xt::xtensor<Real, 2> const & gradient) const
{
    auto const angular_frequency =
        2*M_PI * (
            // Field-related dephasing
            delta_omega
            // FIXME: delta_omega of species 
            // // Species-related dephasing, e.g. chemical shift or susceptibility
            // + species.delta_omega
        )
        // Gradient-related dephasing
        + gamma.magnitude*xt::sum(gradient * this->_positions, {1});
    auto op = this->build_phase_accumulation(duration * angular_frequency);
    op.preMultiply(this->build_relaxation(duration));
    return op;
}

Operator
Model
::build_relaxation(Real duration) const
{
    Operator::Array op = xt::zeros<Operator::value_type>(
        Operator::shape_type{this->_positions.shape()[0], 4, 4});
    auto const E1 = xt::exp(-duration/this->_T1);
    auto const E2 = xt::exp(-duration/this->_T2);
    xt::view(op, xt::all(), 0, 0) = E2;
    xt::view(op, xt::all(), 1, 1) = E2;
    xt::view(op, xt::all(), 2, 2) = E1;
    xt::view(op, xt::all(), 3, 3) = 1;
    xt::view(op, xt::all(), 2, 3) = this->_M0*(1-E1);
    return op;
}

Operator
Model
::build_phase_accumulation(Real angle) const
{
    return this->build_phase_accumulation(xt::xtensor<Real, 1>{angle});
}

Operator
Model
::build_phase_accumulation(xt::xtensor<Real, 1> const & angle) const
{
    Operator::Array op = xt::zeros<Operator::value_type>(
        Operator::shape_type{angle.size(), 4, 4});
    for(std::size_t i=0; i<angle.size(); ++i)
    {
        auto const ca=std::cos(angle[i]), sa=std::sin(angle[i]);
        xt::view(op, i) = xt::xtensor<Real, 2>{
            {ca, -sa, 0, 0},
            {sa,  ca, 0, 0},
            { 0,   0, 1, 0},
            { 0,   0, 0, 1}
        };
    }
    return op;
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

xt::xtensor<Real, 1> const &
Model
::T1() const
{
    return this->_T1;
}

xt::xtensor<Real, 1> const &
Model
::T2() const
{
    return this->_T2;
}

xt::xtensor<Real, 1> const &
Model
::M0() const
{
    return this->_M0;
}

xt::xtensor<Real, 2> const &
Model
::magnetization() const
{
    return this->_magnetization;
}

xt::xtensor<Real, 2> const &
Model
::positions() const
{
    return this->_positions;
}

}

}
