#include "Model.h"

#include <array>
#include <cmath>
#include <stdexcept>
#include <xtensor/xfixed.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

namespace sycomore
{

namespace isochromat
{

Model
::Model(Real T1, Real T2, Magnetization M0, Positions const & positions)
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
::Model(Array T1, Array T2, Magnetizations M0, Positions const & positions)
: _T1(T1), _T2(T2), _M0(M0), _positions(positions)
{
    if(
        T1.size() != T2.size() || T1.size() != M0.size()
        || T1.size() != positions.shape()[0])
    {
        throw std::runtime_error("Size mismatch");
    }
}

Operator
Model
::build_pulse(Real angle, Real phase) const
{
    return this->build_pulse(Array{angle}, Array{phase});
}

Operator
Model
::build_pulse(Array const & angle, Array const & phase) const
{
    if(
        angle.size() != phase.size()
        || angle.size() != 1 && angle.size() != this->_positions.shape()[0])
    {
        throw std::runtime_error("Size mismatch");
    }
    
    Operator op = xt::zeros<Operator::value_type>(
        Operator::shape_type{angle.size(), 4, 4});
    for(std::size_t i=0; i<angle.size(); ++i)
    {
        auto const ca=std::cos(angle[i]), cp=std::cos(phase[i]);
        auto const sa=std::sin(angle[i]), sp=std::sin(phase[i]);
        auto const cp2=std::pow(cp, 2), sp2=std::pow(sp, 2);
        xt::view(op, i) = Matrix{
            {sp2*ca - sp2 + 1, (1-ca)*sp*cp,  sa*sp, 0},
            {    (1-ca)*sp*cp,   sp2+ca*cp2, -sa*cp, 0},
            {          -sa*sp,        sa*cp,     ca, 0},
            {               0,            0,      0, 1}
        };
    }
    return op;
}

void
Model
::build_pulse(Real angle, Real phase, Operator & op) const
{
    this->_combine(this->build_pulse(Array{angle}, Array{phase}), op);
}

void
Model
::build_pulse(
    Array const & angle, Array const & phase, Operator & op) const
{
    this->_combine(this->build_pulse(angle, phase), op);
}

void
Model
::_combine(Operator const & left, Operator & right) const
{
    static thread_local Matrix matrix;
    static thread_local Operator array = xt::empty<Operator::value_type>(
        Operator::shape_type{this->_positions.shape()[0], 4, 4});

    #define MULT_4_4(l, r, dest) \
        for(std::size_t i=0; i<4; ++i) \
        { \
            for(std::size_t j=0; j<4; ++j) \
            { \
                dest.unchecked(i,j) =  \
                    l.unchecked(i,0)*r.unchecked(0,j) \
                    + l.unchecked(i,1)*r.unchecked(1,j) \
                    + l.unchecked(i, 2)*r.unchecked(2,j) \
                    + l.unchecked(i,3)*r.unchecked(3,j); \
            } \
        }
    
    auto const positions_count = this->_positions.shape()[0];
    
    if(left.shape()[0] == 1 && right.shape()[0] == positions_count)
    {
        auto l = xt::view(left, 0);
        for(std::size_t item=0, end=right.shape()[0]; item<end; ++item)
        {
            auto r = xt::view(right, item);
            MULT_4_4(l, r, matrix)
            r = matrix;
        }
    }
    else if(left.shape()[0] == positions_count && right.shape()[0] == 1)
    {
        auto r = xt::view(right, 1);
        for(std::size_t item=0, end=left.shape()[0]; item<end; ++item)
        {
            auto destination = xt::view(array, item);
            auto l = xt::view(left, item);
            MULT_4_4(l, r, destination)
        }
        std::swap(right, array);
    }
    else if(left.shape()[0] == right.shape()[0])
    {
        for(std::size_t item=0, end=left.shape()[0]; item<end; ++item)
        {
            auto l = xt::view(left, item);
            auto r = xt::view(right, item);
            MULT_4_4(l, r, matrix)
            r = matrix;
        }
    }
    else
    {
        throw std::runtime_error("Size mismatch");
    }
    
    #undef MULT_4_4
}

}

}
