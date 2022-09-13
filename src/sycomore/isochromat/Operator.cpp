#include "Operator.h"

#include <stdexcept>
#include <xtensor/xbuilder.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xfixed.hpp>
#include <xtensor/xview.hpp>

#include "sycomore/sycomore.h"

namespace sycomore
{

namespace isochromat
{

Operator
::Operator()
: _array(xt::eye<Real>({1,4,4}))
{
    // Nothing else.
}

Operator
::Operator(xt::xtensor<Real, 3> const & data)
: _array(data)
{
    if(this->_array.shape()[1] != 4 || this->_array.shape()[2] != 4)
    {
        throw std::runtime_error("Invalid shape");
    }
}

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

Operator &
Operator
::operator*=(Operator const & right)
{
    static thread_local xt::xtensor_fixed<Real, xt::xshape<4, 4>> matrix;
    
    if(this->_array.shape()[0] == 1 && right._array.shape()[0] != 1)
    {
        auto temp = xt::empty_like(right._array);
        auto l = xt::view(this->_array, 0);
        for(std::size_t item=0, end=right._array.shape()[0]; item<end; ++item)
        {
            auto destination = xt::view(temp, item);
            auto r = xt::view(right._array, item);
            MULT_4_4(l, r, destination)
        }
        this->_array = std::move(temp);
    }
    else if(this->_array.shape()[0] != 1 && right._array.shape()[0] == 1)
    {
        auto r = xt::view(right._array, 0);
        for(std::size_t item=0, end=this->_array.shape()[0]; item<end; ++item)
        {
            auto l = xt::view(this->_array, item);
            MULT_4_4(l, r, matrix)
            l = matrix;
        }
    }
    else if(this->_array.shape()[0] == right._array.shape()[0])
    {
        for(std::size_t item=0, end=right._array.shape()[0]; item<end; ++item)
        {
            auto l = xt::view(this->_array, item);
            auto r = xt::view(right._array, item);
            MULT_4_4(l, r, matrix)
            l = matrix;
        }
    }
    else
    {
        throw std::runtime_error("Size mismatch");
    }
    
    return *this;
}

Operator &
Operator
::preMultiply(Operator const & left)
{
    static thread_local xt::xtensor_fixed<Real, xt::xshape<4, 4>> matrix;
    
    if(left._array.shape()[0] == 1 && this->_array.shape()[0] != 1)
    {
        auto l = xt::view(left._array, 0);
        for(std::size_t item=0, end=this->_array.shape()[0]; item<end; ++item)
        {
            auto r = xt::view(this->_array, item);
            MULT_4_4(l, r, matrix)
            r = matrix;
        }
    }
    else if(left._array.shape()[0] != 1 && this->_array.shape()[0] == 1)
    {
        auto temp = xt::empty_like(left._array);
        auto r = xt::view(this->_array, 0);
        for(std::size_t item=0, end=left._array.shape()[0]; item<end; ++item)
        {
            auto destination = xt::view(temp, item);
            auto l = xt::view(left._array, item);
            MULT_4_4(l, r, destination)
        }
        std::swap(this->_array, temp);
    }
    else if(left._array.shape()[0] == this->_array.shape()[0])
    {
        for(std::size_t item=0, end=left._array.shape()[0]; item<end; ++item)
        {
            auto l = xt::view(left._array, item);
            auto r = xt::view(this->_array, item);
            MULT_4_4(l, r, matrix)
            r = matrix;
        }
    }
    else
    {
        throw std::runtime_error("Size mismatch");
    }
    
    return *this;
}

#undef MULT_4_4

Operator::Array const &
Operator
::array() const
{
    return this->_array;
}

Operator operator*(Operator left, Operator const & right)
{
    left *= right;
    return left;
}

}

}
