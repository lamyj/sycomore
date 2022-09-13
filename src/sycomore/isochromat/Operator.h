#ifndef _e0796018_5e39_4c59_988a_1e882e463fd4
#define _e0796018_5e39_4c59_988a_1e882e463fd4

#include <xtensor/xtensor.hpp>

#include "sycomore/sycomore.h"

namespace sycomore
{

namespace isochromat
{

/// @brief Isochromat simulation operator, i.e. an array of 4×4 matrices
class SYCOMORE_API Operator
{
public:
    using Array = xt::xtensor<Real, 3>;
    using value_type = Array::value_type;
    using shape_type = Array::shape_type;
    
    /// @brief Build an identity operator.
    Operator();
    
    /// @brief Build an operator from given array.
    Operator(xt::xtensor<Real, 3> const & data);
    
    /// @brief Default copy constructor.
    Operator(Operator const &) = default;
    
    /// @brief Default move constructor.
    Operator(Operator &&) = default;
    
    /// @brief Default destructor.
    ~Operator() = default;
    
    /// @brief Default assignment.
    Operator & operator=(Operator const &) = default;
    
    /// @brief Default move assignment.
    Operator & operator=(Operator &&) = default;
    
    /**
     * @brief In-place chaining of operators: *this will represent right follow
     * by *this.
     *
     * If *this and/or right has shape 1×4×4, it is broadcast to match the other
     * operand. Otherwise, both operand must have shape n×4×4.
     */
    Operator & operator*=(Operator const & right);
    
    /**
     * @brief In-place chaining of operators: *this will represent *this follow
     * by right.
     *
     * If *this and/or right has shape 1×4×4, it is broadcast to match the other
     * operand. Otherwise, both operand must have shape n×4×4.
     */
    Operator & preMultiply(Operator const & left);
    
    /// @brief Numeric representation of the operator.
    Array const & array() const;
private:
    Array _array;
};

/// @brief Operator chaining, representing right followed by left.
Operator operator*(Operator left, Operator const & right);

}

}

#endif // _e0796018_5e39_4c59_988a_1e882e463fd4
