#ifndef _e0796018_5e39_4c59_988a_1e882e463fd4
#define _e0796018_5e39_4c59_988a_1e882e463fd4

#include <xtensor/xtensor.hpp>

#include "sycomore/sycomore.h"

namespace sycomore
{

namespace isochromat
{

/// @brief Isochromat simulation operator, i.e. an array of 4×4 matrices
class Operator
{
public:
    /// @brief Array representation of the operator
    using Array = TensorR<3>;
    
    /// @brief Build an identity operator.
    Operator();
    
    /// @brief Build an operator from given array.
    Operator(Array const & data);
    
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
     * @brief In-place chaining of operators, with right applied first.
     *
     * If *this and/or right has shape 1×4×4, it is broadcast to match the other
     * operand. Otherwise, both operand must have shape n×4×4.
     */
    Operator & operator*=(Operator const & right);
    
    /**
     * @brief In-place chaining of operators, with self applied first.
     *
     * If *this and/or left has shape 1×4×4, it is broadcast to match the other
     * operand. Otherwise, both operand must have shape n×4×4.
     */
    Operator & pre_multiply(Operator const & left);
    
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
