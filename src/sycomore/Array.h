#ifndef _f40986ed_5367_42c0_a613_78150c4df2ae
#define _f40986ed_5367_42c0_a613_78150c4df2ae

#include <vector>

#include "sycomore/sycomore.h"

namespace sycomore
{

/// @brief N-dimensional array.
template<typename TScalar>
class Array
{
public:
    Array();
    Array(Shape const & shape);
    Array(Shape const & shape, TScalar const & value);

    Array(Array<TScalar> const &) = default;
    Array(Array<TScalar> &&) = default;

    Array<TScalar> & operator=(Array<TScalar> const &) = default;
    Array<TScalar> & operator=(Array<TScalar> &&) = default;

    ~Array() = default;

    TScalar const & operator[](Index const & index) const;
    TScalar & operator[](Index const & index);

    unsigned int dimension() const;
    Shape const & shape() const;
    Stride const & stride() const;

    void reshape(Shape const & shape);
    void reshape(Shape const & shape, TScalar const & value);

    TScalar const * data() const;
    TScalar * data();

private:
    Shape _shape;
    Stride _stride;
    std::vector<TScalar> _data;

    static Stride _compute_stride(Shape const & shape);
    void _reshape(Array<TScalar> & new_array);
};

}

#include "Array.txx"

#endif // _f40986ed_5367_42c0_a613_78150c4df2ae
