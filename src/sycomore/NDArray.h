#ifndef _f40986ed_5367_42c0_a613_78150c4df2ae
#define _f40986ed_5367_42c0_a613_78150c4df2ae

#include <vector>

#include "sycomore/sycomore.h"

namespace sycomore
{

/// @brief N-dimensional array.
template<typename TScalar>
class NDArray
{
public:
    NDArray();
    NDArray(Shape const & shape);
    NDArray(Shape const & shape, TScalar const & value);

    NDArray(NDArray<TScalar> const &) = default;
    NDArray(NDArray<TScalar> &&) = default;

    NDArray<TScalar> & operator=(NDArray<TScalar> const &) = default;
    NDArray<TScalar> & operator=(NDArray<TScalar> &&) = default;

    ~NDArray() = default;

    TScalar const & operator[](Index const & index) const;
    TScalar & operator[](Index const & index);

    unsigned int dimension() const;
    Shape const & shape() const;
    Stride const & stride() const;

    void reshape(Shape const & shape);
    void reshape(Shape const & shape, TScalar const & value);

    TScalar const * data() const;
    TScalar * data();

    using iterator = typename std::vector<TScalar>::iterator;
    using const_iterator = typename std::vector<TScalar>::const_iterator;

    iterator begin();
    const_iterator begin() const;
    const_iterator cbegin() const;

    iterator end();
    const_iterator end() const;
    const_iterator cend() const;

private:
    Shape _shape;
    Stride _stride;
    std::vector<TScalar> _data;

    static Stride _compute_stride(Shape const & shape);
    void _reshape(NDArray<TScalar> & new_array);
};

}

#include "NDArray.txx"

#endif // _f40986ed_5367_42c0_a613_78150c4df2ae
