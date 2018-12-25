#ifndef _f40986ed_5367_42c0_a613_78150c4df2ae
#define _f40986ed_5367_42c0_a613_78150c4df2ae

#include <vector>

namespace sycomore
{

using Shape = std::vector<unsigned int>;
using Stride = Shape;
using Index = Shape;

class IndexGenerator
{
public:
    IndexGenerator(Shape const & shape);

    class const_iterator
    {
    public:
        using value_type = Index;

        const_iterator(Shape const & shape, bool is_end);

        bool operator==(const_iterator const & other) const;

        bool operator!=(const_iterator const & other) const;

        const_iterator & operator++();

        const_iterator & operator++(int);

        value_type const & operator*();
    private:
        Shape _shape;
        Index _index;
        Index _end;

        void _next(Index & index);
    };

    const_iterator begin() const;

    const_iterator end() const;
private:
    Shape _shape;
};

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
};

}

#include "Array.txx"

#endif // _f40986ed_5367_42c0_a613_78150c4df2ae
