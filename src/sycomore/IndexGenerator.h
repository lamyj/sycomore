#ifndef _8b5ae53e_83af_4e5b_b359_e9552c756ce5
#define _8b5ae53e_83af_4e5b_b359_e9552c756ce5

#include "sycomore/sycomore.h"

namespace sycomore
{

/// @brief Generate all indices of a grid.
class IndexGenerator
{
public:
    /// @brief Constructor for zero-origin grids.
    IndexGenerator(Shape const & shape);

    /// @brief constructor.
    IndexGenerator(Index const & origin, Shape const & shape);

    class const_iterator
    {
    public:
        using value_type = Index;

        const_iterator(
            Index const & origin, Shape const & shape,
            Index const & last, Index const & end, bool is_end);

        bool operator==(const_iterator const & other) const;

        bool operator!=(const_iterator const & other) const;

        const_iterator & operator++();

        const_iterator & operator++(int);

        value_type const & operator*();
    private:
        Index _origin;
        Shape _shape;
        Index _last;
        Index _end;

        Index _index;
    };

    Index const & origin() const;
    Shape const & shape() const;

    /// @brief Return the one-past-the-end index on all axes.
    Index const & last() const;

    const_iterator begin() const;
    const_iterator end() const;

private:
    Index _origin;
    Shape _shape;

    Index _last;
    Index _end;
};

}

#endif // _8b5ae53e_83af_4e5b_b359_e9552c756ce5
