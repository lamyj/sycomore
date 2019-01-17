#ifndef _06060b47_05d0_4c7a_a1bf_65e13262fb2a
#define _06060b47_05d0_4c7a_a1bf_65e13262fb2a

#include <utility>

#include "sycomore/sycomore.h"

namespace sycomore
{

class GridScanner
{
public:
    using value_type = std::pair<Index, size_t>;

    GridScanner(Index const & origin, Shape const & shape);

    GridScanner(
        Index const & origin, Shape const & shape,
        Index const & region_origin, Shape const & region_shape);

    GridScanner begin() const;

    GridScanner end() const;

    value_type const & operator*() const;
    value_type const * operator->() const;

    bool operator==(GridScanner const & other) const;
    bool operator!=(GridScanner const & other) const;

    GridScanner const & operator++();
private:
    value_type _iterator;

    Index _region_origin;
    Index _region_end;
    Stride _offset_increment;

    Index _end_index;
    size_t _end_offset;
};

Stride compute_stride(Shape const & shape);

Stride
compute_offset_increment(
    Index const & region_origin, Shape const & region_shape,
    Index const & origin, Shape const & shape, Stride const & stride);

void increment(
    Index & index, size_t & offset,
    Index const & region_origin, Index const & region_end,
    Stride const & offset_increment);

}

#endif // _06060b47_05d0_4c7a_a1bf_65e13262fb2a
