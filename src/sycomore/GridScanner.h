#ifndef _06060b47_05d0_4c7a_a1bf_65e13262fb2a
#define _06060b47_05d0_4c7a_a1bf_65e13262fb2a

#include <utility>

#include "sycomore/sycomore.h"
#include "sycomore/sycomore_api.h"

namespace sycomore
{

class SYCOMORE_API GridScanner
{
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = std::pair<Index, std::size_t>;
    using difference_type = std::ptrdiff_t;
    using pointer = value_type const *;
    using reference = value_type const &;

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

    Index const & region_end() const;

private:
    value_type _iterator;

    Index _region_origin;
    Index _region_end;
    Stride _offset_increment;

    Index _end_index;
    std::size_t _end_offset;
};

}

#endif // _06060b47_05d0_4c7a_a1bf_65e13262fb2a
