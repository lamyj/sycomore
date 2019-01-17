#include "GridScanner.h"

#include <algorithm>
#include <utility>

#include "sycomore/sycomore.h"

namespace sycomore
{

GridScanner
::GridScanner(Index const & origin, Shape const & shape)
: GridScanner(origin, shape, origin, shape)
{
    // Nothing else.
}

GridScanner
::GridScanner(
    Index const & origin, Shape const & shape,
    Index const & region_origin, Shape const & region_shape)
: _region_origin(region_origin), _region_end(region_origin.size())
{
    auto const stride(compute_stride(shape));

    this->_iterator.first = region_origin;
    this->_iterator.second = dot(region_origin-origin, stride);
    std::transform(
        region_origin.begin(), region_origin.end(), region_shape.begin(),
        this->_region_end.begin(), std::plus<int>());

    this->_offset_increment = compute_offset_increment(
        region_origin, region_shape, origin, shape, stride);

    this->_end_index = region_origin;
    this->_end_index[this->_end_index.size()-1] =
        region_origin[this->_end_index.size()-1]
        +region_shape[this->_end_index.size()-1];
    this->_end_offset = dot(this->_end_index-origin, stride);
}

GridScanner
GridScanner
::begin() const
{
    return *this;
}

GridScanner
GridScanner
::end() const
{
    auto result(*this);
    result._iterator = std::make_pair(result._end_index, result._end_offset);
    return result;
}

GridScanner::value_type const &
GridScanner
::operator*() const
{
    return this->_iterator;
}

GridScanner::value_type const *
GridScanner
::operator->() const
{
    return &this->_iterator;
}

bool
GridScanner
::operator==(GridScanner const & other) const
{
    return this->_iterator.second == other._iterator.second;
}

bool
GridScanner
::operator!=(GridScanner const & other) const
{
    return !this->operator==(other);
}

GridScanner const &
GridScanner
::operator++()
{
    increment(
        this->_iterator.first, this->_iterator.second,
        this->_region_origin, this->_region_end, this->_offset_increment);
    return *this;
}

Stride compute_stride(Shape const & shape)
{
    if(shape.size() == 0)
    {
        return Stride();
    }

    Stride stride(shape.size()+1);
    stride[0] = 1;
    for(unsigned int i=0; i< shape.size(); ++i)
    {
        stride[i+1] = shape[i]*stride[i];
    }

    return stride;
}

Stride
compute_offset_increment(
    Index const & region_origin, Shape const & region_shape,
    Index const & origin, Shape const & shape, Stride const & stride)
{
    Stride result(stride.size());
    for(std::size_t i=0; i<result.size(); ++i)
    {
        result[i] = (
                (origin[i]+shape[i])-(region_origin[i]+region_shape[i])
                + region_origin[i]-origin[i]
            )*stride[i];
    }
    return result;
}

void increment(
    Index & index, size_t & offset,
    Index const & region_origin, Index const & region_end,
    Stride const & offset_increment)
{
    ++offset;

    for(size_t i=0, end(index.size()); i<end; ++i)
    {
        ++index[i];
        if(index[i] == region_end[i] && i!=end-1)
        {
            index[i] = region_origin[i];
            offset += offset_increment[i];
        }
        else
        {
            break;
        }
    }
}

}
