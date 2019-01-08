#include "Grid.h"

#include <vector>
#include "sycomore/magnetization.h"

namespace sycomore
{

Grid
::Grid()
: Grid(Index(), Shape())
{
    // Nothing else.
}

Grid
::Grid(Index const & origin, Shape const & shape)
: _origin(origin), _array(shape)
{
    // Nothing else.
}

Grid
::Grid(Index const & origin, Shape const & shape, value_type const & value)
: _origin(origin), _array(shape, value)
{
    // Nothing else.
}

Grid::value_type const *
Grid
::data() const
{
    return this->_array.data();
}

Grid::value_type *
Grid
::data()
{
    return this->_array.data();
}

Grid::value_type &
Grid
::operator[](Index const & index)
{
    Index translated_index(index.size());
    std::transform(
        index.begin(), index.end(), this->_origin.begin(),
        translated_index.begin(), std::minus<int>());
    return this->_array[translated_index];
}

Grid::value_type const &
Grid
::operator[](Index const & index) const
{
    Index translated_index(index.size());
    std::transform(
        index.begin(), index.end(), this->_origin.begin(),
        translated_index.begin(), std::minus<int>());
    return this->_array[translated_index];
}

size_t
Grid
::dimension() const
{
    return this->_array.dimension();
}

Index const &
Grid
::origin() const
{
    return this->_origin;
}

Shape const &
Grid
::shape() const
{
    return this->_array.shape();
}

Stride const &
Grid
::stride() const
{
    return this->_array.stride();
}

void
Grid
::reshape(Index const & origin, Shape const & shape)
{
    if(shape.size() != this->dimension())
    {
        throw std::runtime_error("Reshaping must preserve dimension");
    }

    Grid new_grid(origin, shape);
    this->_reshape(new_grid);
}

void
Grid
::reshape(Index const & origin, Shape const & shape, value_type const & value)
{
    if(shape.size() != this->dimension())
    {
        throw std::runtime_error("Reshaping must preserve dimension");
    }

    Grid new_grid(origin, shape, value);
    this->_reshape(new_grid);
}

Grid::iterator
Grid
::begin()
{
    return this->_array.begin();
}

Grid::const_iterator
Grid
::begin() const
{
    return this->_array.begin();
}

Grid::const_iterator
Grid
::cbegin() const
{
    return this->_array.cbegin();
}

Grid::iterator
Grid
::end()
{
    return this->_array.end();
}

Grid::const_iterator
Grid
::end() const
{
    return this->_array.end();
}

Grid::const_iterator
Grid
::cend() const
{
    return this->_array.cend();
}

void
Grid
::_reshape(Grid & new_grid)
{
    Index max_origin(this->_origin.size());
    std::transform(
        this->_origin.begin(), this->_origin.end(), new_grid.origin().begin(),
        max_origin.begin(), [](int x, int y) { return std::max(x, y); });

    auto const old_last = IndexGenerator(this->origin(), this->shape()).last();
    auto const new_last = IndexGenerator(
        new_grid.origin(), new_grid.shape()).last();
    Index min_last(this->_origin.size());
    std::transform(
        old_last.begin(), old_last.end(), new_last.begin(),
        min_last.begin(), [](int x, int y) { return std::min(x, y); });

    // Check whether the two grids are disjoint (on at least one axis,
    // max_origin >= min_last)
    auto mismatch = std::mismatch(
        max_origin.begin(), max_origin.end(), min_last.begin(),
        std::less<int>());
    if(mismatch.first == max_origin.end())
    {
        Shape min_shape(min_last.size());
        std::transform(
            min_last.begin(), min_last.end(),
            max_origin.begin(), min_shape.begin(),
            std::minus<int>());

        // NOTE: we could memcopy line-by-line
        for(auto && index: IndexGenerator(max_origin, min_shape))
        {
            new_grid[index] = (*this)[index];
        }
    }
    // Otherwise the two grids are disjoint, nothing to copy.

    this->_origin = std::move(new_grid._origin);
    this->_array = std::move(new_grid._array);
}

}
