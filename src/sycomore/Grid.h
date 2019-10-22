#ifndef _e22643b8_c194_4aa9_99f9_6effee6222c8
#define _e22643b8_c194_4aa9_99f9_6effee6222c8

#include <vector>

#include "sycomore/sycomore.h"

namespace sycomore
{

/// @brief Discrete, n-dimensional grid of data.
template<typename T>
class Grid
{
public:
    using value_type = T;

    Grid();
    Grid(Index const & origin, Shape const & shape);
    Grid(Index const & origin, Shape const & shape, value_type const & value);

    Grid(Grid const &) = default;
    Grid(Grid &&) = default;

    Grid<T> & operator=(Grid const &) = default;
    Grid<T> & operator=(Grid &&) = default;

    ~Grid() = default;

    value_type const * data() const;
    value_type * data();

    value_type & operator[](std::size_t offset);
    value_type const & operator[](std::size_t offset) const;

    value_type & operator[](Index const & index);
    value_type const & operator[](Index const & index) const;

    std::size_t dimension() const;
    Index const & origin() const;
    Shape const & shape() const;
    Stride const & stride() const;

    void reshape(Index const & origin, Shape const & shape);
    void reshape(
        Index const & origin, Shape const & shape, value_type const & value);

    using iterator = typename std::vector<value_type>::iterator;
    using const_iterator = typename std::vector<value_type>::const_iterator;

    iterator begin();
    const_iterator begin() const;
    const_iterator cbegin() const;

    iterator end();
    const_iterator end() const;
    const_iterator cend() const;

private:
    Index _origin;
    Shape _shape;
    Stride _stride;
    std::vector<value_type> _data;

    static Stride _compute_stride(Shape const & shape);
    void _allocate(Shape);
    void _reshape(Grid & new_grid);
};

}

#include "Grid.txx"

#endif // _e22643b8_c194_4aa9_99f9_6effee6222c8
