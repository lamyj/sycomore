#ifndef _5b699f7b_cbfd_44de_ba50_c1fa8a69f228
#define _5b699f7b_cbfd_44de_ba50_c1fa8a69f228

#include <cstddef>
#include <type_traits>
#include <xsimd/xsimd.hpp>

namespace sycomore
{

/**
 * @brief Low-level container of simple types.
 *
 * Buffer is designed to be more efficient than std::vector by not initializing
 * its contents.
 */
template<
    typename T,
    typename Enable=std::enable_if_t<std::is_standard_layout<T>::value>>
class Buffer
{
public:
    using Self = Buffer<T>;
    using value_type = T;
    using size_type = std::size_t;
    using const_iterator = T const *;
    using iterator = T *;
    
    /// @brief Construct Buffer of given size.
    Buffer(size_type size=0);
    
    /// @brief Copy constructor.
    Buffer(Self const & other);
    
    /// @brief Move constructor.
    Buffer(Self && other);
    
    /// @brief Construct buffer with given values.
    Buffer(std::initializer_list<T> l);
    
    /// @brief Copy assignement.
    Self & operator=(Self const & other);
    
    /// @brief Move assignment.
    Self & operator=(Self && other);
    
    /// @brief Destructor.
    ~Buffer();
    
    /// @brief Number of elements in the buffer.
    size_type size() const;
    
    /// @brief Max number of elements in the buffer before reallocation.
    size_type capacity() const;
    
    /// @brief Pointer to the first element.
    T * data();
    
    /// @brief Pointer to the first element.
    T const * data() const;
    
    /// @brief Access an existing element.
    T const & operator[](std::size_t index) const;
    
    /// @brief Access an existing element.
    T & operator[](std::size_t index);
    
    /// @brief Change the number of elements, reallocate if capacity is exceeded.
    void resize(std::size_t size);
    
    /**
     * @brief Change the number of elements, reallocate if capacity is exceeded,
     *        set newly-created elements to value.
     */
    void resize(std::size_t size, T && x);
    
    /// @brief Swap contents with other buffer.
    void swap(Self & other);
    
    /// @brief Const iterator to the first element.
    const_iterator begin() const;
    /// @brief Const iterator to the first element.
    const_iterator cbegin() const;
    /// @brief Const iterator after the last element.
    const_iterator end() const;
    /// @brief Const iterator after the last element.
    const_iterator cend() const;
    
    /// @brief Iterator to the first element.
    iterator begin();
    /// @brief Iterator after the last element.
    iterator end();
    
private:
    using allocator = xsimd::aligned_allocator<T, 64>;
    allocator _allocator;
    T * _data;
    size_type _capacity, _size;
    
    void _deallocate();
    
    void _allocate(std::size_t n);
};

/// @brief Swap two buffers.
template<typename T>
void swap(Buffer<T> & b1, Buffer<T> & b2);

}

#include "Buffer.txx"

#endif // _5b699f7b_cbfd_44de_ba50_c1fa8a69f228
