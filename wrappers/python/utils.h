#ifndef _492c1887_a16a_4409_a300_9f2d70ba7a3b
#define _492c1887_a16a_4409_a300_9f2d70ba7a3b

#include <pybind11/pybind11.h>

template<typename T>
bool isinstance(pybind11::object o)
{
    try
    {
        o.cast<T>();
    }
    catch(pybind11::cast_error)
    {
        return false;
    }
    return true;
}

#endif // _492c1887_a16a_4409_a300_9f2d70ba7a3b
