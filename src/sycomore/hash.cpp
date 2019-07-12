#include "hash.h"

#include <cstddef>

namespace sycomore
{

void combine_hashes(std::size_t & seed, std::size_t value)
{
    seed ^= value + 0x9e3779b9 + (seed<<6) + (seed>>2);
}

}
