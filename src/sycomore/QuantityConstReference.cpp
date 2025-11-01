#include "QuantityConstReference.h"

#include "sycomore/Quantity.h"

namespace sycomore
{

QuantityConstReference
::operator Quantity() const
{
    return {this->magnitude, this->dimensions};
}

}
