#include "QuantityReference.h"

#include "sycomore/Dimensions.h"
#include "sycomore/Quantity.h"

namespace sycomore
{

QuantityReference &
QuantityReference
::operator=(Quantity const & other)
{
    other.check_dimensions(*this, "Assignment requires same dimensions");
    this->magnitude = other.magnitude;
    
    return *this;
}

QuantityReference
::operator Quantity() const
{
    return {this->magnitude, this->dimensions};
}

}
