#ifndef _846205cf_0438_473b_b1d7_29d36406c34a
#define _846205cf_0438_473b_b1d7_29d36406c34a

#include "sycomore/Dimensions.h"
#include "sycomore/QuantityInterface.h"

namespace sycomore
{

class Quantity;

/// @brief Read/write view to a quantity
class QuantityReference: public QuantityInterface<QuantityReference>
{
public:
    double & magnitude;
    Dimensions & dimensions;
    
    QuantityReference(double & magnitude, Dimensions & dimensions)
    : magnitude(magnitude), dimensions(dimensions)
    {
        // Nothing else
    }
    
    QuantityReference(QuantityReference const & other) = default;
    QuantityReference(QuantityReference && other) = default;
    QuantityReference & operator=(Quantity const & other);
    QuantityReference & operator=(QuantityReference && other) = delete;
    ~QuantityReference() override = default;
};

template<>
struct OwningTypeTrait<QuantityReference> { using Type = Quantity; };

}

#endif // _846205cf_0438_473b_b1d7_29d36406c34a
