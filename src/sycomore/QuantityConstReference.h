#ifndef _eb453df2_9d88_4a08_9015_ed0be16f00ec
#define _eb453df2_9d88_4a08_9015_ed0be16f00ec

#include "sycomore/Dimensions.h"
#include "sycomore/QuantityInterface.h"

namespace sycomore
{

class Quantity;

/// @brief Read-only view to a quantity
class QuantityConstReference: public QuantityInterface<QuantityConstReference>
{
public:
    double const & magnitude;
    Dimensions const & dimensions;
    
    QuantityConstReference(double const & magnitude, Dimensions const & dimensions)
    : magnitude(magnitude), dimensions(dimensions)
    {
        // Nothing else
    }
    
    QuantityConstReference(QuantityConstReference const &) = default;
    QuantityConstReference(QuantityConstReference &&) = default;
    QuantityConstReference & operator=(QuantityConstReference const &) = delete;
    QuantityConstReference & operator=(QuantityConstReference &&) = delete;
    ~QuantityConstReference() override = default;
};

template<>
struct OwningTypeTrait<QuantityConstReference> { using Type = Quantity; };


}

#endif // _eb453df2_9d88_4a08_9015_ed0be16f00ec
