# Time Intervals

Time intervals are defined by their duration and an optional gradient moment. A time interval can be created using either [units](units.md) or raw values: the duration is a time, with a raw value specified in seconds, and the gradient moment is a 3D vector of _1/Length_, with a raw value specified in _1/m_. This latter quantity comes from the fact that the gradient moment is the time integral of the gradient amplitude (in _[T/m]_) multiplied by the gyromagnetic ratio (in _[1/(T*s)]_).

The following code sample creates two equivalent time intervals:

```cpp
#include <sycomore/TimeInterval.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    // Values in base SI values
    sycomore::TimeInterval const idle(0.5, 0.001);
    
    // Values using units
    sycomore::TimeInterval const refocalization(500_ms, 1/mm);
}
```
