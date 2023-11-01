#ifndef _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74
#define _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/units.h"

namespace sycomore
{

/// @brief Time interval, with or without magnetic field gradient.
class TimeInterval
{
public:
    // TODO: gradient shape
    
    /**
     * @brief Shortest possible time interval given 1D gradient area (T/m*s)
     * or dephasing (rad/m) and maximum gradient amplitude
     */
    static TimeInterval shortest(Quantity const & k, Quantity const & G_max);
    
    /**
     * @brief Shortest possible time interval given 3D gradient area (T/m*s)
     * or dephasing (rad/m) and maximum gradient amplitude
     */
    static TimeInterval shortest(Vector3Q const & k, Quantity const & G_max);
    
    /// @brief Default-initialize duration and gradient.
    TimeInterval();
    
    /**
     * @brief Constructor, gradient may be specified as amplitude (in T/m), 
     * area (in T/m*s) or dephasing (in rad/m).
     */
    TimeInterval(
        Quantity const & duration,
        Quantity const & gradient=0*units::T/units::m);
    
    /**
     * @brief Constructor, gradient may be specified as amplitude (in T/m), 
     * area (in T/m*s) or dephasing (in rad/m).
     */
    TimeInterval(Quantity const & duration, Vector3Q const & gradient);
    
    TimeInterval(TimeInterval const &) = default;
    TimeInterval(TimeInterval &&) = default;
    TimeInterval & operator=(TimeInterval const &) = default;
    TimeInterval & operator=(TimeInterval &&) = default;
    ~TimeInterval() = default;
    
    /// @brief Return the duration.
    Quantity const & duration() const;
    
    /// @brief Set the duration.
    void set_duration(Quantity const & q);
    
    /**
     * @brief Set gradient amplitude (in T/m), area (in T/m*s) or dephasing 
     * (in rad/m).
     */
    void set_gradient(Quantity const & q);
    
    /**
     * @brief Set gradient amplitude (in T/m), area (in T/m*s) or dephasing 
     * (in rad/m).
     */
    void set_gradient(Vector3Q const & a);

    /// @brief Return the gradient amplitude
    Vector3Q const & gradient_amplitude() const;
    
    /// @brief Return the gradient area
    Vector3Q gradient_area() const;
    
    /// @brief Return the gradient dephasing.
    Vector3Q gradient_dephasing() const;
    
    /// @brief Equality of duration and gradient
    bool operator==(TimeInterval const & other) const;
    
    /// @brief Difference of duration or gradient
    bool operator!=(TimeInterval const & other) const;

private:
    /// @brief Interval duration
    Quantity _duration;

    /// @brief Gradient amplitude in T/m on x,y,z axes.
    Vector3Q _gradient_amplitude;
};

}

#endif // _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74
