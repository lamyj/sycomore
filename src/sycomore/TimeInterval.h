#ifndef _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74
#define _b442b891_b0a6_4ecd_9de5_e4cfc7e52a74

#include "sycomore/Array.h"
#include "sycomore/Quantity.h"
#include "sycomore/sycomore_api.h"
#include "sycomore/units.h"

namespace sycomore
{

/// @brief Time interval, with or without magnetic field gradient.
class SYCOMORE_API TimeInterval
{
public:
    // TODO: gradient shape
    
    /**
     * @brief Shortest possible time interval given 1D gradient moment and
     * maximum gradient amplitude
     */
    static TimeInterval shortest(
        Quantity const & gradient_moment, Quantity const & G_max);
    
    /**
     * @brief Shortest possible time interval given 3D gradient moment and
     * maximum gradient amplitude
     */
    static TimeInterval shortest(
        Vector3Q const & gradient_moment, Quantity const & G_max);

    /**
     * @brief Constructor, gradient may be specified as amplitude (in T/m), 
     * area (in T/m*s) or dephasing (in rad/m).
     */
    TimeInterval(
        Quantity const & duration=0*units::s,
        Quantity const & gradient=0*units::T/units::m);
    
    /**
     * @brief Constructor, gradient may be specified as amplitude (in T/m), 
     * area (in T/m*s) or dephasing (in rad/m).
     */
    TimeInterval(Quantity const & duration, Vector3Q const & gradient);

    /// @brief Return the duration.
    Quantity const & get_duration() const;
    
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

    /// @brief Return the gradient dephasing.
    Vector3Q get_gradient_moment() const;
    
    /// @brief Set the gradient dephasing.
    void set_gradient_moment(Quantity const & q);
    
    /// @brief Set the gradient dephasing.
    void set_gradient_moment(Vector3Q const & a);
    
    /// @brief Return the gradient amplitude
    Vector3Q const & get_gradient_amplitude() const;
    
    /// @brief Set the gradient amplitude
    void set_gradient_amplitude(Quantity const & q);
    
    /// @brief Set the gradient amplitude
    void set_gradient_amplitude(Vector3Q const & a);
    
    /// @brief Return the gradient area
    Vector3Q get_gradient_area() const;
    
    /// @brief Set the gradient area
    void set_gradient_area(Quantity const & q);
    
    /// @brief Set the gradient area
    void set_gradient_area(Vector3Q const & a);
    
    /// @brief Return the gradient dephasing.
    Vector3Q get_gradient_dephasing() const;
    
    /// @brief Set the gradient dephasing
    void set_gradient_dephasing(Quantity const & q);
    
    /// @brief Set the gradient dephasing
    void set_gradient_dephasing(Vector3Q const & a);
    
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
