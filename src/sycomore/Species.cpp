#include "Species.h"

#include "sycomore/sycomore.h"
#include "sycomore/units.h"

namespace sycomore
{

Species
::Species(Quantity const & R1, Quantity const & R2)
: Species(R1, R2, 0*units::m*units::m/units::s)
{
    // Nothing else. This constructor only exists to disambiguate between
    // the isotropic diffusio version and the anisotropic one.
}


Species
::Species(
    Quantity const & R1, Quantity const & R2, Quantity const & D,
    Quantity const & R2_prime, Quantity const & delta_omega, Real w)
: w(w)
{
    this->set_R1(R1);
    this->set_R2(R2);
    this->set_D(D);
    this->set_R2_prime(R2_prime);
    this->set_delta_omega(delta_omega);
}

Species
::Species(
    Quantity const & R1, Quantity const & R2, Array<Quantity> const & D,
    Quantity const & R2_prime, Quantity const & delta_omega, Real w)
: w(w)
{
    this->set_R1(R1);
    this->set_R2(R2);
    this->set_D(D);
    this->set_R2_prime(R2_prime);
    this->set_delta_omega(delta_omega);
}

Quantity const &
Species
::get_R1() const
{
    return this->_R1;
}

void
Species
::set_R1(Quantity const & q)
{
    if(q.dimensions == Frequency)
    {
        this->_R1 = q;
        this->_T1 = 1/q;
    }
    else if(q.dimensions == Time)
    {
        this->_R1 = 1/q;
        this->_T1 = q;
    }
    else
    {
        std::ostringstream message;
        message << "R1 must be duration or frequency, not " << q.dimensions;
        throw std::runtime_error(message.str());
    }
}

Quantity const &
Species
::get_T1() const
{
    return this->_T1;
}

Quantity const &
Species
::get_R2() const
{
    return this->_R2;
}

void
Species
::set_R2(Quantity const & q)
{
    if(q.dimensions == Frequency)
    {
        this->_R2 = q;
        this->_T2 = 1/q;
    }
    else if(q.dimensions == Time)
    {
        this->_R2 = 1/q;
        this->_T2 = q;
    }
    else
    {
        std::ostringstream message;
        message << "R2 must be duration or frequency, not " << q.dimensions;
        throw std::runtime_error(message.str());
    }
}

Quantity const &
Species
::get_T2() const
{
    return this->_T2;
}

Quantity const &
Species
::get_R2_prime() const
{
    return this->_R2_prime;
}

Array<Quantity> const &
Species
::get_D() const
{
    return this->_D;
}

void
Species
::set_D(Quantity const & q)
{
    using namespace units;

    if(q.dimensions == Diffusion)
    {
        this->_D = {
            q,       0*m*m/s, 0*m*m/s,
            0*m*m/s, q,       0*m*m/s,
            0*m*m/s, 0*m*m/s, q};
    }
    else
    {
        std::ostringstream message;
        message << "D must be a diffusion coefficient, not " << q.dimensions;
        throw std::runtime_error(message.str());
    }
}

void
Species
::set_D(Array<Quantity> const & q)
{
    if(q.size() != 9)
    {
        std::ostringstream message;
        message<< "Diffusion tensor must have 9 elements, not " << q.size();
        throw std::runtime_error(message.str());
    }
    this->_D = Array<Quantity>(q.size());
    for(std::size_t i=0; i<q.size(); ++i)
    {
        if(q[i].dimensions == Diffusion)
        {
            this->_D[i] = q[i];
        }
        else
        {
            std::ostringstream message;
            message
                << "D[" << i << "] must be a diffusion coefficient, not "
                << q[i].dimensions;
            throw std::runtime_error(message.str());
        }
    }
}

void
Species
::set_R2_prime(Quantity const & q)
{
    if(q.dimensions == Frequency)
    {
        this->_R2_prime = q;
        this->_T2_prime = 1/q;
    }
    else if(q.dimensions == Time)
    {
        this->_R2_prime = 1/q;
        this->_T2_prime = q;
    }
    else
    {
        std::ostringstream message;
        message
            << "R2_prime must be duration or frequency, not " << q.dimensions;
        throw std::runtime_error(message.str());
    }
}

Quantity const &
Species
::get_T2_prime() const
{
    return this->_T2_prime;
}

Quantity const &
Species
::get_delta_omega() const
{
    return this->_delta_omega;
}

void
Species
::set_delta_omega(Quantity const & q)
{
    if(q.dimensions == AngularFrequency)
    {
        this->_delta_omega = q;
    }
    else
    {
        std::ostringstream message;
        message
            << "delta_omega must be an angular frequency, not "
            << q.dimensions;
        throw std::runtime_error(message.str());
    }
}

}
