#include <cmath>

#include <sycomore/epg/Regular.h>
#include <sycomore/Species.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    sycomore::Species species(1000*ms, 100*ms);
    
    // Sequence parameters
    auto flip_angle=30*deg, TE=5*ms, TR=25*ms;
    std::vector<sycomore::Quantity> phase_steps{0*deg, 90*deg, 117*deg, 180*deg};
    auto slice_thickness=1*mm, tau_readout=1*ms;
    std::size_t repetitions = std::lround(5*species.T1()/TR);
    
    // Motion to k-space extremity and its associated gradient amplitude
    auto k_max = 0.5 * 2*M_PI / slice_thickness;
    auto G = k_max / sycomore::gamma / (tau_readout/2);
    
    std::vector<sycomore::epg::Regular> models;
    for(std::size_t i=0; i!=phase_steps.size(); ++i)
    {
        models.push_back(sycomore::epg::Regular(species, {0,0,1}, 100, k_max));
    }
    
    sycomore::TensorC<2> signals({models.size(), repetitions}, 0);
    for(std::size_t r=0; r!=repetitions; ++r)
    {
        for(std::size_t index=0; index!=phase_steps.size(); ++index)
        {
            auto & phase_step = phase_steps[index];
            auto & model = models[index];
            
            auto phase = phase_step * 1/2 * (r+1) * r;
            
            // RF-pulse and idle until the readout
            model.apply_pulse(flip_angle, phase);
            model.apply_time_interval(TE-tau_readout);
            
            // Readout prephasing and first half of the readout
            model.apply_time_interval(-G, tau_readout/2);
            model.apply_time_interval(+G, tau_readout/2);
            
            // Echo at the center of the readout, cancel the phase imparted by
            // the RF-spoiling
            signals(index, r) =
                model.echo() * std::exp(sycomore::Complex(0, phase));
            
            // Second half of the readout, idle until the end of the TR
            model.apply_time_interval(+G, tau_readout/2);
            model.apply_time_interval(TR-TE-tau_readout/2);
        }
    }
    
    return 0;
}
