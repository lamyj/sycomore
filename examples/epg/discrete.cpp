#include <chrono>
#include <iostream>

#include <sycomore/epg/Discrete.h>
#include <sycomore/Species.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

using Clock = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

int main()
{
    using namespace sycomore::units;
    
    sycomore::Species species(1000*ms, 100*ms, 1700*std::pow(um, 2)/s);

    // Sequence parameters
    auto flip_angle = 20*deg;
    sycomore::Vector2Q TE{2*ms, 8*ms};
    auto TR = 10*ms, slice_thickness = 1*mm, tau_readout = 1*ms, q = 70/cm;

    // Motion to k-space extremity and its associated gradient amplitude
    auto k_max = 0.5 * 2*M_PI/slice_thickness;
    auto G_readout = k_max / sycomore::gamma / (tau_readout/2);

    // Diffusion gradient
    auto tau_diffusion = (TE[1]-tau_readout/2) - (TE[0]+tau_readout/2);
    auto G_diffusion = q/(sycomore::gamma_bar*tau_diffusion);

    std::vector<sycomore::epg::Discrete>models{
        sycomore::epg::Discrete(species), sycomore::epg::Discrete(species)};
    models[1].threshold = 1e-6;

    std::size_t repetitions = 4*species.T1()/TR;

    sycomore::TensorC<2> S_plus({models.size(), repetitions}, 0);
    sycomore::TensorC<2> S_minus({models.size(), repetitions}, 0);
    for(std::size_t index=0; index!=models.size(); ++index)
    {
        auto & model = models[index];
        
        auto const begin = Clock::now();
        
        for(std::size_t repetition=0; repetition!=repetitions; ++repetition)
        {
            // RF-pulse and idle until the first read-out 
            model.apply_pulse(flip_angle);
            model.apply_time_interval(TE[0]-tau_readout);
            
            // First echo
            model.apply_time_interval(tau_readout/2, -G_readout);
            model.apply_time_interval(tau_readout/2, G_readout);
            S_plus(index, repetition) = model.echo();
            model.apply_time_interval(tau_readout/2, G_readout);
            
            // Diffusion gradient between the two echoes
            model.apply_time_interval(tau_diffusion, G_diffusion);
            
            // Second echo
            model.apply_time_interval(tau_readout/2, G_readout);
            S_minus(index, repetition) = model.echo();
            model.apply_time_interval(tau_readout/2, G_readout);
            model.apply_time_interval(tau_readout/2, -G_readout);
            
            // Idle until the end of the TR
            model.apply_time_interval(TR-TE[1]-tau_readout);
            
            // Make sure the sequence timing is correct
            if(repetition == 0)
            {
                assert((model.elapsed()-TR)/TR < 1e-6);
            }
        }
        
        auto const end = Clock::now();
        
        std::cout
            << "Threshold:" << model.threshold << ", "
            << model.size() << " orders, "
            << 1e3*Duration(end-begin).count() <<  " ms\n";
    }
        
    return 0;
}
