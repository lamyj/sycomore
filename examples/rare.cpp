#include <sycomore/epg/Regular.h>
#include <sycomore/Species.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;
    
    sycomore::Species const species(1000*ms, 100*ms);
    auto const TE = 4*ms;
    int const train_length = 40;
    
    sycomore::epg::Regular model(species);
    sycomore::TensorC<1> signal({train_length}, 0);
    
    model.apply_pulse(90*deg);
    for(int echo=0; echo<train_length; ++echo)
    {
        model.apply_time_interval(TE/2);
        model.apply_pulse(180*deg);
        model.apply_time_interval(TE/2);
    
        signal[echo] = model.echo();
    }
    
    return 0;
}
