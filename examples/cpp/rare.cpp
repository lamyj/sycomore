#include <iostream>

#include <sycomore/magnetization.h>
#include <sycomore/Model.h>
#include <sycomore/Pulse.h>
#include <sycomore/Species.h>
#include <sycomore/sycomore.h>
#include <sycomore/units.h>

int main()
{
    using namespace sycomore::units;

    sycomore::Species const gray_matter(1000_ms, 100_ms, 1_um*um/ms);
    auto const TR = 700_ms;
    auto const TE = 5_ms;
    auto const train_length = 32;
    sycomore::TimeInterval const half_echo(TE/2);
    sycomore::TimeInterval const idle(TR-train_length*TE);
    sycomore::Model model(
        gray_matter, {0,0,1}, {{"half_echo", half_echo}, {"idle", idle}});

    std::vector<std::pair<sycomore::Real, sycomore::Real>> signal;

    sycomore::Real time=0;
    model.apply_pulse({90_deg, 0_deg});
    for(int echo=0; echo<train_length; ++echo)
    {
        model.apply_time_interval("half_echo");
        time += half_echo.duration;

        model.apply_pulse({180_deg, 0_deg});

        model.apply_time_interval("half_echo");
        time += half_echo.duration;

        auto const m = model.isochromat();
        signal.emplace_back(time, m.transversal());
    }
    model.apply_time_interval("idle");
    time += idle.duration;


    for(auto it = signal.begin(); it!=signal.end(); ++it)
    {
        std::cout << it->first << " " << it->second << "\n";
    }
}
