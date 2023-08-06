#include <sycomore/isochromat/Model.h>
#include <sycomore/units.h>

#include <xtensor/xview.hpp>

#include <xtensor/xio.hpp>

int main()
{
    using namespace sycomore::units;
    
    sycomore::Vector3R M0{0,0,1};
    // 1D model at a single position. This can take a more complex form, e.g. for
    // a 3D model at two positions: {{0*mm, 0*mm, 0*mm}, {1*mm, 1*mm, 1*mm}} 
    sycomore::TensorQ<2>positions{{0*mm}};
    auto T1 = 1000*ms, T2 = 100*ms, flip_angle = 60*deg, step = 10*ms;
    
    sycomore::isochromat::Model model(T1, T2, M0, positions);
    
    auto pulse = model.build_pulse(flip_angle);
    auto idle = model.build_time_interval(10*ms);
    
    std::vector<std::pair<sycomore::Quantity, sycomore::Vector3R>> record{
        {0*s, xt::view(model.magnetization(), 0)}};
    
    for(std::size_t i=0; i!=10; ++i)
    {
        model.apply(idle);
        record.emplace_back(
            record.back().first+step, xt::view(model.magnetization(), 0));
    }
    
    model.apply(pulse);
    record.emplace_back(
        record.back().first+step, xt::view(model.magnetization(), 0));
    
    for(std::size_t i=0; i!=100; ++i)
    {
        model.apply(idle);
        record.emplace_back(
            record.back().first+step, xt::view(model.magnetization(), 0));
    }
    
    return 0;
}
