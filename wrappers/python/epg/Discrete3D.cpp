#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "sycomore/epg/Discrete3D.h"
#include "sycomore/magnetization.h"
#include "sycomore/Species.h"

void wrap_epg_Discrete3D(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::epg;

    class_<Discrete3D, Base>(
            m, "Discrete3D",
            "Discrete EPG in which the gradients may be specified in three "
            "dimensions."
        )
        .def(
            init<Species, Magnetization, Quantity>(),
            arg("species"), arg("initial_magnetization")=Magnetization{0,0,1},
            arg("bin_width")=1*units::rad/units::m)
        .def(
            init<
                Species const &, Species const &,
                Magnetization const &, Magnetization const &,
                Quantity const &, Quantity const &, Quantity const &>(),
            arg("species_a"), arg("species_b"), arg("M0_a"), arg("M0_b"),
            arg("k_a"), arg("delta_b")=0*units::Hz,
            arg("bin_width")=1*units::rad/units::m)
        .def(
            init<
                Species const &, Quantity const &,
                Magnetization const &, Magnetization const &,
                Quantity const &, Quantity const &>(),
            arg("species_a"), arg("R1_b_or_T1_b"), arg("M0_a"), arg("M0_b"),
            arg("k_a"), arg("bin_width")=1*units::rad/units::m)
        .def_property_readonly(
            "orders", 
            [](Discrete3D const & model){
                auto const orders_cpp = model.orders();
                list orders_py(model.size());
                for(std::size_t i=0; i<model.size(); ++i)
                {
                    orders_py[i] = make_tuple(
                        orders_cpp[3*i+0], orders_cpp[3*i+1], orders_cpp[3*i+2]);
                }
                return orders_py;
            },
            "Orders of the model.")
        .def_property_readonly("bin_width", &Discrete3D::bin_width)
        .def(
            "state",
            static_cast<
                    std::vector<Complex> (Discrete3D::*)(std::size_t) const
                >(&Discrete3D::state),
            arg("bin"),
            "Magnetization at a given state, expressed by its *index*")
        .def(
            "state", 
            [](Discrete3D const & model, sequence order_py) {
                Discrete3D::Order order_cpp(order_py.size());
                for(std::size_t i=0; i<order_py.size(); ++i)
                {
                    order_cpp[i] = order_py[i].cast<Quantity>();
                }
                return model.state(order_cpp);
            },
            arg("order"), "Access a given state of the model")
        .def_property_readonly("elapsed", &Discrete3D::elapsed)
        .def(
            "apply_time_interval",
            [](
                Discrete3D & model, Quantity const & duration,
                sequence const & gradient)
            {
                Array<Quantity> array(gradient.size());
                for(std::size_t i=0; i<array.size(); ++i)
                {
                    array[i] = gradient[i].cast<Quantity>();
                }
                model.apply_time_interval(duration, array);
            },
            arg("duration"), arg("gradient")=Array<Quantity>{
                0*units::T/units::m, 0*units::T/units::m, 0*units::T/units::m},
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects. States with a population lower than "
            "*threshold* will be removed.")
        .def(
            "apply_time_interval", 
            static_cast<void(Discrete3D::*)(TimeInterval const &)>(
                &Discrete3D::apply_time_interval),
            arg("interval"),
            "Apply a time interval, i.e. relaxation, diffusion, gradient, and "
            "off-resonance effects. States with a population lower than "
            "*threshold* will be removed.")
        .def(
            "shift", 
            [](Discrete3D & model, Quantity duration, sequence gradient_py) {
                Discrete3D::Order gradient_cpp(gradient_py.size());
                for(std::size_t i=0; i<gradient_py.size(); ++i)
                {
                    gradient_cpp[i] = gradient_py[i].cast<Quantity>();
                }
                return model.shift(duration, gradient_cpp);
            },
            arg("duration"), arg("gradient"),
            "Apply a gradient; in discrete EPG, this shifts all orders by "
            "specified value.")
        .def(
            "diffusion", 
            [](Discrete3D & model, Quantity duration, sequence gradient_py) {
                Discrete3D::Order gradient_cpp(gradient_py.size());
                for(std::size_t i=0; i<gradient_py.size(); ++i)
                {
                    gradient_cpp[i] = gradient_py[i].cast<Quantity>();
                }
                return model.diffusion(duration, gradient_cpp);
            },
            arg("duration"), arg("gradient"),
            "Simulate diffusion during given duration with given gradient "
            "amplitude.")
    ;
}
