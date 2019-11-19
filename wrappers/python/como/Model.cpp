#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "sycomore/como/Model.h"

void wrap_como_Model(pybind11::module & m)
{
    using namespace pybind11;
    using namespace sycomore;
    using namespace sycomore::como;
    using namespace sycomore::units;

    class_<Model>(m, "Model")
        .def(init(
            [](
                Species const & species, Magnetization const & magnetization,
                sequence time_intervals_py)
            {
                std::vector<std::pair<std::string, TimeInterval>>
                    time_intervals_cpp;
                for(std::size_t i=0; i<time_intervals_py.size(); ++i)
                {
                    auto item = time_intervals_py[i].cast<sequence>();
                    time_intervals_cpp.emplace_back(
                        item[0].cast<std::string>(),
                        item[1].cast<TimeInterval>());
                }

                return new Model(species, magnetization, time_intervals_cpp);
            }))
        .def_property(
            "epsilon", &Model::get_epsilon, &Model::set_epsilon,
            "Threshold magnetization for clean-up (default to 0).")
        .def(
            "apply_pulse",
            static_cast<void (Model::*)(Pulse const &)>(&Model::apply_pulse),
            "Apply an RF pulse to the model.")
        .def(
            "apply_pulse",
            static_cast<void (Model::*) (HardPulseApproximation const &)>(
                &Model::apply_pulse),
            "Apply a hard pulse approximation to the model.")
        .def(
            "apply_time_interval", &Model::apply_time_interval,
            "Apply a time interval to the model.")
        .def(
            "magnetization", &Model::magnetization,
            return_value_policy::reference_internal,
            "Return the complex magnetizations.")
        .def(
            "isochromat", &Model::isochromat,
            arg("configurations")=set(), arg("position")=Point(),
            arg("relative_frequency")=Quantity(0, AngularFrequency),
            "Return the isochromat for the given configurations."
            "\n"
            "If no configuration is specified, all configurations are used.");
}
