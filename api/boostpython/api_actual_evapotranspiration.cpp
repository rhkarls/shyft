#include "boostpython_pch.h"

#include "core/actual_evapotranspiration.h"

namespace expose {
    using namespace shyft::core::actual_evapotranspiration;
    namespace py=boost::python;


    void actual_evapotranspiration() {
        py::class_<parameter>("ActualEvapotranspirationParameter")
            .def(py::init<py::optional<double>>(py::args("ae_scale_factor"),"a new object with specified parameters"))
            .def_readwrite("ae_scale_factor",&parameter::ae_scale_factor,"typical value 1.5")
            ;
		py::class_<response>("ActualEvapotranspirationResponse")
            .def_readwrite("ae",&response::ae)
            ;
		py::def("ActualEvapotranspirationCalculate_step",calculate_step, (py::arg("water_level"), py::arg("potential_evapotranspiration"), py::arg("scale_factor"), py::arg("snow_fraction")),
             doc_intro("actual_evapotranspiration calculates actual evapotranspiration, returning same unit as input pot.evap")
             doc_intro("based on supplied parameters")
			 doc_parameters()
             doc_parameter("water_level","float","unit[mm]")
			 doc_parameter("potential_evapotranspiration","float","unit[mm/x]")
			 doc_parameter("scale_factor","float","typically 1.5[mm]")
			 doc_parameter("snow_fraction","float"," range  0..1")
             doc_returns( "calculated actual evapotranspiration","float","units [mm/x]")
            );
    }
}
