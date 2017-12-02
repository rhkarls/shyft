#include "boostpython_pch.h"

#include "core/hbv_actual_evapotranspiration.h"

namespace expose {
	using namespace shyft::core::hbv_actual_evapotranspiration;
	using namespace boost::python;


	void hbv_actual_evapotranspiration() {
		class_<parameter>("HbvActualEvapotranspirationParameter")
			.def(init<optional<double>>(args("lp"), "a new object with specified parameters"))
			.def_readwrite("lp", &parameter::lp, "typical value 150")
			;
		class_<response>("HbvActualEvapotranspirationResponse")
			.def_readwrite("ae", &response::ae)
			;
		def("ActualEvapotranspirationCalculate_step", calculate_step, args("soil-moisture", "potential_evapotranspiration", "lp", "snow_fraction"),
			doc_intro(" actual_evapotranspiration calculates actual evapotranspiration, returning same unit as input pot.evap")
			doc_intro(" based on supplied parameters")
			doc_parameters()
			doc_parameter("water_level", "float"," unit [mm]")
			doc_parameter("potential_evapotranspiration","float","unit[mm/x]")
			doc_parameter("soil_moisture threshold","float","typically 150[mm]")
			doc_parameter("lp","float","param snow_fraction 0..1")
			doc_returns("calculated actual evapotranspiration","float","unit[mm/x]")
		);
	}
}
