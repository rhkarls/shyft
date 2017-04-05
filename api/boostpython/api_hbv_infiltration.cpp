#include "boostpython_pch.h"

#include "core/hbv_infiltration.h"

namespace expose {

	void hbv_infiltration() {
		using namespace shyft::core::hbv_infiltration;
		using namespace boost::python;
		using namespace std;

		class_<parameter>("HbvInfiltrationParameter")
			.def(init<optional<double>>(args("Os"), "create parameter object with specifed values"))
			.def_readwrite("Os", &parameter::Os, " Saturated volumetric soil moisture content, default=2.0")
			;

		class_<state>("HbvInfiltrationState")
			.def(init<optional<double>>(args("O0"), "create a state with specified values"))
			.def_readwrite("O0", &state::O0, "Current volumetric soil moisture content []")
			;

		class_<response>("HbvInfiltrationResponse")
			.def_readwrite("Freal", &response::Freal, "from infiltration-routine in [mm]")
			.def_readwrite("Runoff", &response::Runoff, "from infiltration-routine in [mm]")
			;

		typedef  calculator<parameter> HbvInfiltrationCalculator;
		class_<HbvInfiltrationCalculator>("HbvInfiltrationCalculator",
			"tobe done.. \n"
			"\n"
			"\n", no_init
			)
			.def(init<const parameter&>(args("parameter"), "creates a calculator with given parameter"))
			.def("step", &HbvInfiltrationCalculator::step<response, state>, args("state", "response", "t0", "t1", "snow_out"),
				"steps the model forward from t0 to t1, updating state and response")
			;
#if 0
#endif
	}
}
