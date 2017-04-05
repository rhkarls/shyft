#pragma once
#include "utctime_utilities.h"	
namespace shyft {
	namespace core {
		namespace hbv_infiltration {
			using namespace std;


			struct parameter {
				parameter(double Os = 0.434)
					:Os(Os) {
					if (Os < .0)
						throw runtime_error("Os should be > 0.0");
				}
				double Os = 0.434; // fraction
			};

			struct state {
				state(double O0 = 0.0) :O0(O0) {}
				double O0 = 0.0; // mm
				bool operator==(const state&x) const {
					const double eps = 1e-6;
					return fabs(O0 - x.O0)<eps;
				}
				x_serialize_decl();
			};

			struct response {
				double Freal;
				double Runoff;
			};


			/** \brief Hbv soil
			*
			*  reference somewhere..
			*
			* \tparam P Parameter type, implementing the interface:
			* \tparam S State type, implementing the interface:
			* \tparam R Respone type, implementing the interface:
			*/
			template<class P>
			struct calculator {
				P param;
				calculator(const P& p) :param(p) {}
				template <class R, class S>
				void step(S& s, R& r, shyft::core::utctime t0, shyft::core::utctime t1, double snow_out) {
					r.Freal = param.Os*snow_out;
					r.Runoff = 1 - r.Freal;
					s.O0 = s.O0 + r.Freal;
				}
			};
		}
	} // core
} // shyft
  //-- serialization support shyft
x_serialize_export_key(shyft::core::hbv_infiltration::state);
