#include "boostpython_pch.h"
#include <boost/python/make_constructor.hpp>

#include "core/utctime_utilities.h"
#include "core/time_axis.h"

namespace expose {
    namespace time_axis {
        using namespace shyft::core;
        using namespace boost::python;
		namespace py = boost::python;
        using namespace std;
		using namespace shyft::time_axis;

        template <class C_,class PyC>
        PyC& e_time_axis_std(PyC& c) {
                return c.def("total_period",&C_::total_period,
                     doc_returns("total_period","UtcPeriod","the period that covers the entire time-axis")
                )
                .def("size",&C_::size,
                     doc_returns("n","int","number of periods in time-axis")
                )
                .def("time",&C_::time,args("i"),
                     doc_parameters()
                     doc_parameter("i","int","the i'th period, 0..n-1")
                     doc_returns("utctime","int","the start(utctime) of the i'th period of the time-axis")
                )
                .def("period",&C_::period,args("i"),
                     doc_parameters()
                     doc_parameter("i","int","the i'th period, 0..n-1")
                     doc_returns("period","UtcPeriod","the i'th period of the time-axis")
                )
                .def("index_of",&C_::index_of,args("t"),
                     doc_parameters()
                     doc_parameter("t","int","utctime in seconds 1970.01.01")
                     doc_returns("index","int","the index of the time-axis period that contains t, npos if outside range")
                )
                .def("open_range_index_of",&C_::open_range_index_of,args("t"),
                     "returns the index that contains t, or is before t"
                     doc_parameters()
                     doc_parameter("t","int","utctime in seconds 1970.01.01")
                     doc_returns("index","int","the index the time-axis period that contains t, -npos if before first period n-1, if t is after last period")
                )
                .def(self == self)
                .def(self != self)
                    ;

        }
		struct fixed_dt_ext {
			static fixed_dt *create_default() {
				return new fixed_dt{};
			}
			static fixed_dt *create_from_ints(int64_t t0, int64_t dt, int64_t n) {
				return new fixed_dt{ utctime{ seconds(t0) }, seconds(dt), size_t(n) };
			}
			static fixed_dt *create_from_t_ints(utctime t0, int64_t dt, int64_t n) {
				return new fixed_dt{t0, seconds(dt), size_t(n) };
			}
			static fixed_dt *create_from(utctime t0, utctimespan dt,int64_t n) {
				return new fixed_dt{ t0,dt,size_t(n) };
			}
		};

		static void e_fixed_dt() {
            auto f_dt=class_<fixed_dt>("TimeAxisFixedDeltaT",
                    doc_intro("A time-axis is a set of ordered non-overlapping periods,")
                    doc_intro("and this class implements a fixed delta-t time-axis by")
                    doc_intro("specifying the minimal t-start, delta-t and number of consecutive periods.")
                    doc_intro("This class is a wrapper of the Shyft core-library high performance time-axis")
                    doc_see_also("TimeAxisCalendarDeltaT,TimeAxisByPoints,TimeAxis"),no_init
                )
				.def("__init__", make_constructor(&fixed_dt_ext::create_default),
					doc_intro("Construct an empty class")
				)
				.def("__init__", make_constructor(&fixed_dt_ext::create_from_ints,
					default_call_policies(),
					(py::arg("t0"), py::arg("dt"), py::arg("n"))
					),
					doc_intro("creates a time-axis with n intervals, fixed delta_t, starting at start")
					doc_parameters()
					doc_parameter("start", "int", "utc-time 1970 utc based")
					doc_parameter("delta_t", "int", "number of seconds delta-t, length of periods in the time-axis")
					doc_parameter("n", "int", "number of periods in the time-axis")
				)
				.def("__init__", make_constructor(&fixed_dt_ext::create_from_t_ints,
					default_call_policies(),
					(py::arg("t0"), py::arg("dt"), py::arg("n"))
					),
					doc_intro("creates a time-axis with n intervals, fixed delta_t, starting at start")
					doc_parameters()
					doc_parameter("start", "UtcTime", "utc-time 1970 utc based")
					doc_parameter("delta_t", "int", "number of seconds delta-t, length of periods in the time-axis")
					doc_parameter("n", "int", "number of periods in the time-axis")
				)
				.def("__init__", make_constructor(&fixed_dt_ext::create_from,
					default_call_policies(),
					(py::arg("t0"), py::arg("dt"), py::arg("n"))
					),
					doc_intro("creates a time-axis with n intervals, fixed delta_t, starting at start")
					doc_parameters()
					doc_parameter("start", "UtcTime", "utc-time 1970 utc based")
					doc_parameter("delta_t", "TimeSpan", "number of seconds delta-t, length of periods in the time-axis")
					doc_parameter("n", "int", "number of periods in the time-axis")
				)


                .def_readonly("n",&fixed_dt::n,"int,number of periods")
                .def_readonly("start",&fixed_dt::t,"start of the time-axis, in seconds since 1970.01.01 UTC")
                .def_readonly("delta_t",&fixed_dt::dt,"int,time-span of each interval in seconds")
                .def("full_range",&fixed_dt::full_range,"returns a timeaxis that covers [-oo..+oo> ").staticmethod("full_range")
                .def("null_range",&fixed_dt::null_range,"returns a null timeaxis").staticmethod("null_range");
                ;
                e_time_axis_std<fixed_dt>(f_dt);
        }

		struct calendar_dt_ext {
			static calendar_dt *create_default() {
				return new calendar_dt{};
			}
			static calendar_dt *create_from_ints(shared_ptr<calendar>&cal,int64_t t0, int64_t dt, int64_t n) {
				return new calendar_dt{ cal,utctime{ seconds(t0) }, seconds(dt), size_t(n) };
			}
			static calendar_dt *create_from(shared_ptr<calendar>&cal, utctime t0, utctimespan dt, int64_t n) {
				return new calendar_dt{cal, t0,dt,size_t(n) };
			}
		};


		static void e_calendar_dt() {
            auto c_dt=class_<calendar_dt>("TimeAxisCalendarDeltaT",
                    doc_intro("A time-axis is a set of ordered non-overlapping periods,")
                    doc_intro("and this class implements a calendar-unit fixed delta-t time-axis by")
                    doc_intro("specifying the minimal calendar,t-start, delta-t and number of consecutive periods.")
                    doc_intro("This class is particularly useful if you need to work with calendar and daylight-saving time,")
                    doc_intro("or calendar-specific periods like day,week,month,quarters or years.")
                    doc_see_also("TimeAxisFixedDeltaT,TimeAxisByPoints,TimeAxis"),no_init
                )
				.def("__init__", make_constructor(&calendar_dt_ext::create_default),
					doc_intro("Construct an empty class")
				)
				.def("__init__", make_constructor(&calendar_dt_ext::create_from_ints,
					default_call_policies(),
					(py::arg("calendar"), py::arg("start"), py::arg("delta_t"), py::arg("n"))
					),
					doc_intro("creates a time-axis with n intervals, fixed delta_t, starting at start")
					doc_parameters()
					doc_parameter("calendar", "Calendar", "specifies the calendar to be used, keeps the time-zone and dst-arithmetic rules")
					doc_parameter("start", "int", "utc-time 1970 utc based")
					doc_parameter("delta_t", "int", "number of seconds delta-t, length of periods in the time-axis")
					doc_parameter("n", "int", "number of periods in the time-axis")
				)
				.def("__init__", make_constructor(&calendar_dt_ext::create_from,
					default_call_policies(),
					(py::arg("calendar"), py::arg("start"), py::arg("delta_t"), py::arg("n"))
					),
					doc_intro("creates a time-axis with n intervals, fixed delta_t, starting at start")
					doc_parameters()
					doc_parameter("calendar", "Calendar", "specifies the calendar to be used, keeps the time-zone and dst-arithmetic rules")
					doc_parameter("start", "UtcTime", "utc-time 1970 utc based")
					doc_parameter("delta_t", "TimeSpan", "number of seconds delta-t, length of periods in the time-axis")
					doc_parameter("n", "int", "number of periods in the time-axis")
				)

                .def_readonly("n",&calendar_dt::n,"int,number of periods")
                .def_readonly("start",&calendar_dt::t,"start of the time-axis, in seconds since 1970.01.01 UTC")
                .def_readonly("delta_t",&calendar_dt::dt,"int,timespan of each interval,use Calendar.DAY|.WEEK|.MONTH|.QUARTER|.YEAR, or seconds")
                .add_property("calendar",&calendar_dt::get_calendar,"Calendar, the calendar of the time-axis")
                ;
            e_time_axis_std<calendar_dt>(c_dt);
        }

		template<class ta>
		struct point_dt_ext {
			static ta *create_default() {	return new ta{};}
			static ta *create_from_t_end(const vector<utctime>&t, utctime t_end) {return new ta{ t,t_end };};
			static ta *create_from_t_end_int(const vector<utctime>&t, int64_t t_end) {
				return new ta{ t, utctime{utctimespan{seconds(t_end)}} };
			};
			static ta *create_from_tp(const vector<utctime>&t) { return new ta{ t }; };
		};


        static void e_point_dt() {
            auto p_dt=class_<point_dt>("TimeAxisByPoints",
                    doc_intro("A time-axis is a set of ordered non-overlapping periods,")
                    doc_intro("and this class implements this by a set of ")
                    doc_intro("ordered unique time-points. This is the most flexible time-axis representation,")
                    doc_intro("that allows every period in the time-axis to have different length.")
                    doc_intro("It comes at the cost of space&performance in certain cases, so ")
                    doc_intro("avoid use in scenarios where high-performance is important.")
                    doc_see_also("TimeAxisCalendarDeltaT,TimeAxisFixedDeltaT,TimeAxis"),no_init
                )
				.def("__init__", make_constructor(&point_dt_ext<point_dt>::create_default),
					doc_intro("Construct an empty class")
				)
				.def("__init__", make_constructor(&point_dt_ext<point_dt>::create_from_t_end,
						default_call_policies(),
						(py::arg("time_points"), py::arg("t_end"))
					),
					doc_intro("creates a time-axis by specifying the time_points and t-end of the last interval")
					doc_parameters()
					doc_parameter("time_points", "UtcTimeVector", "ordered set of unique utc-time points, the start of each consecutive period")
					doc_parameter("t_end", "UtcTime", "the end of the last period in time-axis, utc-time 1970 utc based, must be > time_points[-1]")

				)

				.def("__init__", make_constructor(&point_dt_ext<point_dt>::create_from_t_end_int,
						default_call_policies(),
						(py::arg("time_points"), py::arg("t_end"))
					),
					doc_intro("creates a time-axis by specifying the time_points and t-end of the last interval")
					doc_parameters()
					doc_parameter("time_points", "UtcTimeVector", "ordered set of unique utc-time points, the start of each consecutive period")
					doc_parameter("t_end", "int", "the end of the last period in time-axis, utc-time 1970 utc based, must be > time_points[-1]")
				)
				.def("__init__", make_constructor(&point_dt_ext<point_dt>::create_from_tp,
						default_call_policies(),
						(py::arg("time_points"))
					),
					doc_intro("create a time-axis supplying n+1 points to define n intervals")
					doc_parameters()
					doc_parameter("time_points", "UtcTimeVector", "ordered set of unique utc-time points, 0..n-2:the start of each consecutive period,n-1: end of last period")
				)
                .def_readonly("t",&point_dt::t,"UtcTimeVector,time_points except end of last period, see t_end")
                .def_readonly("t_end",&point_dt::t_end,"utctime: end of time-axis")
                ;
                e_time_axis_std<point_dt>(p_dt);
        }

        static std::vector<utctime> time_axis_extract_time_points(shyft::time_axis::generic_dt const&ta) {
            std::vector<utctime> r;r.reserve(ta.size() + 1);
            for (size_t i = 0;i < ta.size();++i) {
                r.emplace_back(ta.time(i));
            }
            if (ta.size())
                r.emplace_back(ta.total_period().end);
            return r;
        }
		struct generic_dt_ext:point_dt_ext<generic_dt> {
			static generic_dt *create_default() {return new generic_dt{};}
			static generic_dt *create_f_from_ints( int64_t t0, int64_t dt, int64_t n) {
				return new generic_dt{ utctime{ seconds(t0) }, seconds(dt), size_t(n) };
			}
			static generic_dt *create_c_from_ints(shared_ptr<calendar>&cal, int64_t t0, int64_t dt, int64_t n) {
				return new generic_dt{ cal,utctime{ seconds(t0) }, seconds(dt), size_t(n) };
			}

			static generic_dt *create_f_from(utctime t0, utctimespan dt, int64_t n) {
				return new generic_dt{ t0,dt,size_t(n) };
			}
			static generic_dt *create_c_from(shared_ptr<calendar>&cal, utctime t0, utctimespan dt, int64_t n) {
				return new generic_dt{ cal, t0,dt,size_t(n) };
			}
			template<class T>
			static generic_dt *create_from_t(const T& t) {
				return new generic_dt{ t };
			}
		};

        static void e_generic_dt() {
            enum_<generic_dt::generic_type>("TimeAxisType")
            .value("FIXED",generic_dt::generic_type::FIXED)
            .value("CALENDAR",generic_dt::generic_type::CALENDAR)
            .value("POINT",generic_dt::generic_type::POINT)
            .export_values()
            ;

            auto g_dt=class_<generic_dt>("TimeAxis",
                    doc_intro("A time-axis is a set of ordered non-overlapping periods,")
                    doc_intro("and TimeAxis provides the most generic implementation of this.")
                    doc_intro("The internal representation is selected based on provided parameters")
                    doc_intro("to the constructor.")
                    doc_intro("The internal representation is one of TimeAxis FixedDeltaT CalendarDelataT or ByPoints.")
                    doc_intro("The internal representation type and corresponding realizations are available as properties.")
                    doc_see_also("TimeAxisCalendarDeltaT,TimeAxisFixedDeltaT,TimeAxisByPoints"),no_init
                )
				.def("__init__", make_constructor(&generic_dt_ext::create_default),
					doc_intro("Construct an empty class")
				)

				.def("__init__", make_constructor(&generic_dt_ext::create_f_from_ints,
					default_call_policies(),
					(py::arg("start"), py::arg("delta_t"), py::arg("n"))
					),
					doc_intro("creates a time-axis with n intervals, fixed delta_t, starting at start")
					doc_parameters()
					doc_parameter("start", "int", "utc-time 1970 utc based")
					doc_parameter("delta_t", "int", "number of seconds delta-t, length of periods in the time-axis")
					doc_parameter("n", "int", "number of periods in the time-axis")
				)
				.def("__init__", make_constructor(&generic_dt_ext::create_f_from,
					default_call_policies(),
					(py::arg("start"), py::arg("delta_t"), py::arg("n"))
					),
					doc_intro("creates a time-axis with n intervals, fixed delta_t, starting at start")
					doc_parameters()
					doc_parameter("start", "UtcTime", "utc-time 1970 utc based")
					doc_parameter("delta_t", "TimeSpan", "number of seconds delta-t, length of periods in the time-axis")
					doc_parameter("n", "int", "number of periods in the time-axis")
				)

                
				.def("__init__", make_constructor(&generic_dt_ext::create_c_from_ints,
					default_call_policies(),
					(py::arg("calendar"), py::arg("start"), py::arg("delta_t"), py::arg("n"))
					),
					doc_intro("creates a time-axis with n intervals, fixed delta_t, starting at start")
					doc_parameters()
					doc_parameter("calendar", "Calendar", "specifies the calendar to be used, keeps the time-zone and dst-arithmetic rules")
					doc_parameter("start", "int", "utc-time 1970 utc based")
					doc_parameter("delta_t", "int", "number of seconds delta-t, length of periods in the time-axis")
					doc_parameter("n", "int", "number of periods in the time-axis")
				)
				.def("__init__", make_constructor(&generic_dt_ext::create_c_from,
					default_call_policies(),
					(py::arg("calendar"), py::arg("start"), py::arg("delta_t"), py::arg("n"))
					),
					doc_intro("creates a time-axis with n intervals, fixed delta_t, starting at start")
					doc_parameters()
					doc_parameter("calendar", "Calendar", "specifies the calendar to be used, keeps the time-zone and dst-arithmetic rules")
					doc_parameter("start", "UtcTime", "utc-time 1970 utc based")
					doc_parameter("delta_t", "TimeSpan", "number of seconds delta-t, length of periods in the time-axis")
					doc_parameter("n", "int", "number of periods in the time-axis")
				)


				.def("__init__", make_constructor(&generic_dt_ext::create_from_t_end,
					default_call_policies(),
					(py::arg("time_points"), py::arg("t_end"))
					),
					doc_intro("creates a time-axis by specifying the time_points and t-end of the last interval")
					doc_parameters()
					doc_parameter("time_points", "UtcTimeVector", "ordered set of unique utc-time points, the start of each consecutive period")
					doc_parameter("t_end", "UtcTime", "the end of the last period in time-axis, utc-time 1970 utc based, must be > time_points[-1]")

				)

				.def("__init__", make_constructor(&generic_dt_ext::create_from_t_end_int,
					default_call_policies(),
					(py::arg("time_points"), py::arg("t_end"))
					),
					doc_intro("creates a time-axis by specifying the time_points and t-end of the last interval")
					doc_parameters()
					doc_parameter("time_points", "UtcTimeVector", "ordered set of unique utc-time points, the start of each consecutive period")
					doc_parameter("t_end", "int", "the end of the last period in time-axis, utc-time 1970 utc based, must be > time_points[-1]")
				)
				.def("__init__", make_constructor(&generic_dt_ext::create_from_tp,
					default_call_policies(),
					(py::arg("time_points"))
					),
					doc_intro("create a time-axis supplying n+1 points to define n intervals")
					doc_parameters()
					doc_parameter("time_points", "UtcTimeVector", "ordered set of unique utc-time points, 0..n-2:the start of each consecutive period,n-1: end of last period")
				)




                .def("__init__", make_constructor(&generic_dt_ext::create_from_t<calendar_dt>,
					default_call_policies(),
                    (py::arg("calendar_dt"))
					),
                    doc_intro("create a time-axis from a calendar time-axis")
                    doc_parameters()
                    doc_parameter("calendar_dt","TimeAxisCalendarDeltaT","existing calendar time-axis")
                )
                .def("__init__", make_constructor(&generic_dt_ext::create_from_t<fixed_dt>,
					default_call_policies(),
                    (py::arg("fixed_dt"))
					),
                    doc_intro("create a time-axis from a a fixed delta-t time-axis")
                    doc_parameters()
                    doc_parameter("fixed_dt","TimeAxisFixedDeltaT","existing fixed delta-t time-axis")
                )
                .def("__init__", make_constructor(&generic_dt_ext::create_from_t<point_dt>,
					default_call_policies(),
                    (py::arg("point_dt"))
					),
                    doc_intro("create a time-axis from a a by points  time-axis")
                    doc_parameters()
                    doc_parameter("point_dt","TimeAxisByPoints","existing by points time-axis")
                )


                .def_readonly("timeaxis_type",&generic_dt::gt,"describes what time-axis representation type this is,e.g (fixed|calendar|point)_dt ")
                .def_readonly("fixed_dt",&generic_dt::f,"The fixed dt representation (if active)")
                .def_readonly("calendar_dt",&generic_dt::c,"The calendar dt representation(if active)")
                .def_readonly("point_dt",&generic_dt::p,"The point_dt representation(if active)")
                .def("total_period", &generic_dt::total_period,
                    doc_returns("total_period", "UtcPeriod", "the period that covers the entire time-axis")
                )
                .def("size", &generic_dt::size,
                    doc_returns("n", "int", "number of periods in time-axis")
                )
                .def("time", &generic_dt::time, args("i"),
                    doc_parameters()
                    doc_parameter("i", "int", "the i'th period, 0..n-1")
                    doc_returns("utctime", "int", "the start(utctime) of the i'th period of the time-axis")
                )
                .def("period", &generic_dt::period, args("i"),
                    doc_parameters()
                    doc_parameter("i", "int", "the i'th period, 0..n-1")
                    doc_returns("period", "UtcPeriod", "the i'th period of the time-axis")
                )
                .def("index_of", &generic_dt::index_of, (py::arg("t"),py::arg("ix_hint")=string::npos),
                    doc_parameters()
                    doc_parameter("t", "int", "utctime in seconds 1970.01.01")
                    doc_parameter("ix_hint","int","index-hint to make search in point-time-axis faster")
                    doc_returns("index", "int", "the index of the time-axis period that contains t, npos if outside range")
                )
                .def("open_range_index_of", &generic_dt::open_range_index_of, (py::arg("t"), py::arg("ix_hint") = string::npos),
                    "returns the index that contains t, or is before t"
                    doc_parameters()
                    doc_parameter("t", "int", "utctime in seconds 1970.01.01")
                    doc_parameter("ix_hint", "int", "index-hint to make search in point-time-axis faster")
                    doc_returns("index", "int", "the index the time-axis period that contains t, npos if before first period n-1, if t is after last period")
                )
                .def(self == self)
                .def(self != self)
                ;//.def("full_range",&point_dt::full_range,"returns a timeaxis that covers [-oo..+oo> ").staticmethod("full_range")
                //.def("null_range",&point_dt::null_range,"returns a null timeaxis").staticmethod("null_range");
            //e_time_axis_std<generic_dt>(g_dt);
            def("time_axis_extract_time_points", time_axis_extract_time_points, args("time_axis"),
                doc_intro("Extract all time_axis.period(i).start plus time_axis.total_period().end into a UtcTimeVector")
                doc_parameters()
                doc_parameter("time_axis","TimeAxis","time-axis to extract all time-points from")
                doc_returns("time_points","UtcTimeVector","all time_axis.period(i).start plus time_axis.total_period().end")
            );
        }

    }
    void api_time_axis() {
        time_axis::e_fixed_dt();
        time_axis::e_point_dt();
        time_axis::e_calendar_dt();
        time_axis::e_generic_dt();
    }

}
