#include "boostpython_pch.h"
#include <boost/python/implicit.hpp>
#include <boost/python/dict.hpp>
#include <boost/python/raw_function.hpp>
#include "core/utctime_utilities.h"

namespace expose {
    using namespace shyft::core;
    using namespace boost::python;
    using namespace std;
	namespace py = boost::python;
    typedef std::vector<utcperiod> UtcPeriodVector;


  


    static void e_utcperiod() {
        bool (utcperiod::*contains_t)(utctime) const = &utcperiod::contains;
        bool (utcperiod::*contains_p)(const utcperiod&) const = &utcperiod::contains;
        class_<utcperiod>("UtcPeriod","UtcPeriod defines the open utctime range [start..end> \nwhere end is required to be equal or greater than start")
        .def(init<utctime,utctime>(args("start,end"),"Create utcperiod given start and end"))
		.def(init<int64_t, int64_t>(args("start,end"), "Create utcperiod given start and end"))
        .def("valid",&utcperiod::valid,"returns true if start<=end otherwise false")
        .def("contains",contains_t,args("t"),"returns true if utctime t is contained in this utcperiod" )
        .def("contains",contains_p,args("p"),"returns true if utcperiod p is contained in this utcperiod" )
        .def("overlaps",&utcperiod::overlaps,args("p"), "returns true if period p overlaps this utcperiod" )
        .def("__str__",&utcperiod::to_string,"returns the str using time-zone utc to convert to readable time")
        .def(self == self)
        .def(self != self)
        .def("timespan",&utcperiod::timespan,"returns end-start, the timespan of the period")
        .def_readwrite("start",&utcperiod::start,"Defines the start of the period, inclusive")
        .def_readwrite("end",&utcperiod::end,"Defines the end of the period, not inclusive");
        def("intersection",&intersection,args("a,b"),"Returns the intersection of two utcperiods");
    }
    static bool is_npos(size_t n) {
        return n==string::npos;
    }

	template<class T>
	static T x_arg(const py::tuple& args, size_t i) {
		if (py::len(args) + 1 < (int)i)
			throw std::runtime_error("missing arg #" + std::to_string(i) + std::string(" in UtcTime"));
		py::object o = args[i];
		py::extract<T> xtract_arg(o);
		return xtract_arg();
	}
	static void args_check(const py::tuple& args) {
		if (py::len(args) < 2) 	throw std::runtime_error("compare needs two args");
	}

	struct utctimespan_ext {
		static utctimespan* create_default() {
			return new utctimespan{ 0 };
		}
		static utctimespan* create_from_int(int sec) {
			return new utctimespan(seconds(sec));
		}
		static utctimespan* create_from_double(double sec) {
			return new utctimespan{ from_seconds(sec) };
		}
		static utctimespan x_self(const py::tuple& args) {
			if (py::len(args) == 0)
				throw std::runtime_error("self is null inTimeSpan");
			py::object self = args[0];
			py::extract<utctimespan> xtract_self(self);
			return xtract_self();
		}
		static py::object get_seconds(py::tuple args, py::dict kwargs) {
			utctimespan dt = x_self(args);
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			if (dt_s == dt)
				return py::object(int64_t(dt_s.count()));
			return py::object(to_seconds(dt));
		}
		static py::object get_int_seconds(py::tuple args, py::dict kwargs) {
			utctimespan dt = x_self(args);
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			if (dt_s == dt)
				return py::object(int64_t(dt_s.count()));
			return py::object(int64_t(to_seconds(dt)));
		}

		static py::object repr(py::tuple args, py::dict kwargs) {
			utctimespan dt = x_self(args);
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			char s[100];
			if (dt_s == dt)
#ifndef _WIN32
				sprintf(s, "TimeSpan(%ld)",int64_t(dt_s.count()));
#else
				sprintf(s, "TimeSpan(%lld)", int64_t(dt_s.count()));
#endif
			else
				sprintf(s, "TimeSpan(%0.6lf)", to_seconds(dt));
			return py::str(std::string(s));
		}
		static py::object str(py::tuple args, py::dict kwargs) {
			utctimespan dt = x_self(args);
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			char s[100];
			if (dt_s == dt)
#ifndef _WIN32
				sprintf(s, "%lds", dt_s.count());
#else
				sprintf(s, "%llds", dt_s.count());
#endif
			else
				sprintf(s, "%0.6lfs", to_seconds(dt));
			return py::str(std::string(s));
		}

		static utctimespan abs_timespan(utctimespan x) {
			return x>=x.zero() ? x:-x;
		}
		static double _float_(utctimespan x) {
			return to_seconds(x);
		}
		/// best effort convert from any number or string to utctime 
		/// as long as it is well defined
		static utctimespan as_timespan(const py::object& po) {
			py::extract<utctimespan> o(po);
			if (o.check())
				return o();
			py::extract<int64_t> i(po);
			if (i.check())
				return  from_seconds(i());
			py::extract<double> d(po);
			if (d.check())
				return   from_seconds(d());
			throw std::runtime_error("supplied argument not convertible to TimeSpan");
		}
		// rel ops
		static py::object _lt_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctimespan>(args, 0) < as_timespan(args[1]));}
		static py::object _le_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctimespan>(args, 0) <= as_timespan(args[1]));}
		static py::object _gt_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctimespan>(args, 0) > as_timespan(args[1]));}
		static py::object _ge_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctimespan>(args, 0) >= as_timespan(args[1]));}
		static py::object _eq_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctimespan>(args, 0) == as_timespan(args[1]));}
		static py::object _nq_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctimespan>(args, 0) != as_timespan(args[1]));}
		// math
		static py::object _add_(py::tuple args, py::dict kwargs) { args_check(args); return py::object(x_arg<utctimespan>(args, 0) + as_timespan(args[1])); }
		static py::object _sub_(py::tuple args, py::dict kwargs) { args_check(args); return py::object(x_arg<utctimespan>(args, 0) - as_timespan(args[1])); }
		static py::object _rsub_(py::tuple args, py::dict kwargs) { args_check(args); return py::object(-(x_arg<utctimespan>(args, 0) - as_timespan(args[1]))); }
		static py::object _mult_(py::tuple args, py::dict kwargs) {	args_check(args); return py::object(from_seconds( to_seconds(x_arg<utctimespan>(args, 0)) * to_seconds(as_timespan(args[1])) ));}
		static py::object _div_(py::tuple args, py::dict kwargs) { args_check(args); return py::object(from_seconds( to_seconds(x_arg<utctimespan>(args, 0)) / to_seconds(as_timespan(args[1])) ));	}
		static py::object _mod_(py::tuple args, py::dict kwargs) {	args_check(args);return py::object( x_arg<utctimespan>(args, 0) % as_timespan(args[1]) );}
		static py::object _floordiv_(py::tuple args, py::dict kwargs) {	args_check(args);return py::object(from_seconds(truncf(to_seconds(x_arg<utctimespan>(args, 0)) / to_seconds(as_timespan(args[1])))));}
		static py::object _rfloordiv_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(from_seconds(truncf(to_seconds(as_timespan(args[1]))/to_seconds(x_arg<utctimespan>(args, 0)) )));}
		static py::object _rdiv_(py::tuple args, py::dict kwargs) {	args_check(args);return py::object(from_seconds( to_seconds(as_timespan(args[1]))/ to_seconds(x_arg<utctimespan>(args, 0))));}
		// just for backward(consider change code)
		static py::object _sqrt_(py::tuple args, py::dict kwargs) {return py::object(sqrt(to_seconds(as_timespan(args[0])) / to_seconds(x_arg<utctimespan>(args, 0)))); }
	};

	static void e_utctimespan() {
		class_<utctimespan>("TimeSpan", doc_intro("Is just a number, in unit seconds, with some math-properties"),no_init)
			.def("__init__", make_constructor(&utctimespan_ext::create_default,
					default_call_policies()
				),
				doc_intro("construct a 0s time-span")
			)
			.def("__init__", make_constructor(&utctimespan_ext::create_from_int,
				default_call_policies(),
				(py::arg("seconds"))
				),
				doc_intro("construct a from seconds as integer")
				doc_parameters()
				doc_parameter("seconds", "int", "seconds")
			)
			.def("__init__", make_constructor(&utctimespan_ext::create_from_double,
				default_call_policies(),
				(py::arg("seconds"))
				),
				doc_intro("construct a from seconds as decimal number, where fractions is fractions of second")
				doc_intro("- the resulting time-span preserves 1 micro-second digits")
				doc_parameters()
				doc_parameter("seconds", "int", "seconds")
			)
			.add_property("seconds",raw_function(utctimespan_ext::get_seconds,1),doc_intro("returns timespan in seconds"))
			.def("__abs__",&utctimespan_ext::abs_timespan,(py::arg("self")))
			.def("__float__",&utctimespan_ext::_float_,(py::arg("self")))
			.def("__int__",raw_function(utctimespan_ext::get_int_seconds,1),doc_intro("timespan as int seconds"))
			.def("__long__", raw_function(utctimespan_ext::get_int_seconds, 1), doc_intro("timespan as int seconds"))
			.def("__repr__",raw_function(utctimespan_ext::repr,1),doc_intro("repr of TimeSpan"))
			.def("__str__", raw_function(utctimespan_ext::str, 1), doc_intro("str of TimeSpan"))
			// rel operations
			.def("__lt__", raw_function(utctimespan_ext::_lt_, 1))
			.def("__le__", raw_function(utctimespan_ext::_le_, 1))
			.def("__gt__", raw_function(utctimespan_ext::_gt_, 1))
			.def("__ge__", raw_function(utctimespan_ext::_ge_, 1))
			.def("__eq__", raw_function(utctimespan_ext::_eq_, 1))
			.def("__nq__", raw_function(utctimespan_ext::_nq_, 1))

			// math
			.def("__add__",raw_function(utctimespan_ext::_add_, 1))
			.def("__radd__", raw_function(utctimespan_ext::_add_, 1))
			.def("__sub__", raw_function(utctimespan_ext::_sub_, 1))
			.def("__rsub__", raw_function(utctimespan_ext::_rsub_, 1))
			.def("__truediv__",raw_function(utctimespan_ext::_div_, 1))
			.def("__rtruediv__", raw_function(utctimespan_ext::_rdiv_, 1))
			.def("__floordiv__", raw_function(utctimespan_ext::_floordiv_, 1))
			.def("__rfloordiv__", raw_function(utctimespan_ext::_rfloordiv_, 1))
			.def("__mod__", raw_function(utctimespan_ext::_mod_, 1))

			.def("__mul__", raw_function(utctimespan_ext::_mult_, 1))
			.def("__rmul__", raw_function(utctimespan_ext::_mult_, 1))
			.def("sqrt",raw_function(utctimespan_ext::_sqrt_,1))
			//.def(self % self)
			.def(-self)
			;

	}

	struct utctime_ext {
		static utctime* create_default() {
			return new utctime{};
		}
		static utctime* create_from_int(int sec) {
			return new utctime{ seconds(sec) };
		}
		static utctime* create_from_double(double sec) {
			return new utctime{ from_seconds(sec) };
		}
		static utctime* create_from_timespan(utctimespan sec) {
			return new utctime{ sec };
		}
		// consider spirit for faster and accurate parsing !
		static utctime* create_from_string(const std::string& s) {
			return new utctime{ create_from_iso8601_string(s) };
		}


		static utctime x_self(const py::tuple& args) {
			return x_arg<utctime>(args, 0);
		}

		static py::object get_seconds(py::tuple args, py::dict kwargs) {
			auto dt = x_self(args).time_since_epoch();
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			if (dt_s == dt)
				return py::object(int64_t(dt_s.count()));
			return py::object(to_seconds(dt));
		}

		static py::object get_float_seconds(py::tuple args, py::dict kwargs) {
			return py::object(to_seconds(x_self(args).time_since_epoch()));
		}

		static py::object get_int_seconds(py::tuple args, py::dict kwargs) {
			auto dt = x_self(args).time_since_epoch();
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			if (dt_s == dt)
				return py::object(int64_t(dt_s.count()));
			return py::object(int64_t(to_seconds(dt)));
		}
		static py::object repr(py::tuple args, py::dict kwargs) {
			auto dt = x_self(args).time_since_epoch();
			auto dt_s = std::chrono::duration_cast<std::chrono::seconds>(dt);
			char s[100];
			if (dt_s == dt)
#ifndef _WIN32
				sprintf(s, "UtcTime(%ld)", dt_s.count());
#else
                sprintf(s, "UtcTime(%lld)", dt_s.count());
#endif
			else
				sprintf(s, "UtcTime(%0.6lf)", to_seconds(dt));
			return py::str(std::string(s));
		}
		static py::object str(py::tuple args, py::dict kwargs) {
			return py::str(calendar().to_string(x_self(args)));
		}
		static py::object floor(py::tuple args, py::dict kwargs) {
			utctime t=x_self(args);
			utctimespan dt=x_arg<utctimespan>(args,1);
			return py::object(utctime_floor(t, dt));
		}
		/// best effort convert from any number or string to utctime 
		/// as long as it is well defined
		static utctime as_utctime(const py::object& po, bool throw_on_fail=true) {
			py::extract<utctime> o(po);
			if (o.check())
				return o();
			py::extract<int64_t> i(po);
			if (i.check())
				return utctime{ from_seconds(i()) };
			py::extract<double> d(po);
			if (d.check())
				return  utctime{ from_seconds(d()) };
			py::extract<string> s(po);
			if (s.check()) 
				return create_from_iso8601_string(s);
			if(throw_on_fail)
				throw std::runtime_error("supplied argument not convertible to UtcTime");
			return no_utctime;
		}
		static double as_seconds(const py::object& po) {
			utctime o = as_utctime(po, false);
			if (o != no_utctime)
				return to_seconds(o.time_since_epoch());
			return to_seconds(utctimespan_ext::as_timespan(po));
		}

		static py::object _lt_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args,0) < as_utctime(args[1]));}
		static py::object _le_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) <= as_utctime(args[1]));}
		static py::object _gt_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) > as_utctime(args[1]));}
		static py::object _ge_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) >= as_utctime(args[1]));}
		static py::object _eq_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) == as_utctime(args[1]));}
		static py::object _nq_(py::tuple args, py::dict kwargs) {args_check(args);return py::object(x_arg<utctime>(args, 0) != as_utctime(args[1]));}
		// math (yes, utctime , we want it to be a number)
		static double arg_as_double(py::object o) {
			utctime rhs = as_utctime(o, false);
			double rhs_d = 0;
			if (rhs != no_utctime) {
				rhs_d = to_seconds(rhs.time_since_epoch());
			} else {
				rhs_d = to_seconds(utctimespan_ext::as_timespan(o));
			}
			return rhs_d;
		}

		static py::object _add_(py::tuple args, py::dict kwargs) { args_check(args); return py::object(x_arg<utctime>(args, 0) + utctimespan_ext::as_timespan(args[1])); }
		static py::object _sub_(py::tuple args, py::dict kwargs) { 
			args_check(args); 
			py::extract<utctime> is_utctime(args[1]);
			if (is_utctime.check()) {
				return py::object(x_arg<utctime>(args, 0) - is_utctime()); // -> converts to TimeSpan
			}
			// if other is time -> then result is timespan
			// else other is timespan -> result is time
			return py::object(utctime{ x_arg<utctime>(args, 0).time_since_epoch() - utctimespan_ext::as_timespan(args[1]) });
		}
		static py::object _rsub_(py::tuple args, py::dict kwargs) { 
			args_check(args); 
			py::extract<utctime> is_utctime(args[1]);
			if (is_utctime.check()) {
				return py::object(is_utctime() - x_arg<utctime>(args, 0)); // -> converts to TimeSpan
			}
			// if other is time -> then result is timespan
			// else other is timespan -> result is time
			return py::object(utctime{ -(x_arg<utctime>(args, 0).time_since_epoch() - utctimespan_ext::as_timespan(args[1])) });
		}
		static py::object _mult_(py::tuple args, py::dict kwargs) { 
			args_check(args);
			double rhs_d =arg_as_double(args[1]);
			return py::object(utctime{ from_seconds(to_seconds(x_arg<utctime>(args, 0).time_since_epoch()) * rhs_d) });
		}
		static py::object _div_(py::tuple args, py::dict kwargs) { 
			args_check(args); 
			double rhs_d = arg_as_double(args[1]);
			return py::object(utctime{ from_seconds(to_seconds(x_arg<utctime>(args, 0).time_since_epoch()) / rhs_d) });
		}
		static py::object _rdiv_(py::tuple args, py::dict kwargs) {
			args_check(args);
			double rhs_d = arg_as_double(args[1]);
			return py::object(utctime{ from_seconds(rhs_d/to_seconds(x_arg<utctime>(args, 0).time_since_epoch())) });
		}

		static py::object _mod_(py::tuple args, py::dict kwargs) { 
			args_check(args); 
			double rhs_d = arg_as_double(args[1]);
			return py::object(utctime{ x_arg<utctime>(args, 0).time_since_epoch() % from_seconds(rhs_d) });
		}
		static py::object _floordiv_(py::tuple args, py::dict kwargs) {
			args_check(args);
			return py::object(utctime {from_seconds(truncf(arg_as_double(args[0]) / arg_as_double(args[1])))});
		}
		static py::object _rfloordiv_(py::tuple args, py::dict kwargs) { 
			args_check(args); 
			return py::object(utctime{ from_seconds(truncf(arg_as_double(args[1]) / arg_as_double(args[0]))) });
		}

	};

	

    static void e_utctime() {
		class_<utctime>("UtcTime", doc_intro("utctime is a timepoint, represented as seconds since epoch 1970.01.01 UTC"),no_init)
			.def("__init__", make_constructor(&utctime_ext::create_default,
				default_call_policies()
			    ),
				doc_intro("construct a UtcTime, with 0 offset from epoch")
			)
			.def("__init__", make_constructor(&utctime_ext::create_from_int,
				default_call_policies(),
				(py::arg("seconds"))
				),
				doc_intro("construct a from seconds as integer")
				doc_parameters()
				doc_parameter("seconds", "int", "seconds")
			)
			.def("__init__", make_constructor(&utctime_ext::create_from_double,
				default_call_policies(),
				(py::arg("seconds"))
				),
				doc_intro("construct a from seconds as decimal number since epoch, where fractions is fractions of second")
				doc_intro("- the resulting time-span preserves 1 micro-second digits")
				doc_parameters()
				doc_parameter("seconds", "int", "seconds")
			)
			.def("__init__", make_constructor(&utctime_ext::create_from_timespan,
				default_call_policies(),
				(py::arg("seconds"))
				),
				doc_intro("construct a from seconds as decimal number since epoch, where fractions is fractions of second")
				doc_intro("- the resulting time-span preserves 1 micro-second digits")
				doc_parameters()
				doc_parameter("seconds", "TimeSpan", "seconds since epoch represented as TimeSpan")
			)

			.def("__init__", make_constructor(&utctime_ext::create_from_string,
				default_call_policies(),
				(py::arg("iso8601_str"))
				),
				doc_intro("construct a from iso 8601s string, e.g. '2014-11-12T19:12:14.505Z'")
				doc_intro("- the resulting time-span preserves 1 micro-second digits")
				doc_parameters()
				doc_parameter("seconds", "int", "seconds")
			)

			.add_property("seconds", raw_function(utctime_ext::get_seconds, 1), doc_intro("returns seconds since epoch"))
			.def("__repr__", raw_function(utctime_ext::repr, 1), doc_intro("repr of UtcTime"))
			.def("__str__", raw_function(utctime_ext::str, 1), doc_intro("str of UtcTime"))
			.def("floor",raw_function(utctime_ext::floor,1),
				//(py::arg("self"),py::arg("dt")),
				doc_intro("floor(t,dt)")
			)
			.def("__int__", raw_function(utctime_ext::get_int_seconds, 1),
				doc_intro("seconds since epoch 1970utc")
			)
			.def("__long__", raw_function(utctime_ext::get_int_seconds, 1),
				doc_intro("seconds since epoch 1970utc")
			)
			.def("__float__",raw_function(utctime_ext::get_float_seconds,1),
				doc_intro("seconds since epoch 1970utc")
			)

			// math operations
			.def(self == self)
			.def(self != self)
			.def(self + utctimespan())
			.def(utctimespan()+ self)
			.def(self - utctimespan())
			.def(self - self)
			.def("__lt__",raw_function(utctime_ext::_lt_,1))
			.def("__le__", raw_function(utctime_ext::_le_, 1))
			.def("__gt__", raw_function(utctime_ext::_gt_, 1))
			.def("__ge__", raw_function(utctime_ext::_ge_, 1))
			.def("__eq__", raw_function(utctime_ext::_eq_, 1))
			.def("__nq__", raw_function(utctime_ext::_nq_, 1))

			.def("__add__", raw_function(utctime_ext::_add_, 1))
			.def("__radd__", raw_function(utctime_ext::_add_, 1))
			.def("__sub__", raw_function(utctime_ext::_sub_, 1))
			.def("__rsub__", raw_function(utctime_ext::_rsub_, 1))
			.def("__truediv__", raw_function(utctime_ext::_div_, 1))
			.def("__rtruediv__", raw_function(utctime_ext::_rdiv_, 1))
			.def("__floordiv__", raw_function(utctime_ext::_floordiv_, 1))
			.def("__rfloordiv__", raw_function(utctime_ext::_rfloordiv_, 1))
			.def("__mod__", raw_function(utctime_ext::_mod_, 1))

			.def("__mul__", raw_function(utctime_ext::_mult_, 1))
			.def("__rmul__", raw_function(utctime_ext::_mult_, 1))


			;

        def("utctime_now",utctime_now,"returns utc-time now as seconds since 1970s");
        def("deltahours",deltahours,args("n"),"returns timespan equal to specified n hours");
        def("deltaminutes",deltaminutes,args("n"),"returns timespan equal to specified n minutes");
        def("is_npos",is_npos,args("n"),"returns true if n is npos, - meaning no position");
        scope current;
        current.attr("max_utctime")= max_utctime;
        current.attr("min_utctime")= min_utctime;
        current.attr("no_utctime")=no_utctime;
        current.attr("npos")=string::npos;
    }
	BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(calendar_time_overloads, calendar::time, 1, 6);
	BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(calendar_time_overloads_week, calendar::time_from_week, 1, 6);
	struct cal_ext {
		static py::object cal_units(py::tuple args, py::dict kwargs) {
			py::extract<shared_ptr<shyft::core::calendar>> cx(args[0]);
			if (cx.check()) {
				return py::object( cx()->calendar_units(utctime_ext::as_utctime(args[1])) );
			}
			throw std::runtime_error("Not able to get shared calendar from self");
		}
		static py::object cal_week_units(py::tuple args, py::dict kwargs) {
			py::extract<shared_ptr<shyft::core::calendar>> cx(args[0]);
			if (cx.check()) {
				return py::object(cx()->calendar_week_units(utctime_ext::as_utctime(args[1])));
			}
			throw std::runtime_error("Not able to get shared calendar from self");
		}
	};
	static void e_calendar() {

		std::string(shyft::core::calendar::*to_string_t)(shyft::core::utctime) const = &calendar::to_string;//selects correct ptr.
		std::string(calendar::*to_string_p)(utcperiod) const = &calendar::to_string;
		int64_t(calendar::*diff_units)(utctime, utctime, utctimespan) const = &calendar::diff_units;
		utctime(calendar::*time_YMDhms)(YMDhms) const = &calendar::time;
		utctime(calendar::*time_YWdhms)(YWdhms) const = &calendar::time;
		utctime(calendar::*time_6)(int, int, int, int, int, int) const = &calendar::time;
		utctime(calendar::*time_from_week_6)(int, int, int, int, int, int) const = &calendar::time_from_week;
		class_<calendar, shared_ptr<calendar>>("Calendar",
			doc_intro("Calendar deals with the concept of human calendar")
			doc_intro(" In SHyFT we practice the 'utctime-perimeter' principle,")
			doc_intro("  * so the core is utc-time only ")
			doc_intro("  * we deal with time-zone and calendars at the interfaces/perimeters")
			doc_intro(" in python, this corresponds to timestamp[64], or as the integer version of the time package representation")
			doc_intro(" e.g. the difference between time.time()- utctime_now() is in split-seconds")
			doc_intro("\n")
			doc_intro("Calendar functionality:")
			doc_intro(" -# Conversion between the calendar coordinates YMDhms or iso week YWdhms and utctime, taking  any timezone and DST into account\n")
			doc_intro(" -# Calendar constants, utctimespan like values for Year,Month,Week,Day,Hour,Minute,Second\n")
			doc_intro(" -# Calendar arithmetic, like adding calendar units, e.g. day,month,year etc.\n")
			doc_intro(" -# Calendar arithmetic, like trim/truncate a utctime down to nearest timespan/calendar unit. eg. day\n")
			doc_intro(" -# Calendar arithmetic, like calculate difference in calendar units(e.g days) between two utctime points\n")
			doc_intro(" -# Calendar Timezone and DST handling\n")
			doc_intro(" -# Converting utctime to string and vice-versa\n")
			doc_intro("\n")
			doc_notes()
			doc_note(" Please notice that although the calendar concept is complete")
			doc_note(" we only implement features as needed in the core and interfaces")
			doc_note(" Currently this includes most options, including olson time-zone handling")
			doc_note(" The time-zone support is currently a snapshot of rules ~2014")
			doc_note(" but we plan to use standard packages like Howard Hinnant's online approach for this later.")
			)

			.def(py::init<utctimespan>((py::arg("tz_offset")),
				doc_intro("creates a calendar with constant tz-offset")
				doc_parameters()
				doc_parameter("tz_offset", "TimeSpan", "seconds utc offset, 3600 gives UTC+01 zone")

				)
			)
			.def(py::init<int>((py::arg("tz_offset")),
				doc_intro("creates a calendar with constant tz-offset")
				doc_parameters()
				doc_parameter("tz_offset", "int", "seconds utc offset, 3600 gives UTC+01 zone")
				)
			)

			.def(py::init<string>((py::arg("olson_tz_id")),
				doc_intro("create a Calendar from Olson timezone id, eg. 'Europe/Oslo'")
				doc_parameters()
				doc_parameter("olson_tz_id", "str", "Olson time-zone id, e.g. 'Europe/Oslo'")
				)
			)
			.def("to_string", to_string_t, (py::arg("self"), py::arg("utctime")),
				doc_intro("convert time t to readable iso standard string taking ")
				doc_intro(" the current calendar properties, including timezone into account")
				doc_parameters()
				doc_parameter("utctime", "int", "seconds utc since epoch")
				doc_returns("iso time string", "str", "iso standard formatted string,including tz info")
			)
			.def("to_string", to_string_p, (py::arg("self"), py::arg("utcperiod")),
				doc_intro("convert utcperiod p to readable string taking current calendar properties, including timezone into account")
				doc_parameters()
				doc_parameter("utcperiod", "UtcPeriod", "An UtcPeriod object")
				doc_returns("period-string", "str", "[start..end>, iso standard formatted string,including tz info")
			)
			.def("region_id_list", &calendar::region_id_list,
				doc_intro("Returns a list over predefined Olson time-zone identifiers")
				doc_notes()
				doc_note("the list is currently static and reflects tz-rules approx as of 2014")
			).staticmethod("region_id_list")
			.def("calendar_units",raw_function(cal_ext::cal_units, 1),
				doc_intro("returns YMDhms for specified t, in the time-zone as given by the calendar")
				doc_parameters()
				doc_parameter("t", "int|UtcTime", "timestamp utc seconds since epoch")
				doc_returns("calendar_units", "YMDhms", "calendar units as in year-month-day hour-minute-second")

			)
			.def("calendar_week_units", raw_function(cal_ext::cal_week_units, 1),
				doc_intro("returns iso YWdhms for specified t, in the time-zone as given by the calendar")
				doc_parameters()
				doc_parameter("t", "int", "timestamp utc seconds since epoch")
				doc_returns("calendar_week_units", "YWdms", "calendar units as in iso year-week-week_day hour-minute-second")
			)
			.def("time", time_YMDhms, args("YMDhms"),
				doc_intro("convert calendar coordinates into time using the calendar time-zone")
				doc_parameters()
				doc_parameter("YMDhms", "YMDhms", "calendar cooordinate structure containg year,month,day, hour,minute,second")
				doc_returns("timestamp", "int", "timestamp as in seconds utc since epoch")
			)
			.def("time", time_YWdhms, args("YWdhms"),
				doc_intro("convert calendar iso week coordinates structure into time using the calendar time-zone")
				doc_parameters()
				doc_parameter("YWdhms", "YWdhms", "structure containg iso specification calendar coordinates")
				doc_returns("timestamp", "int", "timestamp as in seconds utc since epoch")
			)

			.def("time", time_6, calendar_time_overloads(
				doc_intro("convert calendar coordinates into time using the calendar time-zone")
				doc_parameters()
				doc_parameter("Y", "int", "Year")
				doc_parameter("M", "int", "Month  [1..12], default=1")
				doc_parameter("D", "int", "Day    [1..31], default=1")
				doc_parameter("h", "int", "hour   [0..23], default=0")
				doc_parameter("m", "int", "minute [0..59], default=0")
				doc_parameter("s", "int", "second [0..59], default=0")
				doc_returns("timestamp", "int", "timestamp as in seconds utc since epoch")
				,
				args("Y", "M", "D", "h", "m", "s")
			)
			)
			.def("time_from_week", time_from_week_6,
				calendar_time_overloads_week(
					doc_intro("convert calendar iso week coordinates into time using the calendar time-zone")
					doc_parameters()
					doc_parameter("Y", "int", "ISO Year")
					doc_parameter("W", "int", "ISO Week  [1..54], default=1")
					doc_parameter("wd", "int", "week_day  [1..7]=[mo..su], default=1")
					doc_parameter("h", "int", "hour   [0..23], default=0")
					doc_parameter("m", "int", "minute [0..59], default=0")
					doc_parameter("s", "int", "second [0..59], default=0")
					doc_returns("timestamp", "int", "timestamp as in seconds utc since epoch")
					,
					args("Y", "W", "wd", "h", "m", "s")
				)
			)

			.def("add", &calendar::add, args("t", "delta_t", "n"),
				doc_intro("calendar semantic add")
				doc_intro("conceptually this is similar to t + deltaT*n")
				doc_intro(" but with deltaT equal to calendar::DAY,WEEK,MONTH,YEAR")
				doc_intro(" and/or with dst enabled time-zone the variation of length due to dst")
				doc_intro(" or month/year length is taken into account")
				doc_intro(" e.g. add one day, and calendar have dst, could give 23,24 or 25 hours due to dst.")
				doc_intro(" similar for week or any other time steps.")
				doc_parameters()
				doc_parameter("t", "int", "timestamp utc seconds since epoch")
				doc_parameter("delta_t", "int", "timestep in seconds, with semantic interpretation of DAY,WEEK,MONTH,YEAR")
				doc_parameter("n", "int", "number of timesteps to add")
				doc_returns("t", "int", "new timestamp with the added time-steps, seconds utc since epoch")
				doc_see_also("diff_units(t1,t2,delta_t),trim(t,delta_t)")
			)
			.def("diff_units", diff_units, args("t1", "t2", "delta_t"),
				doc_intro("calculate the distance t1..t2 in specified units, taking dst into account if observed")
				doc_intro("The function takes calendar semantics when delta_t is calendar::DAY,WEEK,MONTH,YEAR,")
				doc_intro("and in addition also dst if observed.")
				doc_intro("e.g. the diff_units of calendar::DAY over summer->winter shift is 1,")
				doc_intro("even if the number of hours during those days are 23 and 25 summer and winter transition respectively")
				doc_intro("returns: calendar semantics of (t2-t1)/deltaT, where deltaT could be calendar units DAY,WEEK,MONTH,YEAR")
				doc_parameters()
				doc_parameter("t", "int", "timestamp utc seconds since epoch")
				doc_parameter("delta_t", "int", "timestep in seconds, with semantic interpretation of DAY,WEEK,MONTH,YEAR")
				doc_parameter("n", "int", "number of timesteps to add")
				doc_returns("n_units", "int", "number of units, so that t2 = c.add(t1,delta_t,n) + remainder(discarded)")
				doc_see_also("add(t,delta_t,n),trim(t,delta_t)")

			)
			.def("trim", &calendar::trim, args("t", "delta_t"),
				doc_intro("round time t down to the nearest calendar time-unit delta_t")
				doc_intro("taking the calendar time-zone and dst into account")
				doc_parameter("t", "int", "timestamp utc seconds since epoch")
				doc_parameter("delta_t", "int", "timestep in seconds, with semantic interpretation of Calendar.(DAY,WEEK,MONTH,YEAR)")
				doc_returns("t", "int", "new trimmed timestamp, seconds utc since epoch")
				doc_see_also("add(t,delta_t,n),diff_units(t1,t2,delta_t)")
			)
			.def("quarter", &calendar::quarter, args("t"),
				doc_intro("returns the quarter of the specified t, -1 if invalid t")
				doc_parameters()
				doc_parameter("t", "int", "timestamp utc seconds since epoch")
				doc_returns("quarter", "int", "in range[1..4], -1 if invalid time")
			)
			.def_readonly("YEAR", &calendar::YEAR, "defines a semantic year")
			.def_readonly("MONTH", &calendar::MONTH, "defines a semantic calendar month")
			.def_readonly("QUARTER", &calendar::QUARTER, "defines a semantic calendar quarter (3 months)")
			.def_readonly("DAY", &calendar::DAY, "defines a semantic calendar day")
			.def_readonly("WEEK", &calendar::WEEK, "defines a semantic calendar week")
			.def_readonly("HOUR", &calendar::HOUR, "hour, 3600 seconds")
			.def_readonly("MINUTE", &calendar::MINUTE, "minute, 60 seconds")
			.def_readonly("SECOND", &calendar::SECOND, "second, 1 second")
			.add_property("tz_info", &calendar::get_tz_info, "The TzInfo keeping the time-zone name, utc-offset and DST rules (if any)")//,return_value_policy<return_internal_reference>())
			;

		class_<YMDhms>("YMDhms", "Defines calendar coordinates as Year Month Day hour minute second")
			.def(init<int, optional<int, int, int, int, int>>(args("Y", "M", "D", "h", "m", "s"), "Creates calendar coordinates specifying Y,M,D,h,m,s"))
			.def("is_valid", &YMDhms::is_valid, "returns true if YMDhms values are reasonable")
			.def("is_null", &YMDhms::is_null, "returns true if all values are 0, - the null definition")
			.def_readwrite("year", &YMDhms::year)
			.def_readwrite("month", &YMDhms::month)
			.def_readwrite("day", &YMDhms::day)
			.def_readwrite("hour", &YMDhms::hour)
			.def_readwrite("minute", &YMDhms::minute)
			.def_readwrite("second", &YMDhms::second)
			.def("max", &YMDhms::max, "returns the maximum representation").staticmethod("max")
			.def("min", &YMDhms::max, "returns the minimum representation").staticmethod("min")
			.def(self == self)
			.def(self != self)
			;

		class_<YWdhms>("YWdhms", "Defines calendar coordinates as iso Year Week week-day hour minute second")
			.def(init<int, optional<int, int, int, int, int>>(args("Y", "W", "wd", "h", "m", "s"),
				doc_intro("Creates calendar coordinates specifying iso Y,W,wd,h,m,s")
				doc_parameters()
				doc_parameter("Y", "int", "iso-year")
				doc_parameter("W", "int", "iso week [1..53]")
				doc_parameter("wd", "int", "week_day [1..7]=[mo..sun]")
				doc_parameter("h", "int", "hour [0..23]")
				doc_parameter("m", "int", "minute [0..59]")
				doc_parameter("s", "int", "second [0..59]")
				)
			)
			.def("is_valid", &YWdhms::is_valid, "returns true if YWdhms values are reasonable")
			.def("is_null", &YWdhms::is_null, "returns true if all values are 0, - the null definition")
			.def_readwrite("iso_year", &YWdhms::iso_year)
			.def_readwrite("iso_week", &YWdhms::iso_week)
			.def_readwrite("week_day", &YWdhms::week_day, doc_intro("week_day,[1..7]=[mo..sun]"))
			.def_readwrite("hour", &YWdhms::hour)
			.def_readwrite("minute", &YWdhms::minute)
			.def_readwrite("second", &YWdhms::second)
			.def("max", &YWdhms::max, "returns the maximum representation").staticmethod("max")
			.def("min", &YWdhms::max, "returns the minimum representation").staticmethod("min")
			.def(self == self)
			.def(self != self)
			;

		class_<time_zone::tz_info_t, bases<>, time_zone::tz_info_t_, boost::noncopyable>("TzInfo",
			"TzInfo class is responsible for providing information about the\n"
			" time-zone of the calendar.\n"
			"  This include the\n"
			"   * name (olson identifier),\n"
			"   * base_offset\n"
			"   * utc_offset(t) time-dependent\n"
			"The Calendar class provides a shared pointer to it's TzInfo object \n", no_init
			)
			.def(init<utctimespan>(args("base_tz"), "creates a TzInfo with a fixed utc-offset(no dst-rules)"))
			.def("name", &time_zone::tz_info_t::name, "returns the olson time-zone identifier or name for the TzInfo")
			.def("base_offset", &time_zone::tz_info_t::base_offset, "returnes the time-invariant part of the utc-offset")
			.def("utc_offset", &time_zone::tz_info_t::utc_offset, args("t"), "returns the utc_offset at specified utc-time, takes DST into account if applicable")
			.def("is_dst", &time_zone::tz_info_t::is_dst, args("t"), "returns true if DST is observed at given utc-time t")
			;
	}
    void calendar_and_time() {
		e_utctimespan();
        e_utctime();
        e_utcperiod();
        e_calendar();
    }
}
