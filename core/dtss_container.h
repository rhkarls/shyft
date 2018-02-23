#pragma once

#include <string>
#include <vector>

#include <boost/variant.hpp>

#include "time_series_info.h"
#include "utctime_utilities.h"
#include "dtss_db.h"


namespace shyft {
namespace dtss {

/** A CRTP interface for DTSS server containers.
 */
template < class IMPL >
struct container {

    using container_t = IMPL;

    using gta_t = shyft::time_axis::generic_dt;
    using gts_t = shyft::time_series::point_ts<gta_t>;

    /** Save a time-series to the container. */
    void save(
        const std::string & tsid,
        const gts_t & ts,
        bool overwrite = true,
        bool win_thread_close = true
    ) const {
        static_cast<const IMPL*>(this)->save(tsid, ts, overwrite, win_thread_close);
    }

    /** Read a period from a time-series from the container. */
    gts_t read(
        const std::string & tsid,
        core::utcperiod period
    ) const {
        return static_cast<const IMPL*>(this)->read(tsid, period);
    }

    /** Remove a time-series from the container. */
    void remove(const std::string & tsid) const {
        static_cast<const IMPL*>(this)->remove(tsid);
    }

    /** Get minimal info about a time-series stored in the container. */
    ts_info get_ts_info(const std::string & tsid) const {
        return static_cast<const IMPL*>(this)->get_ts_info(tsid);
    }

    /** Find minimal information about all time-series stored in the container matching a regex pattern. */
    std::vector<ts_info> find(const std::string & pattern) const {
        return static_cast<const IMPL*>(this)->find(pattern);
    }

};


template < class ... CIMPLs >
struct container_wrapper {

    using gta_t = shyft::time_axis::generic_dt;
    using gts_t = shyft::time_series::point_ts<gta_t>;
    // -----
    using container_adt = boost::variant<container<CIMPLs>...>;
    // -----
    container_adt _container;

    container_wrapper() = default;
    ~container_wrapper() = default;
    // -----
    template < class CIMPL >
    container_wrapper(container<CIMPL> & c) : _container{ container_adt{ c } } { }
    // -----
    container_wrapper(const container_wrapper &) = default;
    container_wrapper & operator=(const container_wrapper &) = default;
    // -----
    container_wrapper(container_wrapper &&) = default;
    container_wrapper & operator=(container_wrapper &&) = default;

    void save(
        const std::string & tsid,
        const gts_t & ts,
        bool overwrite = true,
        bool win_thread_close = true
    ) const {
        boost::apply_visitor([&](auto && var) {
            var.save(tsid, ts, overwrite, win_thread_close);
        }, _container);
    }
    gts_t read(
        const std::string & tsid,
        core::utcperiod period
    ) const {
        return boost::apply_visitor([&](auto && var) -> gts_t {
            return var.read(tsid, period);
        }, _container);
    }
    void remove(const std::string & tsid) const {
        boost::apply_visitor([&](auto && var) {
            return var.remove(tsid);
        }, _container);
    }
    ts_info get_ts_info(const std::string & tsid) const {
        return boost::apply_visitor([&](auto && var) -> ts_info {
            return var.get_ts_info(tsid);
        }, _container);
    }
    std::vector<ts_info> find(const std::string & pattern) const {
        return boost::apply_visitor([&](auto && var) -> std::vector<ts_info> {
            return var.find(pattern);
        }, _container);
    }

};

}
}
