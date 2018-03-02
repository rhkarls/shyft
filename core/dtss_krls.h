#pragma once

#include <string>
#include <vector>

#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

#include "core/time_series_dd.h"
#include "time_series_info.h"
#include "utctime_utilities.h"


namespace shyft {
namespace dtss {

class krls_pred_db {

    using gta_t = shyft::time_axis::generic_dt;
    using gts_t = shyft::time_series::point_ts<gta_t>;

    std::string root_dir;

public:
    krls_pred_db() = default;
    ~krls_pred_db() = default;
  
    /** Constructs a krls_pred_db with specified container root */
    explicit krls_pred_db(const std::string& root_dir)
        : root_dir{ root_dir }
    {
        if ( ! fs::is_directory(root_dir) ) {
            if ( ! fs::exists(root_dir) ) {
                if ( ! fs::create_directories(root_dir) ) {
                    throw std::runtime_error(std::string("krls_pred_db: failed to create root directory: ") + root_dir);
                }
            } else {
                throw std::runtime_error(std::string("krls_pred_db: designated root directory is not a directory: ") + root_dir);
            }
        }
    }

    krls_pred_db(const krls_pred_db &) = default;
    krls_pred_db(krls_pred_db &&) = default;

    // krls_pred_db & operator=(const krls_pred_db &);  // TODO wait for file handles to close! (windows)
    // krls_pred_db & operator=(krls_pred_db &&);  // TODO wait for file handles to close! (windows)


    /*  Container API
     * =============== */

public:
    void save(const std::string & fn, const gts_t & ts, bool overwrite = true, bool win_thread_close = true) const {
        // TODO
    }

    gts_t read(const std::string & fn, core::utcperiod p) const {
        // TODO
    }

    void remove(const std::string & fn) const {
        // TODO
    }

    ts_info get_ts_info(const std::string & fn) const {
        // TODO
    }

    std::vector<ts_info> find(const std::string & match) const {
        // TODO
    };


    /*  Internal implementation
     * ========================= */

// private:


};

}
}
