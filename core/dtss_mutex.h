#pragma once
#include <string>
#include <memory>
#include <unordered_map>
#include <mutex>
#include <shared_mutex>
#include <functional>

namespace shyft {
namespace dtss {

using std::mutex;
using std::unordered_map;
using std::shared_mutex;
using std::recursive_mutex;
using std::shared_lock;
using std::unique_lock;
using std::lock_guard;
using std::string;
using std::hash;
using std::weak_ptr;
using std::make_shared;
using std::shared_ptr;

/** implements a map<filename,single-writer multi-reader>
 *
 * This helps keeping number of handles/mutex down to a
 * minimum simultaneously used resources.
 *
 * The idea is that when accessing a file, this class
 * provides the shared_mutex (single-writer/multi-reader lock) to be used.
 *
 * If it's already in the map, return that resource, otherwise, create a new one
 * and add it to the map.
 *
 * When the ref.count goes to 1, the shared_mutex is removed.
 *
 * How to use: ref. to dtss_db.h
 * This class .get() and .put() methods should be used through the
 * helper classes writer_file_lock and reader_file_lock below exclusively.
 *
 * Notes:
 * 1.
 *  this only helps for multi-threading the dtss. If several dtss processes, that
 *  share the same common shyft-container files, we need to extend this class to also use
 *  inter-processs file-locks so that the files are locked between processes.
 * 2.
 *  if overhead mutex is to large, then the shared_mutex could be re-cycled into a pool.
 */
struct file_lock_manager {
    using xmap=unordered_map<string,shared_ptr<shared_mutex>>;
    recursive_mutex f_locks_mx; ///< the mutex that protects the map below
    xmap f_locks; ///< map of live dyn_mutexes, TODO: remap dt so that we control destructur/unregister
    file_lock_manager()=default;

    /** returns a shared pointer to dyn_mutex for the specified file_path */
    shared_ptr<shared_mutex> get(const string& file_path) {
        lock_guard<recursive_mutex> sl(f_locks_mx); // protect f_locks exclusively
        auto f= f_locks.find(file_path);
        if(f == f_locks.end()) {
            auto r=make_shared<shared_mutex>();// create a new fresh dyn_mutex
            f_locks[file_path] = r;
            return r;
        } else {
            return f->second;// hand out a shared pointer, there is already a thread accessing the file
        }
    }

    /** if the file_path shared_mutex has ref-count 1, remove it from the map, -release handle to keep minimum live handles */
    void put(const string & file_path) {
        lock_guard<recursive_mutex> sl(f_locks_mx); // protect f_locks exclusively
        auto f= f_locks.find(file_path);
        if(f != f_locks.end()) {
            if(f->second.use_count()==1) {
                f_locks.erase(f);
            }
        }
    }
    //- not implemented:
    file_lock_manager(const file_lock_manager&)=delete;
    file_lock_manager(file_lock_manager&&)=delete;
    file_lock_manager& operator=(const file_lock_manager&)=delete;
    file_lock_manager& operator=(file_lock_manager&&)=delete;

};

/** a scoped shared reader file-locking using file_lock_manager
 *
 * To be used in a scope where filename and lock-manager have longer
 * lifetime than this object. ref. to dtss db for example.
 *
 */
struct reader_file_lock {
    shared_ptr<shared_mutex> sl;///< the shared-mutex that does all the heavy lifting
    const string &fname;///< const ref to filename
    file_lock_manager&dm; ///< reference to file-lock manager that provides the shared-mutex managment
    reader_file_lock(file_lock_manager&dm,const string& fname):sl{dm.get(fname)},fname{fname},dm{dm} {
        sl->lock_shared();
    }
    ~reader_file_lock() {
        sl->unlock_shared();
        sl.reset();
        dm.put(fname);
    }
    //-
    reader_file_lock()=delete;
    reader_file_lock(const reader_file_lock&)=delete;
    reader_file_lock(reader_file_lock&&)=delete;
    reader_file_lock& operator=(const reader_file_lock&)=delete;
    reader_file_lock& operator=(reader_file_lock&&)=delete;
};

/** a scoped exclusive writer_file_lock
 *
 * To be used in a scope where filename and lock-manager have longer
 * lifetime than this object. ref. to dtss db for example.
 *
 */
struct writer_file_lock {
    shared_ptr<shared_mutex> sl;
    const string &fname;
    file_lock_manager& dm;
    writer_file_lock(file_lock_manager&dm,const string& fname):sl{dm.get(fname)},fname{fname},dm{dm} {
        sl->lock();
    }
    ~writer_file_lock() {
        sl->unlock();
        sl.reset();
        dm.put(fname);
    }
    //-
    writer_file_lock()=delete;
    writer_file_lock(const writer_file_lock&)=delete;
    writer_file_lock(writer_file_lock&&)=delete;
    writer_file_lock& operator=(const writer_file_lock&)=delete;
    writer_file_lock& operator=(writer_file_lock&&)=delete;

};

}}
