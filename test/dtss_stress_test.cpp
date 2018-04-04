#include "test_pch.h"

#include "core/dtss.h"
#include "core/dtss_cache.h"
#include "core/dtss_client.h"

#include "core/utctime_utilities.h"
#include "core/time_axis.h"
#include "core/time_series_merge.h"
#include "core/time_series_dd.h"


#include <future>
#include <mutex>
#include <regex>
#include <cstdint>
#include <cstdlib>
#include <shared_mutex>
#include <boost/filesystem.hpp>
namespace  fs=boost::filesystem;

#include <core/dtss_client.h>

using namespace std;
using namespace shyft;
using namespace shyft::core;
using shyft::time_series::dd::apoint_ts;
using shyft::time_series::dd::gta_t;


namespace dtss_stress {

    apoint_ts mk_expression(utctime t, utctimespan dt, int n) {
        std::vector<double> x; x.reserve(n);
        for (int i = 0; i < n; ++i)
            x.push_back(-double(n) / 2.0 + i);
        apoint_ts aa(gta_t(t, dt, n), x);
        auto a = aa*3.0 + aa;
        return a;
    }

    string ts_url(const string&c,const string ts) {return "shyft://"+c+"/"+ts;}
    string server_addr(const string & host,int port) {return host+":"+to_string(port);}
    string mk_ts_url(const string&container,const string& ts_path,int i) {return ts_url(container,ts_path+"/"+to_string(i));}
    struct config {
        string host{"127.0.0.1"};
        int port{20000};
        utctimespan dt_min_read{10000};
        utctimespan dt_read_var{1000};
        string container{"test"};
        string ts_path{"a"};
        size_t min_ts{10};
        size_t max_ts{100};
        size_t min_n{10};
        size_t max_n{100};
        utctimespan dt{10};
        size_t n_pts_min{1};
        int n_pts_var{10000};
    };

    struct reader {
        config c{};
        bool use_cache{ true };
        bool update_cache{ false };
        dtss::ts_vector_t mk_read_expression() {
            dtss::ts_vector_t r;
            int n = c.min_ts + rand()%(c.max_ts-c.min_ts);
            r.reserve(n);for(size_t i=0;i<n;++i) r.push_back(apoint_ts(mk_ts_url(c.container,c.ts_path,i)));
            return r;
        }

        utcperiod mk_read_period(utctime t_now) {
            utctime t_end = t_now - rand()%1;
            utctime t_start = t_end - (c.dt_min_read + rand()%c.dt_read_var);
            return utcperiod{t_start,t_end};
        }

        size_t do_work(utctime t_exit) {
            auto t_now = utctime_now();
            dtss::client d(server_addr(c.host,c.port));
            size_t rr=0;
            size_t rn=0;
            while(t_now < t_exit) {
                t_now =utctime_now();
                auto rp=mk_read_period(t_now);
                int n = c.min_n + rand()%(c.max_n-c.min_n);
                gta_t ta(rp.start,utctimespan{1 + rp.timespan()/n},n);
                try {
                    use_cache = (rand()%100) > 50;  // randomize 
                    auto r = d.percentiles(mk_read_expression(),rp,ta,vector<int64_t>{0,10,50,75,100},use_cache,update_cache);
                    rr++;rn +=r.size();
                    this_thread::sleep_for(chrono::milliseconds(10));
                } catch (const runtime_error &re) {
                    throw runtime_error(string("reader:")+re.what());
                }
            }
            return rr;
        }

    };

    struct writer {
        volatile size_t n{0};
        promise<bool> done;
        bool cache_on_write{ true };
        bool replace_all{ false };
        config c{};
        dtss::ts_vector_t mk_write_ts(utctime t_now) {
            dtss::ts_vector_t r;
            calendar utc;
            size_t n = c.n_pts_min + rand()%c.n_pts_var;
            auto t0= utc.trim(t_now,c.dt)- n*c.dt;
            for(size_t i=0;i<c.max_ts;++i) {
                vector<double> v(n,1.0);
                apoint_ts ts_v(gta_t(t0,c.dt,n),v,time_series::POINT_AVERAGE_VALUE);
                r.push_back(apoint_ts(mk_ts_url(c.container,c.ts_path,i),ts_v));
            }
            return r;
        }

        size_t do_work(utctime t_exit) {
            auto t_now= utctime_now();
            dtss::client d(server_addr(c.host,c.port));
            n = 0;
            while (t_now < t_exit) {
                try {
                    d.store_ts(mk_write_ts(t_now),replace_all, cache_on_write);
                } catch(const runtime_error &ex) {
                    throw runtime_error(string("writer:")+ ex.what());
                }
                if(n==0) {
                    cout<<"w1\n";
                    done.set_value(true);
                }
                ++n;
                cout<<"."; cout.flush();
                this_thread::sleep_for(chrono::milliseconds(20));
                t_now=utctime_now();
            }
            return n;
        }

    };

    struct single_writer_multi_reader {
        vector<reader> r{};
        writer w{};

        void setup(const config &c, size_t n_readers) {
            w.c=c;
            for(size_t i=0;i<n_readers;++i)
                r.push_back(reader{c});
        }

        vector<size_t> run(utctime t_exit) {
            vector<future<size_t>> x;
            //
            //w.do_work(t_exit);
            auto wm =std::async(std::launch::async,[t_exit,this]()->size_t {return w.do_work(t_exit); });
            auto ready=w.done.get_future().wait_for(chrono::seconds(55))!=future_status::timeout;
            if(!ready && w.n< 1)
                throw runtime_error("still missing writer ready");
            std::cout<<"Ok, launch readers\n";
            for(size_t i=0;i<r.size();++i) {
                x.push_back(
                    std::async(std::launch::async, [t_exit,this,i]()->size_t{return r[i].do_work(t_exit);})
                );
            }
            std::cout<<"wait for readers\n";

            vector<size_t> c;

            for (auto &f : x)
                    c.push_back(f.get());
            std::cout<<"FINALE WRITER\n";
            wm.get();
            std::cout<<"DONE\n";
            std::cout.flush();
            return c;
        }
    };
}


TEST_SUITE("dtss") {

TEST_CASE("dtss_stress") {
    using namespace dtss_stress;
    config c;
    auto tmpdir = fs::temp_directory_path()/"shyft.stress";
    dtss::server<dtss::standard_dtss_dispatcher> srv{};
    srv.add_container(c.container,tmpdir.string());
    srv.set_listening_ip(c.host);
    srv.set_listening_port(c.port);
    srv.start_async();
    bool rr=srv.is_running();
    FAST_REQUIRE_EQ(rr,true);
    single_writer_multi_reader swmr{};
    size_t n_readers=1;
    int n_seconds =20;
    if(getenv("SHYFT_DTSS_STRESS")) {
        if(sscanf(getenv("SHYFT_DTSS_STRESS"),"%zu,%d",&n_readers,&n_seconds)!=2) {
            FAST_REQUIRE_EQ(true,false);
        }
    }
    swmr.setup(c,n_readers);
    auto x=swmr.run(utctime_now()+n_seconds);
    auto cs =srv.get_cache_stats();
    srv.clear();
    FAST_CHECK_GT(x.size(),0);
}
}
