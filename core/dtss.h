#pragma once

#include <boost/serialization/vector.hpp>
#include <boost/serialization/shared_ptr.hpp>
#include <boost/serialization/base_object.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

namespace shyft {
    namespace dtss {

        enum message_type {
            SERVER_EXCEPTION,
            EVALUATE_TS_VECTOR,
            EVALUATE_TS_VECTOR_PERCENTILES,
        };

        // TODO: Reformulate this crap using serialization concepts
        namespace msg {
            template <class T>
            message_type read_type(T&in) {
                int32_t mtype;
                in.read((char*)&mtype,sizeof(mtype));
                return (message_type) mtype;
            }

            template <class T>
            void write_type(message_type mt,T&out) {
                int32_t mtype=(int32_t)mt;
                out.write((char const *)&mtype,sizeof(mtype));
            }



            template <class T>
            void write_string(std::string const&s,T& out) {
                int32_t sz=s.size();
                out.write((char const*)&sz,sizeof(sz));
                out.write(s.data(),sz);
            }

            template<class T>
            void write_exception(std::exception const & e,T& out) {
                int32_t sz=strlen(e.what());
                out.write((char const*)&sz,sizeof(sz));
                out.write(e.what(),sz);
            }

            template<class T>
            void send_exception(std::exception const & e,T& out) {
                write_type(message_type::SERVER_EXCEPTION,out);
                int32_t sz=strlen(e.what());
                out.write((char const*)&sz,sizeof(sz));
                out.write(e.what(),sz);
            }

            template<class T>
            std::runtime_error read_exception(T& in) {
                int32_t sz;
                in.read((char*)&sz,sizeof(sz));
                std::string msg(sz,'\0');
                in.read((char*)msg.data(),sz);
                return std::runtime_error(msg);
            }

        }

        typedef std::vector<api::apoint_ts> ts_vector_t;
        typedef std::vector<string> id_vector_t;
        typedef std::function< ts_vector_t (id_vector_t ts_ids,core::utcperiod p)> call_back_t;

        struct server : dlib::server_iostream {
            call_back_t bind_ts_cb;

            template <class CB>
            server(CB&& cb):bind_ts_cb(std::forward<CB>(cb)) {
            }
            ~server() {}


            ts_vector_t
            do_evaluate_ts_vector(core::utcperiod bind_period, ts_vector_t& atsv) {

                std::map<std::string,std::vector<api::ts_bind_info>> ts_bind_map;
                std::vector<std::string> ts_id_list;
                // step 1: bind not yet bound time-series ( ts with only symbol, needs to be resolved using bind_cb)
                for (auto& ats : atsv) {
                    auto ts_refs = ats.find_ts_bind_info();
                    for(const auto& bi:ts_refs) {
                        if (ts_bind_map.find(bi.reference) == ts_bind_map.end()) { // maintain unique set
                            ts_id_list.push_back( bi.reference );
                            ts_bind_map[bi.reference] = std::vector<api::ts_bind_info>();
                        }
                        ts_bind_map[bi.reference].push_back(bi);
                    }
                }

                // step 2: (optional) bind_ts callback should resolve symbol time-series with content
                if(ts_bind_map.size()) {
                    auto bts=bind_ts_cb(ts_id_list,bind_period);
                    if(bts.size()!=ts_id_list.size())
                        throw std::runtime_error(std::string("failed to bind all of ")+std::to_string(bts.size())+std::string(" ts"));

                    for(size_t i=0;i<ts_id_list.size();++i) {
                        for(auto &bi:ts_bind_map[ts_id_list[i]])
                            bi.ts.bind(bts[i]);
                    }
                }
                //-- step 3: evaluate the ts-vector (later: check if percentiles first..)
                ts_vector_t evaluated_tsv;
                for (auto &ats : atsv) { // TODO: in parallel
                    evaluated_tsv.emplace_back(ats.time_axis(), ats.values(), ats.point_interpretation());
                }
                return evaluated_tsv;
            }




            void on_connect(
                std::istream& in,
                std::ostream& out,
                const std::string& foreign_ip,
                const std::string& local_ip,
                unsigned short foreign_port,
                unsigned short local_port,
                dlib::uint64 connection_id
            ) {
                while (in.peek() != EOF) {
                    auto msg_type= msg::read_type(in);
                    try {
                        switch (msg_type) {
                            case EVALUATE_TS_VECTOR: {
                                boost::archive::binary_iarchive ia(in);
                                core::utcperiod bind_period;
                                ts_vector_t rtsv;
                                ia>>bind_period>>rtsv;
                                msg::write_type(message_type::EVALUATE_TS_VECTOR,out);
                                boost::archive::binary_oarchive oa(out);
                                oa<<do_evaluate_ts_vector(bind_period, rtsv);
                            } break;
                            default:
                                throw std::runtime_error(std::string("Got unknown message type:") + std::to_string((int)msg_type));
                        }
                    } catch (std::exception const& e) {
                        msg::send_exception(e,out);
                    }
                }
            }
        };

        inline std::vector<api::apoint_ts> dtss_evaluate(std::string host_port, std::vector<api::apoint_ts> const& tsv, core::utcperiod p) {
            dlib::iosockstream io(host_port);
            msg::write_type(message_type::EVALUATE_TS_VECTOR,io);
            boost::archive::binary_oarchive oa(io);
            oa<<p<<tsv;
            auto response_type= msg::read_type(io);
            if(response_type==message_type::SERVER_EXCEPTION) {
                auto re= msg::read_exception(io);
                throw re;
            } else if(response_type==message_type::EVALUATE_TS_VECTOR) {
                boost::archive::binary_iarchive ia(io);
                ts_vector_t r;
                ia>>r;
                return r;
            }
            throw std::runtime_error(std::string("Got unexpected response:")+std::to_string((int)response_type));
        }
    }
}
