#include "test_pch.h"

#include "core/dtss_mutex.h"

TEST_SUITE("dtss") {
TEST_CASE("dtss_mutex") {
    using namespace shyft::dtss;
    file_lock_manager m;
    string a("a");
    string b("b");
    {
        reader_file_lock f1(m,a);
        FAST_CHECK_EQ(m.f_locks.size(),1);
    }//
    FAST_CHECK_EQ(m.f_locks.size(),0);
    {
        reader_file_lock f2(m,a);
        reader_file_lock f3(m,b);
        FAST_CHECK_EQ(m.f_locks.size(),2);
    }
    FAST_CHECK_EQ(m.f_locks.size(),0);
    {
        writer_file_lock f4(m,a);
        reader_file_lock f5(m,b);
        FAST_CHECK_EQ(m.f_locks.size(),2);
    }
    FAST_CHECK_EQ(m.f_locks.size(),0);
}
}
