
#
# Python install (still some line required here)
#
    set(Python_ADDITIONAL_VERSIONS 3.6)
    find_package(PythonInterp REQUIRED)
    find_package(PythonLibs REQUIRED)
    set(lib_python ${PYTHON_LIBRARIES})
    message("lib_python: " ${lib_python})
    set(include_python ${PYTHON_INCLUDE_DIRS})
    message("include_python: " ${include_python})
#find_package(NumPy REQUIRED)  # yes would be nice but does not exist yet!
    execute_process(
    COMMAND ${PYTHON_EXECUTABLE} "-c" "import numpy as np; print(np.get_include())"
    RESULT_VARIABLE rv
    OUTPUT_VARIABLE include_numpy
    ERROR_VARIABLE  include_numpy_error
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if (rv)
    message("NumPy include:" ${include_numpy})
    message("NumPy import error:" ${include_numpy_error})
    message(FATAL_ERROR "Errors occurred.  Leaving now!")
    return()
    endif()
    message("include_numpy: " ${include_numpy})

#
# define the shyft core api that exposes all common classes and methods
#
    set(core_api  "api")
    set(core_api_sources 
        api_actual_evapotranspiration.cpp     api_kirchner.cpp
        api_cell_environment.cpp                   api_precipitation_correction.cpp
        api_dtss.cpp                                        api_priestley_taylor.cpp
        api_gamma_snow.cpp                         api_pt_gs_k.cpp
        api_geo_cell_data.cpp                         api_region_environment.cpp
        api_geo_point.cpp                               api_routing.cpp
        api_glacier_melt.cpp                           api.cpp
        api_skaugen.cpp                                 api_hbv_actual_evapotranspiration.cpp  
        api_state.cpp                                      api_kalman.cpp
        api_hbv_physical_snow.cpp                api_target_specification.cpp
        api_hbv_snow.cpp                               api_time_axis.cpp
        api_hbv_soil.cpp                                  api_time_series.cpp
        api_hbv_tank.cpp                                api_utctime.cpp
        api_interpolation.cpp                           api_vectors.cpp
     )

    add_library(${core_api} SHARED ${core_api_sources} )
    target_include_directories(${core_api} PRIVATE  ${include_python} ${include_numpy})
    set_target_properties(${core_api} PROPERTIES OUTPUT_NAME ${core_api})
    set_target_properties(${core_api} PROPERTIES PREFIX "_" INSTALL_RPATH "$ORIGIN/../lib") # Python extensions do not use the 'lib' prefix
    set_property(TARGET ${core_api} APPEND PROPERTY COMPILE_DEFINITIONS SHYFT_EXTENSION)
    target_link_libraries(${core_api} shyft_core shyft_api dlib::dlib ${boost_py_link_libraries} ${arma_libs})
    install(TARGETS ${core_api} DESTINATION ${CMAKE_SOURCE_DIR}/shyft/api)

#
# define each shyft method-stack that exposes complete stacks
#
    foreach(shyft_stack  "pt_gs_k" "pt_ss_k" "pt_hs_k" "hbv_stack" "pt_hps_k")
        add_library(${shyft_stack} SHARED  ${shyft_stack}.cpp)
        target_include_directories(${shyft_stack} PRIVATE  ${include_python} ${include_numpy} )
        set_target_properties(${shyft_stack} PROPERTIES OUTPUT_NAME ${shyft_stack})
        set_target_properties(${shyft_stack} PROPERTIES PREFIX "_" INSTALL_RPATH "$ORIGIN/../../lib") # Python extensions do not use the 'lib' prefix
        set_property(TARGET ${shyft_stack} APPEND PROPERTY COMPILE_DEFINITIONS SHYFT_EXTENSION)
        target_link_libraries(${shyft_stack} shyft_core shyft_api  dlib::dlib ${boost_py_link_libraries} ${ARMADILLO_LIBRARIES} )
        install(TARGETS ${shyft_stack} DESTINATION ${CMAKE_SOURCE_DIR}/shyft/api/${shyft_stack})
    endforeach(shyft_stack)

#
# install 3rd party .so files ( but we still require blas/lapack + gcc 7.x runtime)
# we place them in the lib directory parallell to the shyft directory, and notice the INSTALL_RPATH
# for the targets do refer to this location. If you change any of these, just keep in sync 
#  (if not the you will get reminded when runtime fail to load referenced shared libraries)
#
    # step 1: pull out abs-paths from imported targets
    get_property(dlib_lib TARGET dlib::dlib PROPERTY IMPORTED_LOCATION_RELEASE) 
    # not yet supported on cmake 3.9.x:
    #get_property(boost_system_lib TARGET Boost::system PROPERTY IMPORTED_LOCATION_RELEASE)
    #get_property(boost_filesystem_lib TARGET Boost::filesystem PROPERTY IMPORTED_LOCATION_RELEASE)
    #get_property(boost_serialization_lib  TARGET Boost::serialization PROPERTY IMPORTED_LOCATION_RELEASE)
    #get_property(boost_python3_lib TARGET Boost::python3 PROPERTY IMPORTED_LOCATION_RELEASE)
    #get_target_property(arma_lib armadillo IMPORTED_LOCATION_NOCONFIG)
    # step 2: get the real-name (we want libboost_system.so.1.66.0, not the symlink libboost_system.so)
    get_filename_component(a_lib ${ARMADILLO_LIBRARIES} REALPATH CACHE)
    get_filename_component(d_lib ${dlib_lib} REALPATH CACHE)
    get_filename_component(b_filesystem ${boost_filesystem_lib} REALPATH CACHE)
    get_filename_component(b_system ${boost_system_lib} REALPATH CACHE)
    get_filename_component(b_serialization ${boost_serialization_lib} REALPATH CACHE)
    get_filename_component(b_python3 ${boost_python3_lib} REALPATH CACHE)

    # step 3: install to the shyft/lib (that is referenced by the RPATH relative to the extensions libraries on linux)
    install(FILES 
        ${a_lib}
        ${d_lib}
        ${b_system}
        ${b_filesystem}
        ${b_python3}
        ${b_serialization}
        DESTINATION ${CMAKE_SOURCE_DIR}/shyft/lib
    )
    #install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ${$CMAKE_SOURCE_DIR}/shyft/lib/libdlib.so.8)")
