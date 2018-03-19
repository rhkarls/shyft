# CMake file for compiling the C++ core library
# library sources
set(SOURCES   dtss_client.cpp dtss.cpp experimental.cpp dream_optimizer.cpp sceua_optimizer.cpp 
    time_series_dd.cpp core_serialization.cpp time_series_serialization.cpp
    utctime_utilities.cpp expression_serialization.cpp
)
message(STATUS "Generating Makefile for the SHyFT static library...")
add_library(shyftcore STATIC ${SOURCES} )
set_target_properties(shyftcore PROPERTIES OUTPUT_NAME shyftcore)
if (MSVC)
  set_target_properties(shyftcore PROPERTIES PREFIX lib)
endif()
target_link_libraries(shyftcore ${LIBS})
set(core_headers
    actual_evapotranspiration.h      hbv_snow_common.h           pt_hs_k.h
    bayesian_kriging.h               hbv_snow.h                  pt_ss_k_cell_model.h
    cell_model.h                     hbv_soil.h                  pt_ss_k.h
    core_archive.h                   hbv_stack_cell_model.h      region_model.h
    core_serialization.h             hbv_stack.h                 routing.h
    dream_optimizer.h                hbv_tank.h                  sceua_optimizer.h
    dtss_cache.h                     inverse_distance.h          skaugen.h
    dtss_client.h                    kalman.h                    time_axis.h
    dtss_db.h                        kirchner.h                  time_series_average.h
    dtss.h                           kriging.h                   time_series_common.h
    dtss_msg.h                       model_calibration.h         time_series_dd.h
    dtss_mutex.h                     model_state_tuning.h        time_series.h
    dtss_url.h                       optimizer_utils.h           time_series_info.h
    experimental.h                   precipitation_correction.h  time_series_merge.h
    expression_serialization.h       predictions.h               time_series_point_merge.h
    gamma_snow.h                     priestley_taylor.h          time_series_qm.h
    geo_cell_data.h                  pt_gs_k_cell_model.h        time_series_statistics.h
    geo_point.h                      pt_gs_k.h                   unit_conversion.h
    glacier_melt.h                   pt_hps_k_cell_model.h       utctime_utilities.h
    hbv_actual_evapotranspiration.h  pt_hps_k.h
    hbv_physical_snow.h              pt_hs_k_cell_model.h
)
install(TARGETS shyftcore ARCHIVE DESTINATION lib COMPONENT shyft_dist EXCLUDE_FROM_ALL)
install(FILES ${core_headers}  DESTINATION include/shyft/core COMPONENT shyft_dist EXCLUDE_FROM_ALL)
