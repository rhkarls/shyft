# Compute the installation prefix relative to this file.(is this really needed ?)
get_filename_component(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
if(_IMPORT_PREFIX STREQUAL "/")
  set(_IMPORT_PREFIX "")
endif()

set(shyft_core_INCLUDE_DIRS ${_IMPORT_PREFIX}/include)
set(shyft_core_LIBRARIES ${_IMPORT_PREFIX}/lib/libshyft_core.a)

