cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH ) ### MIGRATE-ACTION-REQUIRED (migrate-3.20.00) - avoid: use find_package() with custom Find<pkg>.cmake for Spack compatibility
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH ) ### MIGRATE-ACTION-REQUIRED (migrate-3.20.00) - avoid: use find_package() with custom Find<pkg>.cmake for Spack compatibility

add_subdirectory(Common)
if(SBND)
  add_subdirectory(SBND)
endif()
if(ICARUS)
  add_subdirectory(ICARUS)
endif()
