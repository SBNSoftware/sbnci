cet_find_library( BOOST_SERIALIZATION NAMES boost_serialization PATHS ENV BOOST_LIB NO_DEFAULT_PATH )
cet_find_library( BOOST_DATE_TIME     NAMES boost_date_time     PATHS ENV BOOST_LIB NO_DEFAULT_PATH )

execute_process(
    COMMAND bash "
        if [[ ${ref_pfp_hist} == *sbnd* ]]
         then
            echo 'sbnd'
        elif [[ ${ref_pfp_hist} == *icarus* ]]
         then
            echo 'icarus'
        fi"
    OUTPUT_VARIABLE expName
)

add_subdirectory(Common)
#if(${nsbnd} EQUAL 1)
if("${expName}" EQUAL "sbnd")
  add_subdirectory(SBND)
endif()
#if(${nicarus} EQUAL 1)
if("${expName}" EQUAL "icarus")
  add_subdirectory(ICARUS)
endif()
