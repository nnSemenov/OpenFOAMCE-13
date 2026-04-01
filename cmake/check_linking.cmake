cmake_minimum_required(VERSION 3.28)

set(linking_tool "ldd" CACHE STRING "Tool for link examination")
set(search_dirs "" CACHE PATH "Directories to be searched")

set(root_lib_names "libfiniteVolume;libOpenFOAM;libparallel;libODE;libtracking;libspecie;libsampling;libDSMC;libchemistryModel;libDSMC" CACHE STRING "whitelist of root libs")

function(try_linking exe_loc working_dir)
    execute_process(
        COMMAND ${linking_tool} -r -d ${exe_loc}
        WORKING_DIRECTORY ${working_dir}
        OUTPUT_VARIABLE output
    )

    if((${output} MATCHES "not found") OR (${output} MATCHES "undefined symbol"))
        message(FATAL_ERROR "${exe_loc} failed to be loaded: \n${output}")
    endif ()
endfunction()

foreach (dir ${search_dirs})
    message(STATUS "Checking executable/libraries under ${dir} . Some libraries are skipped.")
    file(GLOB_RECURSE exe_libs LIST_DIRECTORIES false "${dir}/*")

    set(dylib_only ON)
    cmake_path(GET dir FILENAME stem)
    if(${stem} STREQUAL "bin")
        set(dylib_only OFF)
    endif ()

    foreach (exe ${exe_libs})
        cmake_path(GET exe EXTENSION exe_extension)
        cmake_path(GET exe STEM exe_stem)
        cmake_path(GET exe FILENAME exe_filename)

        if(exe_extension STREQUAL ".sh")
            continue()
        endif ()

        # Check if this file should be skipped
        if(${dylib_only}) # lib
            if(NOT ${exe_extension} MATCHES ".so")
                continue()
            endif ()

            if(${exe_stem} MATCHES "^lib[A-Za-z0-9_]+Solver") # in white list, check
#                message(STATUS "Found ${exe_stem}")
            elseif (${exe_stem} IN_LIST root_lib_names) # is solver module, check

            else () # Not a lib that could be dynamically loaded
                continue()
            endif ()
        else ()
            # not dylib only, check this
        endif ()



        message(STATUS "Checking ${exe_filename} ...")
        try_linking(${exe} ${dir})

    endforeach ()
endforeach ()