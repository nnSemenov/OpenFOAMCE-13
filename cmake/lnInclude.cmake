function(lnInclude dir)
    cmake_parse_arguments(arg "" "RESULT_DIR;HEADER_LIST;SOURCE_LIST" "EXCLUDE_DIRS" ${ARGN})

    cmake_path(ABSOLUTE_PATH dir)

    set(lnInclude_dir "${dir}/lnInclude")
    if(arg_RESULT_DIR)
        set(${arg_RESULT_DIR} ${lnInclude_dir} PARENT_SCOPE)
    endif ()

    cmake_path(RELATIVE_PATH dir BASE_DIRECTORY ${CMAKE_SOURCE_DIR} OUTPUT_VARIABLE disp_dir)


    file(GLOB_RECURSE headers "${dir}/*.H")
    file(GLOB_RECURSE h_headers "${dir}/*.h")
    file(GLOB_RECURSE sources "${dir}/*.C")
    set(full_sources "${headers};${h_headers};${sources}")

    if(arg_HEADER_LIST)
        set(${arg_HEADER_LIST} "${headers};${h_headers}" PARENT_SCOPE)
    endif ()
    if(arg_SOURCE_LIST)
        set(${arg_SOURCE_LIST} ${sources} PARENT_SCOPE)
    endif ()

    if(IS_DIRECTORY ${lnInclude_dir})
        message(STATUS "lnInclude ${disp_dir} ... Already finished.")
        return()
    else ()
        message(STATUS "lnInclude ${disp_dir}")
    endif ()

    file(MAKE_DIRECTORY ${lnInclude_dir})

    # Some projects have sub projects. Don't lnInclude files from subprojects
    set(exclude_dir_list)
    foreach (exclude_dir ${arg_EXCLUDE_DIRS})
        cmake_path(ABSOLUTE_PATH exclude_dir)
        list(APPEND exclude_dir_list ${exclude_dir})
    endforeach ()
    #    message("Absolute exclude dir: ${exclude_dir_list}")

    foreach (file ${full_sources})
        if(${file} MATCHES "lnInclude")
            continue()
        endif ()

        set(skip_this_file OFF)
        foreach (exclude_dir ${exclude_dir_list})
            cmake_path(IS_PREFIX exclude_dir ${file} is_prefix_result)
            if (${is_prefix_result})
                set(skip_this_file ON)
                break()
            endif ()
        endforeach ()
        if (${skip_this_file})
            continue()
        endif ()

        cmake_path(GET file FILENAME link_name)
        set(link_name "${lnInclude_dir}/${link_name}")
        cmake_path(RELATIVE_PATH file BASE_DIRECTORY ${lnInclude_dir} OUTPUT_VARIABLE link_target)
#        message("${link_name} -> ${link_target}")
        file(CREATE_LINK ${link_target} ${link_name} SYMBOLIC RESULT result)
#        if(NOT result)
#            message(WARNING "Failed to create symlink ${link_name} -> ${file} with result ${result}")
#        endif ()
    endforeach ()
endfunction()