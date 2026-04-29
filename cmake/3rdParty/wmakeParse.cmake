find_package(wmakeParse CONFIG REQUIRED)

function(wmake_parse_files file OUT_VAR)
    unset(${OUT_VAR} PARENT_SCOPE)

    get_target_property(parser_loc wmakeParse::parse-wmake LOCATION)
#    message(STATUS "Found parse-wmake at ${parser_loc}")

    execute_process(
        COMMAND ${parser_loc} ${file} --strict --print-files
        COMMAND_ERROR_IS_FATAL ANY
        OUTPUT_VARIABLE output
    )

    set(${OUT_VAR} ${output} PARENT_SCOPE)
endfunction()