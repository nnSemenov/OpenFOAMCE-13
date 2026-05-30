find_package(wmakeParse
    0.2.2
    COMPONENTS wmakeParse parse-wmake
    CONFIG
    REQUIRED
)

function(wmake_parse_files file OUT_VAR)
    unset(${OUT_VAR} PARENT_SCOPE)

    get_target_property(parser_loc wmakeParse::parse-wmake LOCATION)

    execute_process(
        COMMAND ${parser_loc} ${file} --strict --print-files
        COMMAND_ERROR_IS_FATAL ANY
        OUTPUT_VARIABLE output
    )

    string(REPLACE "\n" ";" output ${output})
    set(${OUT_VAR} ${output} PARENT_SCOPE)
endfunction()