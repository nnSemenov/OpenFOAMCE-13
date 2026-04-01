
target_compile_options(OpenFOAM_Defines INTERFACE
    -Wall
    -Wextra
    -Wold-style-cast
    -Wnon-virtual-dtor
    -Wno-unused-parameter
    -Wno-invalid-offsetof
    -Wno-attributes
    -ftemplate-depth-256
    -Wno-array-bounds

    -frounding-math
    -ftrapping-math
)

#set(CMAKE_LINKER_TYPE BFD)

target_link_options(OpenFOAM_Defines INTERFACE
#    -fuse-ld=bfd
#    -Xlinker --add-needed
    -Xlinker --no-as-needed
)

list(APPEND Mikeno_less_warn_options
    -Wno-old-style-cast
    -Wno-unused-local-typedefs
    -Wno-tautological-undefined-compare
    -Wno-shift-negative-value
)
