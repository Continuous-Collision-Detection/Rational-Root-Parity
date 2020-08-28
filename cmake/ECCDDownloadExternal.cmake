################################################################################
include(DownloadProject)

# With CMake 3.8 and above, we can hide warnings about git being in a
# detached head by passing an extra GIT_CONFIG option
if(NOT (${CMAKE_VERSION} VERSION_LESS "3.8.0"))
    set(ECCD_EXTRA_OPTIONS "GIT_CONFIG advice.detachedHead=false")
else()
    set(ECCD_EXTRA_OPTIONS "")
endif()

# Shortcut function
function(eccd_download_project name)
    download_project(
        PROJ         ${name}
        SOURCE_DIR   ${THIRD_PARTY_DIR}/${name}
        DOWNLOAD_DIR ${THIRD_PARTY_DIR}/.cache/${name}
        QUIET
        ${ECCD_EXTRA_OPTIONS}
        ${ARGN}
    )
endfunction()

################################################################################

# eigen
function(eccd_download_eigen)
    eccd_download_project(eigen
        GIT_REPOSITORY  https://gitlab.com/libeigen/eigen
        GIT_TAG         25424d91f60a9f858e7dc1c7936021cc1dd72019
    )
endfunction()


function(eccd_download_catch2)
    eccd_download_project(Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v2.10.1
    )
endfunction()