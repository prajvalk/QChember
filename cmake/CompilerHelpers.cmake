# Helper module to configure compiler flags across multiple toolchains

# Public function: configure_compiler_flags()
function(configure_compiler_flags)
  # Normalize build types and set defaults
  set(VALID_BUILD_TYPES Release Debug RelWithDebInfo CACHE STRING "")
  if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build" FORCE)
  endif()

  # Base flags we will modify per-compiler
  set(_cflags "")
  set(_cxxflags "")
  set(_ldflags "")

  # Detect compiler
  message(STATUS "Detected C compiler: ${CMAKE_C_COMPILER_ID}")
  message(STATUS "Detected CXX compiler: ${CMAKE_CXX_COMPILER_ID}")

  # Default per-config values (will be appended/overridden per-compiler)
  set(CMAKE_C_FLAGS_DEBUG "-O1 -g")
  set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g")

  set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g")

  set(CMAKE_C_FLAGS_RELEASE "-O3 -Ofast")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -Ofast")

  # Compiler-specific tuning
  if(CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
    # Use -march=native where appropriate (user may override)
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -flto")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -flto")

    set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -march=native -flto")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -flto")

    # Debug should keep optimization low but still include debug info
    set(CMAKE_C_FLAGS_DEBUG "-O1 -g -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -fno-omit-frame-pointer")

    # LTO linker flags (gcc/clang)
    set(_ldflags "-flto")

  elseif(CMAKE_C_COMPILER_ID STREQUAL "IntelLLVM")
    # Intel compilers prefer -O3 and -ipo for link-time optimization
    set(CMAKE_C_FLAGS_RELEASE "-O3 -xHost -ipo")
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -xHost -ipo")

    set(CMAKE_C_FLAGS_RELWITHDEBINFO "-O3 -g -xHost -ipo")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -g -xHost -ipo")

    set(CMAKE_C_FLAGS_DEBUG "-O1 -g")
    set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g")

    set(_ldflags "-ipo")

  elseif(CMAKE_C_COMPILER_ID STREQUAL "MSVC")
    # MSVC uses different flag syntax
    set(CMAKE_C_FLAGS_RELEASE "/O2 /DNDEBUG /GL")
    set(CMAKE_CXX_FLAGS_RELEASE "/O2 /DNDEBUG /GL")

    set(CMAKE_C_FLAGS_RELWITHDEBINFO "/O2 /Zi /GL")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/O2 /Zi /GL")

    # For Debug, request low optimization but include debugging info
    set(CMAKE_C_FLAGS_DEBUG "/O1 /Zi")
    set(CMAKE_CXX_FLAGS_DEBUG "/O1 /Zi")

    # MSVC LTO happens at link stage via /LTCG (CMake will generally add it when /GL present)
    set(_ldflags "")

  elseif(CMAKE_C_COMPILER_ID STREQUAL "AMD")
    # AMD AOCC and clang-based compilers are similar to clang/gcc in flags
    set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -march=native -flto")
    set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -flto")

    set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} -march=native -flto")
    set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} -march=native -flto")

    set(CMAKE_C_FLAGS_DEBUG "-O1 -g -fno-omit-frame-pointer")
    set(CMAKE_CXX_FLAGS_DEBUG "-O1 -g -fno-omit-frame-pointer")

    set(_ldflags "-flto")

  else()
    message(WARNING "Unknown compiler: ${CMAKE_C_COMPILER_ID}. Falling back to conservative flags.")
  endif()

  # Apply any linker flags (global)
  if(_ldflags)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${_ldflags}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${_ldflags}")
  endif()

  # Setup OpenMP: find package but do not enable automatically for Debug
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_FLAGS}")
    # We don't force adding flags here — we provide apply_common_target_flags() to opt-in per-target.
    set(ENABLE_OPENMP_IN_RELEASE_ONLY TRUE)
  else()
    message(STATUS "OpenMP not found; continuing without OpenMP support")
    set(ENABLE_OPENMP_IN_RELEASE_ONLY FALSE)
  endif()

  # Expose convenience function to apply flags to targets
  function(apply_common_target_flags target)
    # Apply generic compile options and per-config flags using generator expressions
    if(TARGET ${target})
      # Generic language standards (example)
      target_compile_features(${target} PUBLIC cxx_std_17)

      # Per-configuration flags: use generator expressions to choose
      target_compile_options(${target} PRIVATE
        $<$<CONFIG:Debug>:${CMAKE_CXX_FLAGS_DEBUG}>
        $<$<CONFIG:RelWithDebInfo>:${CMAKE_CXX_FLAGS_RELWITHDEBINFO}>
        $<$<CONFIG:Release>:${CMAKE_CXX_FLAGS_RELEASE}>
      )

      target_link_options(${target} PRIVATE
        $<$<CONFIG:Debug>:${CMAKE_EXE_LINKER_FLAGS}>
        $<$<CONFIG:RelWithDebInfo>:${CMAKE_EXE_LINKER_FLAGS}>
        $<$<CONFIG:Release>:${CMAKE_EXE_LINKER_FLAGS}>
      )

      # Optionally add OpenMP only for Release and RelWithDebInfo
      if(OpenMP_CXX_FOUND)
        target_compile_options(${target} PRIVATE
          $<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>:${OpenMP_CXX_FLAGS}>
        )
        if(OpenMP_CXX_LIB_NAMES)
          target_link_libraries(${target} PRIVATE
            $<$<OR:$<CONFIG:Release>,$<CONFIG:RelWithDebInfo>>:${OpenMP_CXX_LIB_NAMES}>
          )
        endif()
      endif()

      # Expose macros for code
      target_compile_definitions(${target} PRIVATE
        $<$<CONFIG:Debug>:BUILD_DEBUG>
        $<$<CONFIG:RelWithDebInfo>:BUILD_RELWITHDEBINFO>
        $<$<CONFIG:Release>:BUILD_RELEASE>
      )

    else()
      message(FATAL_ERROR "apply_common_target_flags: target ${target} does not exist")
    endif()
  endfunction()

  # Expose chosen flags as cache variables so they can be inspected or overridden
  set(PROJECT_C_FLAGS "${CMAKE_C_FLAGS}" CACHE STRING "Project C flags")
  set(PROJECT_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING "Project CXX flags")
endfunction()
