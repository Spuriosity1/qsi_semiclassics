
common_lib = static_library('common_lib',
  files(
    'defect_strategy.cc',
    'direction.cc'
    ),
  include_directories: inc,
  dependencies: [fftw3_dep, gsl_dep]
  )

lib += [common_lib]
