struct_lib = static_library('struct_lib',
  files(
    'dynamic_correlator.cc',
    'fourier.cc',
    'plaq.cc',
    'ptetra.cc',
    'qsi.cc',
    'spin.cc',
    'static_correlator.cc',
    'tetra.cc',
    'util.cc'
    ),
  include_directories: inc,
  dependencies: [gsl_dep, fftw3_dep]
  )

lib += [struct_lib]
