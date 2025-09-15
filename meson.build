project('qsi_classical_dynamics', ['cpp', 'c'],
  version: '1.2.0',
  meson_version : '>= 1.5.0',
  default_options : [ 
    'warning_level=3',
    'buildtype=release',
    'cpp_std=c++2a',
    'prefix=' + meson.project_source_root(), 
    ],
)



if get_option('buildtype') == 'release'
  add_project_arguments(['-Ofast', '-march=native'],  language : 'cpp')
endif

inc = []
subdir('include')

gsl_dep = dependency('gsl', required: true)
fftw3_dep = dependency('fftw3', required: true)
json_dep = dependency('nlohmann_json', required: true)


g_deps = [
  gsl_dep,
  fftw3_dep
  ]


subdir('src')

