lib = []

subdir('common')
subdir('manager')
subdir('struct')

fs = import('fs')

# Common sources (deps = g_deps)
exe_sources = [
  'anneal.cc',
  'average.cc',
  'bb_plot.cc',
  'print_positions.cc',
  'qsi_defect_disorder.cc',
  'qsi_defect_dynamic.cc',
  'qsi_defect_dynamic_from_file.cc',
  'qsi_defect_static.cc',
  'qsi_defect_static_from_file.cc',
  'qsi_dynamic.cc',
  'qsi_static.cc',
]

# Special case: needs json_dep in addition
exe_sources_json = [
  'couplingfile_from_diluted.cc',
]

foreach src : exe_sources
  name = fs.stem(src)   # removes .cc
  executable(name,
    src,
    dependencies: g_deps,
    include_directories: inc,
    link_with: lib,
    install: true,
  )
endforeach

foreach src : exe_sources_json
  name = fs.stem(src)   # removes .cc
  executable(name,
    src,
    dependencies: [g_deps, json_dep],
    include_directories: inc,
    link_with: lib,
    install: true,
  )
endforeach

