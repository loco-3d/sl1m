CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

# Project properties
SET(PROJECT_ORG loco-3d)
SET(PROJECT_NAME sl1m)
SET(PROJECT_DESCRIPTION "Implementation of the sl1m solver for multi contact planning")
SET(PROJECT_URL "https://github.com/${PROJECT_ORG}/${PROJECT_NAME}")

# Project configuration
SET(PROJECT_USE_CMAKE_EXPORT TRUE)

# JRL-cmakemodule setup
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/python.cmake)

# Project definition
COMPUTE_PROJECT_ARGS(PROJECT_ARGS LANGUAGES CXX)
PROJECT(${PROJECT_NAME} ${PROJECT_ARGS})

FINDPYTHON()

SET(${PROJECT_NAME}_SOURCES
  __init__.py
  constants_and_tools.py
  solver.py
  generic_solver.py
  fix_sparsity.py
  problem_data.py
  constraints_gait.py
  constraints.py
  constraints_biped.py
  planner_generic_gait.py
  planner_generic.py
  planner_biped.py
  problem_definition.py
  problem_definition_gait.py
  )

SET(${PROJECT_NAME}_PLANNER_SCENARIOS_SOURCES
  __init__.py
  )

SET(${PROJECT_NAME}_PLANNER_SCENARIOS_HRP2_SOURCES
  __init__.py
  problem_definition_hrp2.py
  complex.py
  path_files/lp_complex_path.py
  path_files/lp_stair_bauzil_hrp2_path.py
  path_files/lp_stair_bauzil_hrp2_path_10.py
  )

SET(${PROJECT_NAME}_PLANNER_SCENARIOS_TALOS_SOURCES
  __init__.py
  complex.py
  stairs.py
  rubbles.py
  path_files/lp_complex_path.py
  path_files/lp_ramp_path.py
  path_files/lp_rubbles_path.py
  path_files/lp_stairs_path.py
  )

SET(${PROJECT_NAME}_RBPRM_SOURCES
  __init__.py
  narrow_convex_hull.py
  constants_and_tools.py
  surfaces_from_planning.py
  )

SET(${PROJECT_NAME}_STANDALONE_SCENARIOS_SOURCES
  __init__.py
  hrp2_complex.py
  hrp2_stairs.py
  talos_ramp.py
  talos_rubble_stairs.py
  problem_definition_talos.py
  problem_definition_hrp2.py
  )

  SET(${PROJECT_NAME}_STANDALONE_SURFACES_SOURCES
  __init__.py
  stair_surfaces.py
  complex_surfaces.py
  flat_ground.py
  rubble_stair_surfaces.py
  ramp_surfaces.py
  )
SET(${PROJECT_NAME}_TOOLS_SOURCES
  __init__.py
  geom_utils.py
  obj_to_constraints.py
  plot_utils.py
  plot_tools.py
  heightmap_tools.py
  transformations.py
  )

FOREACH(src ${${PROJECT_NAME}_SOURCES})
  PYTHON_INSTALL_ON_SITE(${PROJECT_NAME} ${src})
ENDFOREACH(src ${PROJECT_NAME}_SOURCES)

FOREACH(src ${${PROJECT_NAME}_PLANNER_SCENARIOS_SOURCES})
  PYTHON_INSTALL_ON_SITE(${PROJECT_NAME}/planner_scenarios ${src})
ENDFOREACH(src ${PROJECT_NAME}_PLANNER_SCENARIOS_SOURCES)

FOREACH(src ${${PROJECT_NAME}_PLANNER_SCENARIOS_TALOS_SOURCES})
  PYTHON_INSTALL_ON_SITE(${PROJECT_NAME}/planner_scenarios/talos ${src})
ENDFOREACH(src ${PROJECT_NAME}_PLANNER_SCENARIOS_TALOS_SOURCES)

FOREACH(src ${${PROJECT_NAME}_PLANNER_SCENARIOS_HRP2_SOURCES})
  PYTHON_INSTALL_ON_SITE(${PROJECT_NAME}/planner_scenarios/hrp2 ${src})
ENDFOREACH(src ${PROJECT_NAME}_PLANNER_SCENARIOS_HRP2_SOURCES)

FOREACH(src ${${PROJECT_NAME}_RBPRM_SOURCES})
  PYTHON_INSTALL_ON_SITE(${PROJECT_NAME}/rbprm ${src})
ENDFOREACH(src ${PROJECT_NAME}_RBPRM_SOURCES)

FOREACH(src ${${PROJECT_NAME}_STANDALONE_SURFACES_SOURCES})
  PYTHON_INSTALL_ON_SITE(${PROJECT_NAME}/stand_alone_scenarios/surfaces ${src})
ENDFOREACH(src ${PROJECT_NAME}_STANDALONE_SURFACES_SOURCES)

FOREACH(src ${${PROJECT_NAME}_STANDALONE_SCENARIOS_SOURCES})
  PYTHON_INSTALL_ON_SITE(${PROJECT_NAME}/stand_alone_scenarios ${src})
ENDFOREACH(src ${PROJECT_NAME}_STANDALONE_SCENARIOS_SOURCES)

FOREACH(src ${${PROJECT_NAME}_TOOLS_SOURCES})
  PYTHON_INSTALL_ON_SITE(${PROJECT_NAME}/tools ${src})
ENDFOREACH(src ${PROJECT_NAME}_TOOLS_SOURCES)

INSTALL(DIRECTORY sl1m/stand_alone_scenarios/constraints_files
  DESTINATION ${PYTHON_SITELIB}/sl1m/stand_alone_scenarios)

ADD_LIBRARY(${PROJECT_NAME} INTERFACE)
INSTALL(TARGETS ${PROJECT_NAME} EXPORT ${TARGETS_EXPORT_NAME} DESTINATION lib)
