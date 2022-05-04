# sl1m

[![Pipeline status](https://gitlab.laas.fr/$ORG/sl1m/badges/master/pipeline.svg)](https://gitlab.laas.fr/loco-3d/sl1m/commits/master)
[![Coverage report](https://gitlab.laas.fr/loco-3d/sl1m/badges/master/coverage.svg?job=doc-coverage)](http://projects.laas.fr/gepetto/doc/loco-3d/sl1m/master/coverage/)


## Python dependencies:
	
	python -m pip install scipy
	sudo apt install robotpkg-py38-quadprog
	python -m pip install matplotlib
	python -m pip install pycddlib
	python -m pip install cvxpy

:warning: this package requires gurobi :warning:
A free academic license can be obtained at https://www.gurobi.com

## Demo dependencies:

These are the dependencies for the robot and rbprm server.
```
sudo apt install -y robotpkg-py38-hpp-rbprm \
					robotpkg-py38-hpp-rbprm-corba \
					robotpkg-py38-anymal-rbprm \
					robotpkg-py38-solo-rbprm \
					robotpkg-py38-talos-rbprm \
```

## Run demo:

```
# Source you install prefix.
source ~/devel/workspace_sl1m/install/setup.bash
# Locate where the Solo12 specific data are.
export SOLO3D_ENV_DIR=/home/mnaveau/Documents/cloud_laas/Shared/Solo3D/
# Point toward you HPP installation folder
export INSTALL_HPP_DIR=/opt/openrobots/share
# Add the data folder to the ROS environment for urdf parsing.
export ROS_PACKAGE_PATH=/opt/openrobots/share/example-robot-data/robots:${SOLO3D_ENV_DIR}:$ROS_PACKAGE_PATH
```

Execute the scripts:

```bash
python sl1m/stand_alone_scenarios/anymal_stairs.py
python sl1m/stand_alone_scenarios/anymal_trot.py
python sl1m/stand_alone_scenarios/data.pickle
python sl1m/stand_alone_scenarios/hrp2_complex.py
python sl1m/stand_alone_scenarios/hrp2_stairs.py
python sl1m/stand_alone_scenarios/problem_definition_hrp2.py
python sl1m/stand_alone_scenarios/problem_definition_talos.py
python sl1m/stand_alone_scenarios/solo_stairs.py
python sl1m/stand_alone_scenarios/solo_trot.py
python sl1m/stand_alone_scenarios/talos_ramp.py
python sl1m/stand_alone_scenarios/talos_rubble_stairs.py
```
