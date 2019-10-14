#! /bin/bash


gepetto-gui &
hpp-rbprm-server &
cp lp_complex.py lp_complex_path.py /media/data/dev/linux/hpp/src/hpp-rbprm-corba/script/scenarios/demos
#~ ipython -i --no-confirm-exit $DEVEL_HPP_DIR/src/multicontact-locomotion-planning/scripts/run_mlp.py talos_circle
ipython -i --no-confirm-exit $DEVEL_HPP_DIR/src/multicontact-locomotion-planning/scripts/run_mlp.py lp_complex

pkill -f  'gepetto-gui'
pkill -f  'hpp-rbprm-server'

