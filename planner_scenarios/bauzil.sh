#! /bin/bash


gepetto-gui &
hpp-rbprm-server &
cp lp_stair_bauzil.py lp_stair_bauzil_hrp2_path.py /media/data/dev/linux/hpp/src/hpp-rbprm-corba/script/scenarios/demos
#~ ipython -i --no-confirm-exit $DEVEL_HPP_DIR/src/multicontact-locomotion-planning/scripts/run_mlp.py talos_circle
ipython -i --no-confirm-exit $DEVEL_HPP_DIR/src/multicontact-locomotion-planning/scripts/run_mlp.py lp_stair_bauzil

pkill -f  'gepetto-gui'
pkill -f  'hpp-rbprm-server'

