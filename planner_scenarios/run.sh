#!/bin/bash         

gepetto-gui &
hpp-rbprm-server &
ipython3 -i --no-confirm-exit ./$1

pkill -f  'gepetto-gui'
pkill -f  'hpp-rbprm-server'
