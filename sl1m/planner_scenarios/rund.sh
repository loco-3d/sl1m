#!/bin/bash         

gepetto-gui &
ipython -i --no-confirm-exit ./$1

pkill -f  'gepetto-gui'
