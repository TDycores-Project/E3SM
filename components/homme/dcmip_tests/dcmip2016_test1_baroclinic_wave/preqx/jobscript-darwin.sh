#!/bin/bash
#
#   Jobscript for launching dcmip2016 test 1 on a mac running Darwin
#
# usage: ./jobscript-...

# launch the simulation
EXEC=../../../test_execs/preqx-nlev30-interp/preqx-nlev30-interp        # set name of executable
openmpiexec -n 6 $EXEC < ./namelist-lowres.nl                           # launch simulation
