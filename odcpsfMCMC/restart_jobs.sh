#!/bin/bash                                                                                                                    
if [ $1 ]; then
    walltimeString='-l walltime='$1
else
    walltimeString=''
fi

find ~/projectspace/odcpsf/jobs -type f -name "myBatchFile.sh" -exec qsub {}  $walltimeString \;
