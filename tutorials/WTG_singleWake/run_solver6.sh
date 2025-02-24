#!/bin/bash
cp -r 0.orig 0
blockMesh
setExprFields
topoSet > log.topoAD
decomposePar 
mpirun -np 6 renumberMesh -overwrite -parallel 
mpirun -np 6 simpleFoam -parallel > log.solver1
reconstructPar -latestTime
for i in {0..5}; do
  if [ -d "processor${i}" ]; then
    rm -rf "processor${i}"
  fi
done
if [ ! -d "simLogs" ]; then
  mkdir simLogs
fi
mv log.topoAD ./simLogs
mv log.solver1 ./simLogs
foamToVTK -cellSet disk1


