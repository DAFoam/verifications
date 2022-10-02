#!/bin/bash

# Check if the OpenFOAM enviroments are loaded
if [ -z "$WM_PROJECT" ]; then
  echo "OpenFOAM environment not found, forgot to source the OpenFOAM bashrc?"
  exit
fi

# generate mesh
echo "Generating mesh.."
blockMesh
topoSet
echo "Generating mesh.. Done!"

# copy initial and boundary condition files
cp -r 0.orig 0
