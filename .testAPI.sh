#!/usr/bin/env bash
  
if [ -z "$WM_PROJECT" ]; then
  echo "OpenFOAM environment not found, forgot to source the OpenFOAM bashrc?"
  exit 1
fi

echo "Running.."
#find */runScript* -type f -exec sed -i '/"primalMinResTol"/c\    "primalMinResTol": 0.9,' {} \;
cd DARhoSimpleCFoam && ./preProcessing.sh && python runScript.py --mode=reverse --task=runAdjoint && python runScript.py --mode=forward --task=runForwardAD --dvName="shape" --seedIndex=0 && cd - || exit 1
cd DARhoSimpleFoam && ./preProcessing.sh && python runScript.py --mode=reverse --task=runAdjoint && python runScript.py --mode=forward --task=runForwardAD --dvName="shape" --seedIndex=0 && cd - || exit 1
cd DASimpleFoam && ./preProcessing.sh && python runScript.py --mode=reverse --task=runAdjoint && python runScript.py --mode=forward --task=runForwardAD --dvName="shape" --seedIndex=0 && cd - || exit 1
cd DASimpleFoamField && ./preProcessing.sh && python runScript.py --mode=reverse --task=runAdjoint && python runScript.py --mode=forward --task=runForwardAD --dvName="beta" --seedIndex=0 && cd - || exit 1 
