#!/usr/bin/env python
"""
DAFoam run script to verify the adjoint derivative accuracy
"""

# =============================================================================
# Imports
# =============================================================================
import os
import argparse
from mpi4py import MPI
from dafoam import PYDAFOAM, optFuncs
from pygeo import *
from pyspline import *
from idwarp import USMesh
import numpy as np

# =============================================================================
# Input Parameters
# =============================================================================
parser = argparse.ArgumentParser()
parser.add_argument("--task", help="type of run to do", type=str, default="runAdjoint")
parser.add_argument("--mode", help="AD mode: either reverse or forward", type=str, default="reverse")
parser.add_argument("--dvName", help="design variable name for forward AD", type=str, default="shape")
parser.add_argument("--seedIndex", help="which design variable index to set seeds", type=int, default=0)
args = parser.parse_args()
gcomm = MPI.COMM_WORLD

MRF0 = -500
UIn = 99.9
# Set the parameters for optimization
daOptions = {
    "solverName": "DARhoSimpleFoam",
    "designSurfaces": ["blade"],
    "primalMinResTol": 1e-12,
    "useAD": {"mode": args.mode},
    "primalBC": {
        "U0": {"variable": "U", "patches": ["inlet"], "value": [0.0, 0.0, UIn]},
        "MRF": MRF0,
    },
    "objFunc": {
        "CMX": {
            "part1": {
                "type": "moment",
                "source": "patchToFace",
                "patches": ["blade"],
                "axis": [0.0, 0.0, 1.0],
                "center": [0.0, 0.0, 0.0],
                "scale": 1.0,
                "addToAdjoint": True,
            }
        },
        "FX": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["blade"],
                "directionMode": "fixedDirection",
                "direction": [0.0, 0.0, 1.0],
                "scale": 1.0,
                "addToAdjoint": True,
            }
        },
        "PL": {
            "part1": {
                "type": "totalPressure",
                "source": "patchToFace",
                "patches": ["inlet"],
                "scale": 1.0,
                "addToAdjoint": True,
            },
            "part2": {
                "type": "totalPressure",
                "source": "patchToFace",
                "patches": ["outlet"],
                "scale": -1.0 ,
                "addToAdjoint": True,
            },
        },
    },
    "primalVarBounds": {
        "UMax": 1000.0,
        "UMin": -1000.0,
        "pMax": 500000.0,
        "pMin": 20000.0,
        "eMax": 500000.0,
        "eMin": 100000.0,
        "rhoMax": 5.0,
        "rhoMin": 0.2,
    },
    "normalizeStates": {"U": 10.0, "p": 10000.0, "nuTilda": 1e-3, "phi": 1.0, "T": 300.0},
    "adjPartDerivFDStep": {"State": 1e-6},
    "adjEqnOption": {"gmresRelTol": 1.0e-10, "pcFillLevel": 1, "jacMatReOrdering": "rcm"},
    "adjPCLag": 5,
    # Design variable setup
    "designVar": {
        "shapey": {"designVarType": "FFD"},
        "shapez": {"designVarType": "FFD"},
        "MRF": {"designVarType": "BC"},
        "U0": {"designVarType": "BC", "patches": ["inlet"], "variable": "U", "comp": 0}
    },
    "decomposeParDict": {"preservePatches": ["per1", "per2"]}
}

# mesh warping parameters, users need to manually specify the symmetry plane
meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "OpenFOAM",
    # point and normal for the symmetry plane
    "symmetryPlanes": [],
}


# =============================================================================
# Design variable setup
# =============================================================================
FFDFile = "./FFD/localFFD.xyz"
DVGeo = DVGeometry(FFDFile)
DVGeo.addRefAxis("bodyAxis", xFraction=0.25, alignIndex="k")
# select points
pts = DVGeo.getLocalIndex(0)
indexList = pts[1:3, :, 1].flatten()
PS = geo_utils.PointSelect("list", indexList)
DVGeo.addLocalDV("shapey", lower=-1.0, upper=1.0, axis="y", scale=1.0, pointSelect=PS)
DVGeo.addLocalDV("shapez", lower=-1.0, upper=1.0, axis="z", scale=1.0, pointSelect=PS)


def MRF(val, geo):
    DASolver.setOption("primalBC", {"MRF": float(val[0])})
    DASolver.updateDAOption()

DVGeo.addGlobalDV("MRF", [MRF0], MRF, lower=-1000.0, upper=1000.0, scale=1.0)

def U0(val, geo):
    inletU = float(val[0])
    DASolver.setOption("primalBC", {"U0": {"variable": "U", "patches": ["inlet"], "value": [0.0, 0.0, inletU]}})
    DASolver.updateDAOption()
DVGeo.addGlobalDV("U0", [UIn], U0, lower=-1000.0, upper=1000.0, scale=1.0)

# =============================================================================
# DAFoam initialization
# =============================================================================
DASolver = PYDAFOAM(options=daOptions, comm=gcomm)
DASolver.setDVGeo(DVGeo)
mesh = USMesh(options=meshOptions, comm=gcomm)
DASolver.addFamilyGroup(DASolver.getOption("designSurfaceFamily"), DASolver.getOption("designSurfaces"))
DASolver.printFamilyList()
DASolver.setMesh(mesh)
evalFuncs = []
DASolver.setEvalFuncs(evalFuncs)

# =============================================================================
# Constraint setup
# =============================================================================
DVCon = DVConstraints()
DVCon.setDVGeo(DVGeo)
DVCon.setSurface(DASolver.getTriangulatedMeshSurface(groupName=DASolver.getOption("designSurfaceFamily")))

# =============================================================================
# Initialize optFuncs for optimization
# =============================================================================
optFuncs.DASolver = DASolver
optFuncs.DVGeo = DVGeo
optFuncs.DVCon = DVCon
optFuncs.evalFuncs = evalFuncs
optFuncs.gcomm = gcomm

# =============================================================================
# Task
# =============================================================================
if args.task == "runAdjoint":

    optFuncs.runAdjoint()

elif args.task == "runForwardAD":

    optFuncs.runForwardAD(args.dvName, args.seedIndex)

else:
    print("task arg not found!")
    exit(0)
