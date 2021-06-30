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
parser.add_argument("--dvName", help="design variable name for forward AD", type=str, default="shape")
parser.add_argument("--seedIndex", help="which design variable index to set seeds", type=int, default=0)
args = parser.parse_args()
gcomm = MPI.COMM_WORLD

# Define the global parameters here
U0 = 10.0
p0 = 0.0
nuTilda0 = 4.5e-5
alpha0 = 3.0
A0 = 0.1

adjJacOpt = "JacobianFree"
mode = "reverse"

if args.task == "runForwardAD":
    adjJacOpt = "JacobianFD"
    mode = "forward"

# Set the parameters for optimization
daOptions = {
    "designSurfaces": ["wing"],
    "solverName": "DASimpleFoam",
    "primalMinResTol": 1.0e-14,
    "adjJacobianOption": adjJacOpt,
    "useAD": {"mode": mode, "dvName": args.dvName, "seedIndex": args.seedIndex},
    # don't bound nuTilda because it will degrade the adjoint derivative accuracy
    "primalVarBounds": {"nuTildaMin": -1e16},
    "primalBC": {
        "U0": {"variable": "U", "patches": ["inout"], "value": [U0, 0.0, 0.0]},
        "p0": {"variable": "p", "patches": ["inout"], "value": [p0]},
        "nuTilda0": {"variable": "nuTilda", "patches": ["inout"], "value": [nuTilda0]},
        "useWallFunction": True,
    },
    "fvSource": {
        "disk1": {
            "type": "actuatorDisk",
            "source": "cylinderAnnulusSmooth",
            "center": [-0.5, 0.0, 0.05],
            "direction": [1.0, 0.0, 0.0],
            "innerRadius": 0.01,
            "outerRadius": 0.4,
            "rotDir": "right",
            "scale": 10.0,
            "POD": 0.8,
            "eps": 0.1,
            "expM": 1.0,
            "expN": 0.5,
            "adjustThrust": 1,
            "targetThrust": 0.2,
        },
    },
    "objFunc": {
        "CD": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wing"],
                "directionMode": "parallelToFlow",
                "alphaName": "alpha",
                "scale": 1.0 / (0.5 * U0 * U0 * A0),
                "addToAdjoint": True,
            }
        },
        "CL": {
            "part1": {
                "type": "force",
                "source": "patchToFace",
                "patches": ["wing"],
                "directionMode": "normalToFlow",
                "alphaName": "alpha",
                "scale": 1.0 / (0.5 * U0 * U0 * A0),
                "addToAdjoint": True,
            }
        },
    },
    "adjStateOrdering": "cell",
    "adjEqnOption": {"gmresRelTol": 1.0e-14, "pcFillLevel": 1, "jacMatReOrdering": "natural"},
    "normalizeStates": {
        "U": U0,
        "p": U0 * U0 / 2.0,
        "nuTilda": nuTilda0 * 10.0,
        "phi": 1.0,
    },
    "adjPartDerivFDStep": {"State": 1e-6, "FFD": 1e-3},
    "designVar": {},
}

# mesh warping parameters, users need to manually specify the symmetry plane and their normals
meshOptions = {
    "gridFile": os.getcwd(),
    "fileType": "openfoam",
    # point and normal for the symmetry plane
    "symmetryPlanes": [[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], [[0.0, 0.0, 0.1], [0.0, 0.0, 1.0]]],
}

# =============================================================================
# Design variable setup
# =============================================================================
DVGeo = DVGeometry("./FFD/wingFFD.xyz")
nTwists = DVGeo.addRefAxis("bodyAxis", xFraction=0.25, alignIndex="k")


def alpha(val, geo):
    aoa = val[0] * np.pi / 180.0
    inletU = [float(U0 * np.cos(aoa)), float(U0 * np.sin(aoa)), 0]
    DASolver.setOption("primalBC", {"U0": {"variable": "U", "patches": ["inout"], "value": inletU}})
    DASolver.updateDAOption()


def pitch(val, geo):
    for i in range(nTwists):
        geo.rot_z["bodyAxis"].coef[i] = -val[0]


def actuator(val, geo):
    actX = float(val[0])
    actY = float(val[1])
    actZ = float(val[2])
    actR1 = float(val[3])
    actR2 = float(val[4])
    actScale = float(val[5])
    actPOD = float(val[6])
    actExpM = float(val[7])
    actExpN = float(val[8])
    DASolver.setOption(
        "fvSource",
        {
            "disk1": {
                "type": "actuatorDisk",
                "source": "cylinderAnnulusSmooth",
                "center": [actX, actY, actZ],
                "direction": [1.0, 0.0, 0.0],
                "innerRadius": actR1,
                "outerRadius": actR2,
                "rotDir": "right",
                "scale": actScale,
                "POD": actPOD,
                "eps": 0.1,  # eps should be of cell size
                "expM": actExpM,
                "expN": actExpN,
                "adjustThrust": 1,
                "targetThrust": 0.2,
            },
        },
    )
    DASolver.updateDAOption()


# select points
iVol = 0
pts = DVGeo.getLocalIndex(iVol)
indexList = pts[:, :, :].flatten()
PS = geo_utils.PointSelect("list", indexList)
# shape
DVGeo.addGeoDVLocal("shape", lower=-1.0, upper=1.0, axis="y", scale=1.0, pointSelect=PS)
daOptions["designVar"]["shape"] = {"designVarType": "FFD"}
# pitch
DVGeo.addGeoDVGlobal("pitch", np.zeros(1), pitch, lower=-10.0, upper=10.0, scale=1.0)
daOptions["designVar"]["pitch"] = {"designVarType": "FFD"}
# AOA
DVGeo.addGeoDVGlobal("alpha", value=[alpha0], func=alpha, lower=0.0, upper=10.0, scale=1.0)
daOptions["designVar"]["alpha"] = {"designVarType": "AOA", "patches": ["inout"], "flowAxis": "x", "normalAxis": "y"}
# Actuator
DVGeo.addGeoDVGlobal(
    "actuator",
    value=[-0.5, 0.0, 0.05, 0.01, 0.4, 10.0, 0.8, 1.0, 0.5],
    func=actuator,
    lower=-100.0,
    upper=100.0,
    scale=1.0,
)
daOptions["designVar"]["actuator"] = {"actuatorName": "disk1", "designVarType": "ACTD"}
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

    DASolver()

elif args.task == "testAPI":

    DASolver.setOption("primalMinResTol", 1e-2)
    DASolver.updateDAOption()
    optFuncs.runPrimal()

else:
    print("task arg not found!")
    exit(0)
