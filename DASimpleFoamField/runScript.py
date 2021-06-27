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
nCells = 4032

adjJacOpt = "JacobianFree"
mode = "reverse"

if args.task == "runForwardAD":
    adjJacOpt = "JacobianFD"
    mode = "forward"

# Set the parameters for optimization
daOptions = {
    "designSurfaces": ["wing"],
    "solverName": "DASimpleFoam",
    "primalMinResTol": 1.0e-10,
    "adjJacobianOption": adjJacOpt,
    "useAD": {"mode": mode, "dvName": args.dvName, "seedIndex": args.seedIndex},
    "primalBC": {
        "U0": {"variable": "U", "patches": ["inout"], "value": [U0, 0.0, 0.0]},
        "p0": {"variable": "p", "patches": ["inout"], "value": [p0]},
        "nuTilda0": {"variable": "nuTilda", "patches": ["inout"], "value": [nuTilda0]},
        "useWallFunction": True,
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
    "adjEqnOption": {"gmresRelTol": 1.0e-10, "pcFillLevel": 1, "jacMatReOrdering": "rcm"},
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


def betaSA(val, geo):
    for idxI, v in enumerate(val):
        DASolver.setFieldValue4GlobalCellI(b"betaSA", v, idxI)
        DASolver.updateBoundaryConditions(b"betaSA", b"scalar")


def alphaPorosity(val, geo):
    for idxI, v in enumerate(val):
        DASolver.setFieldValue4GlobalCellI(b"alphaPorosity", v, idxI)
        DASolver.updateBoundaryConditions(b"alphaPorosity", b"scalar")


# select points
iVol = 0
pts = DVGeo.getLocalIndex(iVol)
indexList = pts[:, :, :].flatten()
PS = geo_utils.PointSelect("list", indexList)
# AOA
DVGeo.addGeoDVGlobal("alpha", value=[alpha0], func=alpha, lower=0.0, upper=10.0, scale=1.0)
daOptions["designVar"]["alpha"] = {"designVarType": "AOA", "patches": ["inout"], "flowAxis": "x", "normalAxis": "y"}
# Beta
beta0 = np.ones(nCells, dtype="d")
DVGeo.addGeoDVGlobal("beta", value=beta0, func=betaSA, lower=0.0, upper=10.0, scale=1.0)
daOptions["designVar"]["beta"] = {"designVarType": "Field", "fieldName": "betaSA", "fieldType": "scalar"}
# AlphaPorosity
alphaPorosity0 = np.zeros(nCells, dtype="d")
DVGeo.addGeoDVGlobal("alphaPorosity", value=alphaPorosity0, func=alphaPorosity, lower=0, upper=100.0, scale=1.0)
daOptions["designVar"]["alphaPorosity"] = {
    "designVarType": "Field",
    "fieldName": "alphaPorosity",
    "fieldType": "scalar",
}
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
