# --------------------------------------------------------------------------- #
#    p    #     version: 0.1
#    y    #     date: 07/12/2020
#    F    #     author: Martin Saravia
#    S    #     description: Base class for all solvers
#    I    #
# --------------------------------------------------------------------------- #

import importlib
from abc import ABCMeta, abstractmethod
import numpy as np
#from tqdm import tqdm
from pyFSI.mesh.fsiMesh1D import fsiMesh1D


class solverBase:
    def __init__(self, fsi, odb):
        # Private Attributes
        self._fsi = fsi
        self._odb = odb
        self._time = fsi.time()
        self._execution = fsi.execution()

        # Public Attributes
        self.control = fsi.execution()['solver']

        # # Output
        # bufferSize = 1
        # self.output = []
        # self.output.append(open(self._execution['paths']['fsiPath'] / 'time.out', 'a+', buffering=bufferSize))
        
    def residual(self, x0, x1):
        dx = x0 - x1
        res = np.linalg.norm(dx, 2) / np.linalg.norm(x0, 2)
        return res
        
    def createSolidIntegrator(self):
        # Choose the solid solver
        solidIntegratorType = self.control['integrators']['solid']['type']
        if solidIntegratorType in ["Radau", "RK45", "RK23", "DOP853", "BDF", "LSODA"]:
            solidIntegratorModule = importlib.import_module('pyFSI.math.ode.' + "scipyIVPSolid")
            self.solidIntegrator = getattr(solidIntegratorModule, "scipyIVPSolid")(self._fsi.solid(), self._time)
        elif solidIntegratorType in ["Newmark", "generalizedAlpha"]:
            solidIntegratorModule = importlib.import_module('pyFSI.math.ode.' + solidIntegratorType)
            self.solidIntegrator = getattr(solidIntegratorModule, solidIntegratorType)(self._fsi.solid(), self._time)

    def createFlowIntegrator(self):
        # Choose the flow Integrator
        flowIntegratorType = self.control['integrators']['flow']['type']
        if flowIntegratorType in ["Radau", "RK45", "RK23", "DOP853", "BDF", "LSODA"]:
            flowIntegratorModule = importlib.import_module('pyFSI.math.ode.' + "scipyIVPFlow")
            self.flowIntegrator = getattr(flowIntegratorModule, "scipyIVPFlow")(self._fsi.flow(), self._time)



    # Abstract methods
    @abstractmethod
    def solve(self):
        pass

    # Getters
    def execution(self):
        return self._execution

    def odb(self):
        return self._odb
    
    def fsi(self):
        return self._fsi

