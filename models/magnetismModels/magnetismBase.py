"""@package docstring
Base class for the fluid models

Only one boundary can be flexible
"""

from abc import ABCMeta, abstractmethod
from pyFSI.models.properties.materialProperties import fluids as db
from pyFSI.models.properties.dimensionlessNumbers import *
from pyFSI.mesh.region.fsiRegion1D import fsiRegion1D


class magnetismModel(metaclass=ABCMeta):
    """ Base class for the magnetic models"""
    def __repr__(self):
        return 'magneticModel Abstract Class'

    def __init__(self, execution, control, mesh, time):
        """Constructor"""

        # ----- Public attribues ----- #
        self.base = "magnetism"
        self.control = control
        self.name = control['name']
        self.output = []
        self.path = execution['paths']['magnetismPath']  # Associated path
        self.dimNumbers = {}  # Dimensionless numbers
        self.varMap = {
            "numbers":      "dimNumbers"
        }

        # ----- Private attributes ----- #
        self._execution = execution
        self._mesh = mesh
        self._updated = False
        self._fluid = None
        self._time = time
        if execution['debug'] == 'yes':
            self._debug = True
        else:
            self._debug = False

        # ----- Procedures ----- #

    # Calculate the dimensionless numbers
    def calcNumbers(self):
        pass

    # Pure virtual methods

    # Getters
    def execution(self):
        return self._execution

    def control(self):
        return self.control

    def mesh(self):
        return self._mesh

    def time(self):
        return self._time

    def closeOutput(self):
        for i in self.output:
            i.close()  # Close all files

    def finish(self):
        self.closeOutput()

    def updated(self):
        return self._updated






