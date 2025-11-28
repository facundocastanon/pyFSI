"""@package docstring
Base class for the fluid models

Only one boundary can be flexible
"""

from abc import ABCMeta, abstractmethod
from pyFSI.models.properties.materialProperties import fluids as db
from pyFSI.models.properties.dimensionlessNumbers import *
from pyFSI.mesh.region.fsiRegion1D import fsiRegion1D
import matplotlib.pyplot as plt


class flowModel(metaclass=ABCMeta):
    """ Base class for the flow models"""
    def __repr__(self):
        return 'fluidModel Abstract Class'

    def __init__(self, execution, control, mesh, boundary, time):
        """Constructor"""

        # ----- Public attribues ----- #
        self.base = "flow"
        self.control = control
        self.name = control['name']
        self.dof = None  # Number of DOF of the model
        self.regions = []  # Regions comprising the domain
        self.vRef = None  # Reference velocity
        self.lRef = None  # Reference length
        self.dRef = None  # Reference inlet size
        self.output = []
        self.path = execution['paths']['flowPath']  # Associated path
        self.dimNumbers = {}  # Dimensionless numbers
        self.varMap = {
            "numbers":      "dimNumbers"
        }

        # ----- Private attributes ----- #
        self._execution = execution
        self._mesh = mesh
        self._boundary = boundary
        self._updated = False
        self._fluid = None
        self._time = time
        self._debug = (execution.get('debug', False) == 'yes')


        # ----- Procedures ----- #
        # Get the fluid properties
        if "db" in control["fluid"]:  # Specify the fluid from the database
            self._fluid = db[control["fluid"]["db"]]
        else:
            self._fluid = control["fluid"]  # Specify the properties directly

        # Create the region objects
        [self.regions.append(fsiRegion1D(i, mesh, boundary)) for i in control['regions']]

        # Plot the boundaries
        if 'plot' in control:
            if control['plot'] == 'yes':
                self.plot()

    # Calculate the dimensionless numbers
    def calcNumbers(self):
        """ Calculates dimensionless numbers"""
        self.dimNumbers["Re"] = ReynoldsNumber(self)
        self.dimNumbers["Rd"] = ReynoldsNumber(self, type="Rd")
        self.dimNumbers["Fr"] = FroudeNumber(self)

    # Plot function
    def plot(self):
        for region in self.regions:
            for boundary in region.boundaries():
                plt.plot(boundary.x, boundary.y)
        plt.title('Flow Domain')
        plt.xlabel("x")
        plt.ylabel("y")
        # Setting major and minor grid lines for both x and y axes
        plt.grid(which='both', linestyle='-', linewidth=0.5)
        plt.minorticks_on()
        plt.grid(which='minor', linestyle=':', linewidth=0.5)
        plt.show()



    # Getters
    def execution(self):
        return self._execution

    def mesh(self):
        return self._mesh

    def boundary(self):
        return self._boundary

    def fluid(self):
        return self._fluid

    def closeOutput(self):
        for i in self.output:
            i.close()  # Close all files

    def finish(self):
        self.closeOutput()

    def updated(self):
        return self._updated

    def time(self):
        return self._time
