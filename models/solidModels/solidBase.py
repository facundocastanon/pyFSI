from abc import ABCMeta
from pyFSI.models.properties.materialProperties import solids as dataBase
from pyFSI.models.properties.dimensionlessNumbers import *
from pyFSI.execution.utilities import setParameter
from pyFSI.execution.utilities import getObjectbyName
from colorama import Fore, Back, Style

class solidModel(metaclass=ABCMeta):
    def __repr__(self):
        return 'solidModel Abstract Class'

    def __init__(self, execution, control, mesh, boundary, time):

        # ----- Public attribues ----- #
        self.base = "solid"
        self.control = control
        self.name = control['name']
        self.boundary = []
        self.dof = None
        self.lRef = None  # Reference length
        self.uRef = None  # Reference displacement
        self.vRef = None  # Reference velocity
        self.dimNumbers = {}  # Dimensionless numbers
        self.output = []  # List of file objects
        self.path = execution['paths']['solidPath']  # Associated path
        self.varMap = {}
        if "parameter" in control:
            self.varMap["parameter"] = "control" + control['parameter']['command']
        
        # ----- Private attributes ----- #'
        self._execution = execution
        self._mesh = mesh
        self._boundaryDict = boundary
        self._updated = False
        self._time = time
        self._debug = (execution.get('debug', False) == 'yes')
        
        # ----- Procedures ----- #
        # Get the material properties if present
        if 'material' in control:
            if 'db' in control['material']:  # Specify the material from the database
                self._material = dataBase[control['material']['db']]
            else:
                self._material = control['material']  # Specify the properties directly

        # Create the boundary objects
        if 'boundary' in control:
            self.makeBoundary(boundary)

            
    def parametrize(self, time):
        setParameter(self.control, time)
        self.__init__(self._execution, self.control, self._mesh, self.boundary, time)

    # Calculate the dimensionless numbers
    def calcNumbers(self):
        self.dimNumbers["Dp"] = displacementNumber(self)
        self.dimNumbers["Eg"] = elastoGravityNumber(self)
        
    def makeBoundary(self, boundary):
        for bname in self.control['boundary']:
            self.boundary.append(getObjectbyName(bname, boundary))
    


    # Getters
    def execution(self):
        return self._execution

    def mesh(self):
        return self._mesh
    
    def material(self):
        return self._material

    def closeOutput(self):
        for i in self.output:
            i.close()  # Close all files

    def finish(self):
        self.closeOutput()

    def updated(self):
        return self._updated

    def time(self):
        return self._time