from pyFSI.models.properties.dimensionlessNumbers import *
from abc import ABC, abstractmethod

class fsiBase(ABC):
    def __init__(self, execution, control, TIME, OREG):
        
        # ----- Public attributes ----- #
        self.base = "fsi"
        self.control = control
        self.name = control['name']
        self.dof = None  # Number of degrees of freedom
        self.name = control['name']
        self.dimNumbers = {}  # Dimensionless numbers
        
        # Output control
        self.path = execution['paths']['fsiPath']  # Associated path
        self.output = []
        self.varMap = {
            "numbers": "dimNumbers"
        }

        # ----- Private attributes ----- #
        if 'participants' in control:
            for name in control['participants']:
                obj = OREG.get(name)
                if obj.base == 'solid': 
                    self._solid = obj
                elif obj.base == 'flow':
                    self._flow = obj
                elif obj.base == 'magnetism':
                    self._magnetism = obj
        self._execution = execution
        self._time = TIME
        self._debug = (execution.get('debug', False) == 'yes')
        


    def calcNumbers(self):
        self.dimNumbers["Cy"] = CauchyNumber(self)
        self.dimNumbers["Ms"] = massNumber(self)
        self.dimNumbers["Vr"] = reducedVelocityNumber(self)

    # Getters
    def execution(self):
        return self._execution

    def solid(self):
        return self._solid

    def flow(self):
        return self._flow

    def magnetism(self):
        return self._magnetism

    def time(self):
        return self._time

    def numbers(self):
        return self.dimNumbers

