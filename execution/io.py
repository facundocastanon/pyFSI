""" Input-Output Module
This module holds the classes used for writing input and output data.
"""
import filecmp
import sys
from abc import ABC, abstractmethod
import numpy as np
from pyFSI.execution.errors import error
from pyFSI.execution.utilities import writeVersionHeader

class IOFile:
    """
    Class for output files. It takes and object and a variable and creates
     an output file that is written when the database
    write method is called.
    """
    def __init__(self, obj, variable,  mode='a+', bufferSize=8192):
        filename = variable + ".out"
        self.location = obj.path / filename
        self.obj = obj  # Object from which we read the data
        
        # Flag to update only if it is not an initial variable. Not flagged 
        # variables do not go to the ODB, they are flushed and closed.
        self.update = True
        if "initial" in obj.varMap:
            if variable in obj.varMap["initial"]:
                self.update = False
            
        # Get the variable
        try:
            self.var = obj.varMap[variable]
        except KeyError:
            error("ERROR: Variable: " + variable + " is not available for output for object named: " + obj.name + "\n" +
                  ". Valid variables are: " + str(list(obj.varMap.keys()))[:])

        try:
            self.expr = "self.obj." + self.var  # Expression to execute to evaluate the variable value
        except AttributeError:
            msg = ["ERROR: The variable " + obj.varMap[variable] + 
                   " is not a valid output variable for object " + self.obj]
            error(msg)

        self.file = open(self.location,
                         mode,
                         buffering=bufferSize)

        writeVersionHeader(self.file)

        # Choose the writer
        value = eval(self.expr)
        if (isinstance(value, np.floating) or isinstance(value, float) or 
            isinstance(value, int) or isinstance(value, np.float64)):
            self.writer = NumWriter(obj, self.file)
        elif isinstance(value, np.ndarray):
            if isinstance (value[0], np.ndarray):
                self.writer = ArrayOfArraysWriter(obj, self.file)
            elif value.ndim == 1:
                self.writer = ArrayWriter1D(obj, self.file)
            elif value.ndim == 2:
                self.writer = ArrayWriter2D(obj, self.file)
            else:
                msg = "---> ERROR: pyFSI cannot find a suitable writer for file"
                error(msg)
        elif isinstance(value, dict):  # If list, split it)
            self.writer = DictWriter(obj, self.file, eval(self.expr))
        elif isinstance(value, list):
            try: # Test if it is a list of np arrays
                if isinstance (value[0], np.ndarray):
                    self.writer = ListOfArraysWriter(obj, self.file)
            except: # if not a list of arrays then is a list of numbers
                self.writer = ListWriter(obj, self.file)
        else:
            msg = "---> ERROR: pyFSI cannot find a suitable writer for file"
            error(msg)
            
            

    def write(self):
        value = eval(self.expr)
        self.writer.line(self.file, value)

    def close(self):
        self.file.close()

class Writer(ABC):
    """ Base class for writing data to files"""
    def __init__(self, obj, file):
        self.header(obj, file)

    def header(self, obj, file, data=None):
        try:
            name = obj.name
        except:
            name = "None"
        file.write('#   Results file for object: ' + name + '\n')

    @abstractmethod
    def line(self, file, value):
        pass
    
class ListOfArraysWriter(Writer):
    def __init__(self, obj, file):
        super().__init__(obj, file)
    def line(self, file, value):
        for i in value:
            file.write(" ".join(map(str, i[:])) + '\n')
            
class ArrayOfArraysWriter(Writer):
    def __init__(self, obj, file):
        super().__init__(obj, file)
    def line(self, file, value):
        for i in value:
            file.write(" ".join(map(str, i[:])) + '\n') 

class ListWriter(Writer):
    def __init__(self, obj, file):
        super().__init__(obj, file)

    def line(self, file, value):
        file.write(" ".join(map(str, value)) + '\n')  

class ArrayWriter1D(Writer):
    def __init__(self, obj, file):
        super().__init__(obj, file)

    def line(self, file, value):
        file.write(" ".join(map(str, value[:])) + '\n') 
            
class ArrayWriter2D(Writer):
    def __init__(self, obj, file):
        super().__init__(obj, file)

    def line(self, file, value):
        for i in value:
            file.write(" ".join(map(str, i[:])) + '\n')  
        
class NumWriter(Writer):
    def __init__(self, obj, file):
        super().__init__(obj, file)

    def line(self, file, value):
        file.write(str(value) + '\n')

class DictWriter:
    def __init__(self, obj, file, dict):
        self.header(obj, file, dict)

    def header(self, obj, file, dict):
        file.write('# ')
        for k, v in dict.items():
            file.write(k + 17 * ' ')
        file.write('\n')

    def line(self, file, dict):
        for k, v in dict.items():
            file.write(str(v) + 4 * " ")
        file.write('\n')


class IODataBase:
    """
    Input-Output data base. It holds all the output files (of type IOFile).
    Note: f the file is not flagged to be updated, write and close without 
          adding it to the ODB
    """
    def __init__(self, caseDict, registry):
        self.files = []
        # Add current time time to the database
        time = registry.get('TIME')
        self.files.append(IOFile(time, 'time'))
        # Add the rest of the variables to the database
        for key, val in caseDict['execution']['output'].items():
            obj = registry.get(key)  # Find the object by name
            for i in val:  # Iterate through variable names of the current object
                file = IOFile(obj, i)
                if file.update: 
                    self.files.append(file)
                else:
                    file.write()
                    file.close()
                
            
        
    def write(self):
        for file in self.files:
            file.write()

    def close(self):
        for file in self.files:
            print("---> Results written to file: " + str(file.location))
            file.close()

class ObjectRegistry:
    """
    A class for registering all objects. It has a get method for retrieving an object by name.
    """
    def __init__(self):
        self.objects = []
        self.names = []

    def append(self, obj):   # Add an object to the registry
        self.objects.append(obj)
        self.names.append(obj.name)

    def get(self, name):  # Get and object from the registry
        for i in self.objects:
            if i.name == name:
                return i
        # Raise error if object is not found
        err_msg = ("Object named " + name + " not found. Available objects: " 
                   + str(self.names))
        raise RuntimeError(err_msg)


class Logger:
    """
    Class for logging simultaneously to the console and to a file.
    From https://stackoverflow.com/questions/616645/how-to-duplicate-sys-stdout-to-a-log-file
    """
    def __init__(self, name='stdout.log', mode='w', file=False, console=True):
        if file:
            self.file = open(name, mode)
        else:
            self.file = file
        self.stdout = sys.stdout
        sys.stdout = self
        
        self.console = console

    def __del__(self):
        sys.stdout = self.stdout
        if self.file:
            self.file.close()

    def write(self, data):
        if self.file:
            self.file.write(data)
        if self.console:
            self.stdout.write(data)

    def flush(self):
        if self.file:
            self.file.flush()
