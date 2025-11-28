# --------------------------------------------------------------------------- #
#    p    #     version: 0.2.0
#    y    #     date: 10/4/2021
#    F    #     author: Martin Saravia
#    S    #     description: routines for case management
#    I    #     return:
# --------------------------------------------------------------------------- #
# Notes:
# Available functions are readCase, cleanCase and runCase
# --------------------------------------------------------------------------- #
import importlib
import json
import pathlib
import shutil
import sys
import os
from colorama import Fore, Back, Style
from copy import deepcopy
from pyFSI.execution.utilities import banner, section
from pyFSI.execution import io, time
from pyFSI.mesh.fsiMesh1D import fsiMesh1D


class MFSICase:
    def __init__(self, caseName):
        banner("pyFSI v0.4 - Martin Saravia 2023")        # Timing starts
        # Public attributes
        self.DICT = None  # The input dictionary
        self.IODB = None  # The Input-Output Data Base
        self.OREG = None  # The object registry
        self.SOLU = None  # The solver object (yes, weird name)
        self.MFSI = None  # The Magnetic-Fluid-Structure-Interaction object
        self.TIME = None  # The time database

        # Get the path to the caller location
        casePath = pathlib.Path(os.path.abspath(caseName).rsplit('/', 1)[0])  # split cut the name of the file from the Path (added by os)
        print("Running pyFSI case from: ", casePath)

        # Clean the directory tree if present.
        self.clean(casePath)

        # Read the case file and create the fsi object list from the original dictionary (input file)
        self._readInput(casePath, caseName)

        # Create all objects
        self._createObjects()

    def solve(self):
        # Create the case list (contains solver objects, which have a solution)
        section("Calling the solver...")

        self.SOLU.solve()
    
        # End run
        print(f"Problem solved in {self.TIME.elapsed()} seconds ")
        print("************************************************************************")
        

    # Read the case file
    def _readInput(self, casePath, caseName):
        # casePath = pathlib.Path(__file__).parent.absolute() / ".." / "cases" / caseName
        fileName = caseName + ".json"
        caseFile = casePath / fileName
    
        with open(caseFile, 'r') as f:
            self.DICT = json.load(f)  # Load the input file
    
        # Add case paths to the execution dictionary
        paths = {}
        paths['caseName'] = caseName
        paths['casePath'] = casePath
        
        # Copy or create the output dir structure
        paths['executionPath'] = casePath / "execution"
        if os.path.isdir(paths['casePath'] / "execution.0"):
            shutil.copytree(paths['casePath'] / "execution.0", paths['executionPath'])
        else:
            paths['executionPath'].mkdir(parents=True, exist_ok=True)
        if "solid" in self.DICT: 
            paths['solidPath'] = casePath / "solid"
            if os.path.isdir(paths['casePath'] / "solid.0"):
                shutil.copytree(paths['casePath'] / "solid.0", paths['solidPath'])
            else:
                paths['solidPath'].mkdir(parents=True, exist_ok=True)
        if "boundary" in self.DICT: 
            paths['boundaryPath'] = casePath / "boundary"
            if os.path.isdir(paths['casePath'] / "boundary.0"):
                shutil.copytree(paths['casePath'] / "boundary.0", paths['boundaryPath'])
            else:
                paths['boundaryPath'].mkdir(parents=True, exist_ok=True)
        if "flow" in self.DICT:
            paths['flowPath'] = casePath / "flow"
            if os.path.isdir(paths['casePath'] / "flow.0"):
                shutil.copytree(paths['casePath'] / "flow.0", paths['flowPath'])
            else:
                paths['flowPath'].mkdir(parents=True, exist_ok=True)
        if "fsi" in self.DICT:
            paths['fsiPath'] = casePath / "fsi"
            if os.path.isdir(paths['casePath'] / "fsi.0"):
                shutil.copytree(paths['casePath'] / "fsi.0", paths['fsiPath'])
            else:
                paths['fsiPath'].mkdir(parents=True, exist_ok=True)
        if "magnetism" in self.DICT:
            paths['magnetismPath'] = casePath / "magnetism"
            if os.path.isdir(paths['casePath'] / "magnetism.0"):
                shutil.copytree(paths['casePath'] / "magnetism.0", paths['magnetismPath'])
            else:
                paths['magnetismPath'].mkdir(parents=True, exist_ok=True)

        self.DICT['execution']['paths'] = paths

    def clean(self, casePath, boundary=True, solid=True, flow=True, fsi=True, magnetism=True, execution=True):
        # Get the path to the caller location
        print("Cleaning the case at:", casePath)
        if execution:
            shutil.rmtree(casePath / "execution", ignore_errors=True)
        if boundary:
            shutil.rmtree(casePath / "boundary", ignore_errors=True)            
        if solid:
            shutil.rmtree(casePath / "solid", ignore_errors=True)
        if flow:
            shutil.rmtree(casePath / "flow", ignore_errors=True)
        if fsi:
            shutil.rmtree(casePath / "fsi", ignore_errors=True)
        if magnetism:
            shutil.rmtree(casePath / "magnetism", ignore_errors=True)

    def _createObjects(self):
        print("-> Preprocessing...")

        # Create the object registry
        self.OREG = io.ObjectRegistry()

        self.TIME = time.time(self.DICT["execution"])
        self.OREG.append(self.TIME)

        # Create the mesh object
        mesh = getattr(fsiMesh1D, self.DICT["mesh"]["type"])(self.DICT["mesh"])
        self.OREG.append(mesh)
        mesh.addObjectRegistry(self.OREG) # Add the object registy to the mesh
        
            
        # Create the boundaries
        if "boundary" in self.DICT:
            boundary = {}
            print("---> Creating the boundaries...")
            # Create the dict of boundary objects
            for b in self.DICT["boundary"]:
                # Import the module (file .py)
                module = importlib.import_module('pyFSI.mesh.region.' + b['type'])                
                
                if "method" in b:
                    # get the class name inside the module
                    obj = getattr(module, b['type'])(self.DICT['execution'], b, mesh)
                    # Create the object and store in in the boundary dictionary
                else:
                    # get the class name inside the module
                    obj = getattr(module, b['type'])(self.DICT['execution'], b, solid)
    
                self.OREG.append(obj)
    
                boundary[b['name']] = obj
    
        # Create the magnetism
        if "magnetism" in self.DICT:
            print("---> Creating the magnetism...")
            magnetism = {}
            for m in self.DICT["magnetism"]:
                # Load the module
                moduleName = 'pyFSI.models.magnetismModels.'
                className =  m['formulation']
                magnetismModule = importlib.import_module(moduleName 
                                                         + className 
                                                         + 'Model')
                # Create the object
                obj = getattr(magnetismModule, className)(self.DICT['execution'], m, mesh, self.TIME)
                magnetism[m['name']] = obj
                self.OREG.append(obj)
        else:
            magnetism = [] # Remove this line after removing the fsi dependence
        
        # Create the solids
        if "solid" in self.DICT:
            print("---> Creating the solids...")
            for solidDict in self.DICT['solid']:
                # Create the solid model
                solidModule = importlib.import_module('pyFSI.models.solidModels.' +
                                                    solidDict['formulation'] +
                                                    'Model')
                solid = getattr(solidModule,
                                solidDict['formulation'])(self.DICT['execution'],
                                                        solidDict,
                                                        mesh,
                                                        boundary,
                                                        self.TIME)
                self.OREG.append(solid)

    
        # Create the flows
        if "flow" in self.DICT:
            print("---> Creating the flow...")
            # Create the flow object
            flowModule = importlib.import_module('pyFSI.models.flowModels.' +
                                                 self.DICT['flow']['formulation'] +
                                                 'Model')
            flow = getattr(flowModule,
                           self.DICT['flow']['formulation'])(self.DICT['execution'],
                                                             self.DICT['flow'],
                                                             mesh,
                                                             boundary,
                                                             self.TIME)
            self.OREG.append(flow)

    
        # Create a list of fsi objects
        if "fsi" in self.DICT:
            print("---> Building the FSI object...")
            fsiModule = importlib.import_module('pyFSI.models.fsiModels.' +
                                                self.DICT['fsi']['formulation'] +
                                                'Model')
            mfsi = getattr(fsiModule,
                           self.DICT['fsi']['formulation'])(self.DICT['execution'],
                                                            self.DICT['fsi'],
                                                            self.TIME, 
                                                            self.OREG)

            self.MFSI = mfsi
            self.OREG.append(mfsi)


        # Create the output data base
        self.IODB = io.IODataBase(self.DICT, self.OREG)
    
        # Create the solver object
        print("---> Building the solver object...")
        solverType = self.DICT['execution']['solver']['type']
        solverModule = importlib.import_module('pyFSI.solvers.' + solverType + 'Solver')
        if self.MFSI is not None:
            solver = getattr(solverModule, solverType)(self.MFSI, self.IODB)
        else: # Solid or Fluid solvers, no interaction
            control = self.DICT['execution']['solver']
            if solverType == "transientSolid":
                solver = getattr(solverModule, solverType)(control, solid, self.TIME, self.IODB)
            if solverType == "transientFlow":
                solver = getattr(solverModule, solverType)(control, flow, self.TIME, self.IODB)
        self.SOLU = solver
    
        print("Preprocessing finished successfully...", "\n")
