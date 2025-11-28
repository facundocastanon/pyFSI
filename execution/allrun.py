#!/usr/bin/env python3
from pyFSI.execution.case import *
from pyFSI.post.prints import *
import sys
import os

def helpMessage():
    print ("")
    print ("Run a pyFSI case. You must be in a pyFSI case directory in order to run this command,")
    print ("since it reads the .json file in the root.")
    print ("")
    print ("Usage: fsiRun [OPTIONS]")
    print ("")
    print ("Options:")
    print ("    -p: the plots are opened after the case is done running.")
    print ("    -f: a .log file is created after the case is done running.")
    print ("    -c: console output enabled.")
    print ("    -l: both a .log file is created as well as console output is enabled.")
    print ("    -pf: for pyFSI-preCICE cases only. Run the /flow/ participant in the root directory.")
    print ("    -ps: for pyFSI-preCICE cases only. Run the /solid/ participant in the root directory.")
    print ("    -h: display the help message.")
    print ("")
    print ("For more information on pyFSI core commands visit the documentation page")
    print ("https://gitlab.com/martinsaravia/pyfsi")

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
# Check the arguments
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
validArguments = ['-p', '-f', '-c', '-l', '-h', '-pf', '-ps']

if validArguments[4] in sys.argv:
    helpMessage()
    sys.exit(1)
elif any(arg not in validArguments for arg in sys.argv[1:]):
    print("Invalid argument.")
    helpMessage()
    sys.exit(1)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
# For pyFSI-preCICE cases only
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
if validArguments[5] in sys.argv:
    print("Running the fluid participant...")
    os.chdir('flow')
elif validArguments[6] in sys.argv:
    print("Running the solid participant...")
    os.chdir('solid')

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
# Define case name
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
# Files in the current directory
files = os.listdir(os.getcwd())

# Empty variable for the case name
caseName = ""

for file in files:
    # Find the file that ends in .json
    if file.endswith(".json"):
        # Extract the file name without the extension
        caseName = os.path.splitext(file)[0]
        print("Running case", caseName, "...")
        break

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
# Create a log file? The user decides
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
if validArguments[1] in sys.argv:
    print("  Creating log file...")
    sys.stdout = io.Logger(console=False, file=True)
elif validArguments[2] in sys.argv:
    print("  Creating console log...")
    sys.stdout = io.Logger(console=True, file=False)
elif validArguments[3] in sys.argv:
    print("  Creating both console log and log file...")
    sys.stdout = io.Logger(console=True, file=True)
else:
    print("  Neither creating log files nor printing on console...")
    sys.stdout = io.Logger(console=False, file=False)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
# Create the case object
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
case = MFSICase(caseName)

# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
# Solve the case
# ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
case.solve()
