
import sys, matplotlib.pyplot as plt
from pyFSI.post.anim import *
from pyFSI.execution.case import *
from pyFSI.post.prints import *

caseName = "Tosi"

# Create a log file
#sys.stdout = io.Logger()

# Create the case object
case = MFSICase(caseName).solve()

# Run the plots file
os.system("python3 plots.py " + caseName)



