import pathlib
import pickle
from pyFSI.post.anim import *
from pyFSI.util.dampingFitting import *
import sys


plotName = "dyna50"

try:
	replot = sys.argv[1]
except:
    replot = 'no'
    
if replot == 'replot':
	print('---> Loading existing plot object named ' + plotName + '.obj')
	with open(plotName + '.obj', 'rb') as file:
		pickle.load(file)
	plt.show()
	sys.exit()
else:
	print('---> Creating/overwriting plot object named ' + plotName + '.obj')
	 

# Parse input commands
if len(sys.argv) > 1:
	ls = int(sys.argv[1]) - 5
else:
	ls = -10

def file(name):
	filePath = pathlib.Path(__file__).parent.absolute()/name
	return filePath


fig, axs = plt.subplots(2, 2)

# Tip displacement evolution
p1 = plotFromFile(file('execution/TIME_time.out'),
				file('solid/SOLID_displacements.out'),
				xIndexes=np.s_[0:ls],
				yIndexes=np.s_[0:ls, -1],
				fig=fig,
				axe=axs[0,0],
				scale=1000)
p1.set(title="Beam tip displacement", xLabel="t (s)", yLabel="U (mm)")

# Flow rate evolution
p2 = plotFromFile(file('execution/TIME_time.out'),
				file('flow/FLOW_flowRates.out'),
				scale=840,
				xIndexes=np.s_[0:ls],
				yIndexes=np.s_[0:ls, 0:2],
				fig=fig,
				axe=axs[0,1])
p2.set(title="Inlet flow rate", xLabel="t (s)", yLabel="Q0 (L/min)")

# Pressure center of top beam boundary
p3 = plotFromFile(file('execution/TIME_time.out'),
				file('flow/BEAMTOP_pressureCenter.out'),
				xIndexes=np.s_[0:ls],
				yIndexes=np.s_[0:ls],
				fig=fig,
				axe=axs[1,0],
				scale=1000)
p3.set(title="Pressure center", xLabel="t (s)", yLabel="c (mm)")

# Pressure distribution animation
tfilePath = pathlib.Path(__file__).parent.absolute()/'execution'/'TIME_time.out'
xfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'SOLID_mesh.out'
yfilePath = pathlib.Path(__file__).parent.absolute()/'flow'/'BEAMTOP_pressure.out'
p4 = pShape(xfilePath, yfilePath, time=-100, fig=fig, axe=axs[1,1], scale=1000)
p4.set(title="BEAMTOP Pressure Distribution at ts=-100", xLabel="x (mm)", yLabel="p (Pa)")


plt.show()
