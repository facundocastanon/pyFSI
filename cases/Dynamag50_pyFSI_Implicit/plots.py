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



# # Calculate the damping fit
# cutoff = 14500
# time = loadFiles([file('execution/TIME_time.out')])[0][0:cutoff]
# signal = loadFiles([file('solid/displacements.out')])[0][0:cutoff,0]
# dampingFit(time, signal, 145, 0)


fig, axs = plt.subplots(2, 2)
# Tip displacement evolution
p1 = plotFromFile(file('execution/TIME_time.out'),
				file('solid/SOLID_displacements.out'),
				xIndexes=np.s_[0:ls],
				yIndexes=np.s_[0:ls, 0],
				fig=fig,
				axe=axs[0,0],
				scale=1000)
p1.set(title="Beam root displacement", xLabel="t (s)", yLabel="U (mm)")



# Flow rate evolution
# p2 = plotFromFile(file('execution/TIME_time.out'),
# 				file('flow/flowRates.out'),
# 				scale=840,
# 				xIndexes=np.s_[0:ls],
# 				yIndexes=np.s_[0:ls, 0:2],
# 				fig=fig,
# 				axe=axs[0,1])
# p2.set(title="Evolution", xLabel="Time (s)", yLabel="Flow Rate (L/min)")

# Flow rate evolution
p2 = plotFromFile(file('execution/TIME_time.out'),
				file('flow/FLOW_flowRates.out'),
				scale=840,
				xIndexes=np.s_[0:ls],
				yIndexes=np.s_[0:ls, 0:2],
				fig=fig,
				axe=axs[0,1])
p2.set(title="Inlet flow rate", xLabel="t (s)", yLabel="Q0 (L/min)")
# p2.axe.set_title("Implicit")
# p2.axe.set_xlabel("Time(s)")
# p2.axe.set_ylabel("Flow Rate (m3/s)")



# Pressure center of top beam boundary
p3 = plotFromFile(file('execution/TIME_time.out'),
				file('flow/BEAMTOP_pressureCenter.out'),
				xIndexes=np.s_[0:ls],
				yIndexes=np.s_[0:ls],
				fig=fig,
				axe=axs[1,0],
				scale=1000)
p3.set(title="Pressure center", xLabel="t (s)", yLabel="u (mm)")

# Pressure distribution animation
tfilePath = pathlib.Path(__file__).parent.absolute()/'execution'/'TIME_time.out'
xfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'SOLID_mesh.out'
yfilePath = pathlib.Path(__file__).parent.absolute()/'flow'/'BEAMTOP_pressure.out'
p4 = pShape(xfilePath, yfilePath, time=-100, fig=fig, axe=axs[1,1])
p4.set(title="BEAMTOP Pressure Distribution", xLabel="x (m)", yLabel="p (Pa)", xLim=[0, 0.05], yLim=[-3E4, 3E4])

# anim = pAnimation(p4, tfilePath, name='response', startFrame=30000, endFrame=31200, frames=1200)
#anim.save(plotName, fps=30)

# p4.axe.set_title("Tosi Harvester")
# p4.axe.set_xlabel("x (m)")
# p4.axe.set_ylabel("Displacement (m)")
# p4.axe.set_xlim([0, 0.047])
# p4.axe.set_ylim([-0.01, 0.01])
# anim = pAnimation(p4, tfilePath, name='response', interval=1000)





# # Tip displacement
# p3 = plotFromFile(file('execution/TIME_time.out'),
# 				file('solid/SOLID_displacements.out'),
# 				xIndexes=np.s_[0:ls],
# 				yIndexes=np.s_[0:ls, -1],
# 				fig=fig,
# 				axe=axs[1,0],
# 				scale=1000)
# p3.set(title="Beam tip displacement", xLabel="t (s)", yLabel="u (mm)")

# # Velocities
# p4 = plotFromFile(file('execution/TIME_time.out'),
# 				file('solid/CIRCUIT_power.out'),
# 				xIndexes=np.s_[0:ls],
# 				yIndexes=np.s_[0:ls],
# 				fig=fig,
# 				axe=axs[1,1])
# p4.set(title="Damper dissipated power", xLabel="Time (s)", yLabel="Power (W)")




# Phase portrait
# p4= plotFromFile(file('solid/displacements.out'),
# 				file('solid/velocities.out'),
# 				xIndexes=np.s_[0:ls, -1],
# 				yIndexes=np.s_[0:ls, -1],
# 				fig=fig,
# 				axe=axs[1,1])
# p4.set(title="Phase Portrait", xLabel="Displacement (m)", yLabel="Velocity (m/s)")


# p4= pModes(file('solid/mesh.out'),
# 		   file('solid/eigenValues.out'),
#      		file('solid/eigenVectors.out'),
# 				fig=fig,
# 				axe=axs[1,1])
# p4.set(title="Phase Portrait", xLabel="Displacement (m)", yLabel="Velocity (m/s)")

# Beam deformed shape animation
# tfilePath = pathlib.Path(__file__).parent.absolute()/'execution'/'time.out'
# xfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'mesh.out'
# yfilePath = pathlib.Path(__file__).parent.absolute()/'solid'/'displacements.out'
# p4 = pShape(xfilePath, yfilePath, time=-100, fig=fig, axe=axs[1,1])
# p4.set(title="Deformed shape", xLabel="x (m)", yLabel="y (m)", xLim=[0, 0.05], yLim=[-0.002, 0.002])

# anim = pAnimation(p4, tfilePath, name='response', startFrame=30000, endFrame=31200, frames=1200)
#anim.save(plotName, fps=30)

# p4.axe.set_title("Tosi Harvester")
# p4.axe.set_xlabel("x (m)")
# p4.axe.set_ylabel("Displacement (m)")
# p4.axe.set_xlim([0, 0.047])
# p4.axe.set_ylim([-0.01, 0.01])
# anim = pAnimation(p4, tfilePath, name='response', interval=1000)
# anim.save('magneticTransientLFBeam', fps=30)

# plt.savefig("DynamagSN50_SweepUp.pdf")

with open(plotName + '.obj', 'wb') as file:
	pickle.dump(fig, file)

plt.show()
plt.close(fig)
