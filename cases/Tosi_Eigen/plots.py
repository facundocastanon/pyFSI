import pathlib, sys, os
from pyFSI.post.anim import *


ss = 0

plotName = "Tosi_Eigen"

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
	print('---> Existing plot object named ' + plotName + '.obj' +  ' was not found. Building it.')
 
 
# Parse input commands
if len(sys.argv) > 1:
	ls = -10
else:
	ls = -10

def file(name):
	filePath = pathlib.Path(__file__).parent.absolute()/name
	return filePath

fig, axs = plt.subplots(1, 1)

p = plotFromFile(file('fsi/eigenValues.out'), 
				file('fsi/eigenValues.out'),
				color = file('flow/flowRates.out'),
				xIndexes=np.s_[:, 0:10], 
				yIndexes=np.s_[:, 0:10],
				fig=fig, 
				axe=axs, 
				type="Eigenvalues",
				scatter=True)
p.set(xLabel="Real Part", yLabel="Imaginary Part", title="Eingenvalue Progression")
p.axe.set_xlim([-2000, 2000])
p.axe.set_ylim([-1000, 60000])
#plt.savefig('Eigenvalues' + sys.argv[1] + '.pdf')
plt.show()

#fig2, axs2 = plt.subplots(1, 1)


#p = pScatter(os.system("python3 allrun.py"))
#p.plot()
#p.update()

with open(plotName + '.obj', 'wb') as file:
	pickle.dump(fig, file)

plt.show()
plt.close(fig)


