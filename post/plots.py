# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt, numpy as np
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.font_manager
import scipy.signal as sps
from matplotlib.animation import FuncAnimation
import pandas as pd
from multiprocessing import Pool
import scipy as sp
from colorama import Fore, Style


plt.rcParams.update({
    "text.usetex":              True,  
    "font.family":              "serif",
    "font.serif":               "Computer Modern Roman",
    "font.size":                20,
    "figure.figsize":           [16, 9],
    "axes.grid":                True,
    "savefig.format":           "pdf",
    "animation.writer":         "ffmpeg",
    "animation.ffmpeg_args":    ['-c:v', 'libx264', '-crf', '0'],
    "animation.bitrate":        24000,
    "animation.frame_format":   "png",
    "legend.loc":               "lower left"

})


def parLoadFiles(fileList):
    n = len(fileList)
    pool = Pool(2)
    result = pool.map(np.loadtxt, fileList)
    pool.close()
    return result

def loadFiles(fileList):
    n = len(fileList)
    result = []
    for file in fileList:
        try:
            print(Fore.YELLOW + "---> Loading file: ", file)
            data = np.loadtxt(file)
        except:
            data = np.loadtxt(file, dtype=complex)
        result.append(data)

    return result


# Load files in parallel
class fsiPlot():
    def __init__(self, fig=None, axe=None):
        self.fig = fig
        self.axe = axe
        self.xData = None
        self.yData = None
        self.P = [] # list of plots
        # plt.rcParams.update({
        #     "text.usetex": True,
        #     "font.family": "sans-serif",
        #     "font.sans-serif": ["Helvetica"]})
        ## for Palatino and other serif fonts use:
        # plt.rcParams.update({
        # "text.usetex": True,
        # "font.family": "sans-serif",
        # "font.sans-serif": ["Helvetica"]})
        # for Palatino and other serif fonts use:
        plt.rcParams.update({
            "text.usetex": False,
            "font.family": "serif",
            "font.serif": ["Palatino"],
            "font.size": 12})
        # It's also possible to use the reduced notation by directly setting font.family:
        # plt.rcParams.update({
        # "text.usetex": True,
        # "font.family": "Helvetica"
        # })

    def set(self,
            fig=False,
            axe=False,
            title=False,
            xLim=False,
            yLim=False,
            xLabel=False,
            yLabel=False):

        if fig:
            self.fig = fig
        if axe:
            self.axe = axe
        if xLim:
            self.axe.set_xlim(xLim[0], xLim[1])
        if yLim:
            self.axe.set_ylim(yLim[0], yLim[1])
        if title:
            self.axe.set_title(title)
        if xLabel: 
            self.axe.set_xlabel(xLabel)
        if yLabel:
            self.axe.set_ylabel(yLabel)


class plotFromFile(fsiPlot):
    def __init__(self, xFile, yFile, xIndexes=None, yIndexes=None, 
                scatter=False, type=None, fig=False, axe=False, scale=False,
                color=None):
        super().__init__(fig, axe)

        if self.fig == False and self.axe == False:
            self.fig, self.axe = plt.subplots(1)

        # Get the data
        data = loadFiles([xFile, yFile])

        if type == "Eigenvalues":  # Reading complex number files
            data[0] = np.real(data[0])
            data[1] = np.imag(data[1])
        if type == "parametricReal":
            data[1] = np.real(data[1])
        if type == "parametricImag":
            data[1] = np.imag(data[1])
        self.P = []
        
        # Build the color scheme
        if color is not None:
            tmp = loadFiles([color])[0]
            if tmp.ndim == 1:
                colorArray = tmp
            elif tmp.ndim == 2: # Build the color with the first column
                colorArray = tmp[:, 0]
            cmap = plt.get_cmap('jet', 10)
            norm = mpl.colors.Normalize(vmin=np.min(colorArray), vmax=np.max(colorArray))
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm) 
            sm.set_array([])
        

        # Number of plots
        try:
            nPlots = abs(yIndexes[1].stop - yIndexes[1].start)
        except:
            nPlots = 1

        # Number of x Axes
        try:
            xAxes = abs(xIndexes[1].stop - xIndexes[1].start)
        except:
            xAxes = 1

        # Assign the data to the local attributes
        self.xData = np.array(data[0])[xIndexes]
        self.yData = np.array(data[1])[yIndexes]

        # Scale the y axis
        if scale:
            self.yData *= scale

        # Solve the problem of different lengths for an incomplete run
        y1D = False
        try:
            xLength = len(self.xData[:, 0])
        except:
            x1D = True
            xLength = len(self.xData[:])
        try:
            yLength = len(self.yData[:, 0])
        except:
            y1D = True
            yLength = len(self.yData[:])

        if xLength > yLength:
            if x1D:
                self.xData = self.xData[0:yLength]
            else:
                self.xData = self.xData[0:yLength,:]
        if yLength > xLength:
            if y1D:
                self.yData = self.yData[0:xLength]
            else:
                self.yData = self.yData[0:xLength,:]

       #
        #
        # if xAxes == 1:
        #     xLength = len(self.xData[:])
        #     yLength = len(self.yData[:, 0])
        #     if xLength > yLength:
        #         self.xData = self.xData[0:yLength]
        #     if yLength > xLength:
        #         self.yData = self.yData[0:xLength, :]
        # else:
        #     xLength = len(self.xData[:, 0])
        #     yLength = len(self.yData[:, 0])
        #     if xLength > yLength:
        #         self.xData = self.xData[0:yLength, :]
        #     if yLength > xLength:
        #         self.yData = self.yData[0:xLength, :]

        if xAxes == 1 and nPlots == 1:
            if scatter:
                curve = self.axe.scatter(self.xData, self.yData, c=colorArray)
                plt.colorbar(label="Flow Rate", orientation="horizontal")
                #self.ax.set_title("scale")
            else:
                curve, = self.axe.plot(self.xData, self.yData)
            self.P.append(curve)

        else:
            if xAxes == 1:
                for i in range(nPlots):
                    if scatter:
                        curve = self.axe.scatter(self.xData, self.yData[:, i], s=20)
                    else:
                        curve, = self.axe.plot(self.xData, self.yData[:, i])
                    self.P.append(curve)
            else:
                for i in range(nPlots):
                    if scatter:
                        curve = self.axe.scatter(self.xData[:, i], self.yData[:, i], c=colorArray, cmap='jet')

                    else:
                        curve, = self.axe.plot(self.xData[:, i], self.yData[:, i], linewidth=4)
                    self.P.append(curve)
        if scatter:
            plt.colorbar(sm)


    def update(self, i):
        for c in self.yIndexes:
            self.P[c].set_xdata(self.xData[self.xIndexes[0]:i])
            self.P[c].set_ydata(self.yData[self.xIndexes[0]:i, c])
        return self.P

# Deformed shape plot
class pShape(fsiPlot):
    def __init__(self, xfile, yfile, time=-1, fig=False, axe=False, scale=False,
                shadows=False, freq=None, dt=None):
        super().__init__(fig, axe)
        if self.fig == False and self.axe == False:
            self.fig, self.axe = plt.subplots(1)

        self.xData = np.loadtxt(xfile)
        self.yData = np.loadtxt(yfile)
        if scale:
            self.yData *= scale

        if shadows:
            istart = time
            ifinal = istart - 1 / (freq * dt)
            instants = np.linspace(ifinal, istart, shadows).astype(int)
            # Draw the shape at t=0
            self.P = self.axe.plot(self.xData[0, :],
                        self.yData[0, :],
                        linewidth=3,
                        color='black')
            # Draw the shadow shapes
            for i in instants:
                self.P.append(self.axe.plot(self.xData[:],
                                            self.yData[i, :],
                                            linewidth=2,
                                            color='grey'))
        else:
            self.P = self.axe.plot(self.xData[:],
                                    self.yData[time, :],
                                    linewidth=3,
                                    color='black')



    def update(self, i):
        self.P[0].set_ydata(self.yData[i, :])
        return self.P[0]


class pFlow(fsiPlot):
    def __init__(self, xfile, yfile, time):
        super().__init__()
        self.fig, self.axe = plt.subplots(1)
        self.xData = np.loadtxt(xfile)
        self.yData = np.loadtxt(yfile)
        self.P = self.axe.plot(self.xData[time, :],
                               self.yData[time, :],
                               linewidth=5,
                               color='black')

    def update(self, i):
        self.P[0].set_ydata(self.yData[i, :])
        return self.P[0]

# Scatter plot
class pScatter(fsiPlot):
    def __init__(self, solution):
        # pidx is the parameter index of the solution list
        # tidx is the time  index
        super().__init__()
        if solution[0][0].execution()['isParametric']:
            para = solution[0][0].execution()['parameters']
        else:
            para = {'p': [1], 'pi': [1], 'steps': 1}

        time = solution[0][0].execution()['time']
        esize = solution[0][0].ES.size()
        self.colors = cm.rainbow(np.linspace(0, 1, esize))
        self.data = np.zeros((para['steps'], time['steps'], esize, 2))
        self.evec = np.zeros((para['steps'], time['steps'], esize, esize),
                             dtype=complex)
        self.P = self.axe.scatter([], [], edgecolors='k', facecolors='r', s=10)
        for pi, ppf in enumerate(para['p']):  # Parameter loop
            for ti, t in enumerate(time['t']):  # Time loop
                evalues = solution[pi][ti].ES.evalues()
                evectors = solution[pi][ti].ES.evectors()
                # Store real and imaginary parts of eigenvalues
                self.data[pi, ti, :, 0] = np.real(evalues)
                self.data[pi, ti, :, 1] = np.imag(evalues)
                self.evec[pi, ti, :, :] = evectors

        self.esize = esize
        self.tsize = len(time['t'])

    def plot(self, i):
        x = self.data[i, :, :, 0].flatten()
        y = self.data[i, :, :, 1].flatten()
        self.P = self.axe.scatter(x, y, edgecolors='k', facecolors='r', s=10)
        self.fig.show()

    def update(self, i):
        x = self.data[i, :, :, 0].flatten()
        y = self.data[i, :, :, 1].flatten()
        data = np.array([x, y]).T  # The transpose MUST be called
        self.P.set_offsets(data)
        self.fig.show()

class pMat(fsiPlot):
    def __init__(self, files):
        super().__init__()
        self.fig, self.axe = plt.subplots(2)
        for file in files:
            self.data = np.genfromtxt(file, delimiter="\t")
            mu0 = 1.256E-6
            self.fig.suptitle('Material properties for '  + file.name)
            # Plot the B-H curve
            self.axe[0].grid(b=True)
            self.axe[0].plot(self.data[1:, 1], self.data[1:, 0])
            #self.axe[0].set_title('B-H Curve')
            self.axe[0].set_xlabel('H [A/m]')
            self.axe[0].set_ylabel('B [T]')


            # Plot the B-mur curve
            mur = (self.data[1:, 0] / self.data[1:, 1]) / mu0
            self.axe[1].grid(b=True)
            self.axe[1].plot(self.data[1:, 0], mur)
            #self.axe[1].set_title('B-mur Curve')
            self.axe[1].set_xlabel('B [T]')
            self.axe[1].set_ylabel('mur')

        plt.show()
        
class pFFT(fsiPlot):
    def __init__(self, time, signal, xIndexes, yIndexes, fig=False, axe=False,):
        super().__init__(fig, axe)
        if self.fig == False and self.axe == False:
            self.fig, self.axe = plt.subplots(1)
        self.xData = np.loadtxt(time)[xIndexes]
        self.yData = np.loadtxt(signal)[yIndexes]
        N = len(self.xData)
        # sample spacing
        T = (self.xData[1] - self.xData[0])
        yf = sp.fft.fft(self.yData)
        xf = sp.fft.fftfreq(N, T)[:N//2]
        self.P = self.axe.plot(xf, 2.0/N * np.abs(yf[0:N//2]))
        
class pModes(fsiPlot):
    def __init__(self, meshFile, evaluesFile, evectorsFile, fig=False, axe=False):
        super().__init__(fig, axe)
        data = loadFiles([meshFile, evaluesFile, evectorsFile])
        mesh = data[0]
        evalues = data[1]
        evectors = data[2]
        for i, value in enumerate(evalues):
            freq = value / 6.28
            label = str(freq.round(2)) + " Hz"
            curve, = self.axe.plot(mesh, evectors[i], label=label, linewidth=4)
            self.P.append(curve)
        self.axe.legend()
        self.fig.show()
        
# Class for plotting magnetic flux function taken from FEMM with lua script
class pFlux(fsiPlot):
    def __init__(self, filePath, coor=1, deriv=True, filter=True):
        super().__init__()
        self.data = np.genfromtxt(filePath, delimiter=",", skip_header=3)

        # Flux
        self.aflx = np.zeros((len(self.data[:, 0]), 3))
        self.aflx[:, 0] = self.data[:, 0]
        self.aflx[:, 1] = sps.savgol_filter(self.data[:, coor], 11, 2, deriv=0)
        self.aflx[:, 2] = self.data[:, coor]
        # Derivative of the flux
        self.dflx = np.zeros((len(self.aflx[:, 0]), 2))  # Derivative of average flux
        self.dflx[:, 0] = self.aflx[:, 0]  # LENGTH COORDINATE

        if filter:
            temp = sps.savgol_filter(self.data[:, 1], 51, 2, deriv=1)
            # Necesito el dx para dividir la derivada porqeu sale calculada tomando dx=1 por defecto
            dx = self.aflx[1, 0] - self.aflx[0, 0]
            # Ojo que tira la derivada cambiada de signo
            self.dflx[:, 1] = -sps.savgol_filter(temp, 51, 2, deriv=0) / dx
        else:
            self.dflx[:, 1] = np.gradient(self.aflx[:, 1], self.aflx[:, 0], edge_order=2)

        # Choose the flux function or the derivative
        if deriv:
            curve, = plt.plot(self.dflx[:, 0], self.dflx[:, 1])
        else:
            curve, = plt.plot(self.aflx[:, 0], self.aflx[:, 2])
        self.axe = plt.gca()
        self.P = []
        self.P.append(curve)


class pBifurcation(fsiPlot):
    def __init__(self, solution, mode0, mode1, var='freq'):
        super().__init__()
        #self.P = self.axe.scatter([], [], edgecolors='k', facecolors='r', s=10)
        self.bif = []  # Bifurcation points list
        size = solution.fsi().ES.size

        for sol in solution:  # Loop through the parameters
            oldReals = -np.ones(size)  # Assume stability
            searchModes = range(mode0, mode1)
            for fsi in sol:
                newReals = np.real(fsi.ES.values())
                newImags = np.imag(fsi.ES.values())
                signs = np.sign(newReals) + np.sign(oldReals)

                # Extract all unstable states
                for mode in searchModes:
                    if newReals[mode] > 0 and newImags[mode] > 0:
                        xVar = fsi.execution()['parameters']['pi']
                        if var == 'frequency':
                            yVar = newImags[mode]
                        if var == 'flowRate':
                            yVar = fsi.flow().qx0
                        elif var == 'velocity':
                            yVar = fsi.flow().vRef
                        else:
                            yVar = None
                            print(Fore.GREEN + "---> ERROR: No yVar defined...")
                        self.bif.append(bifurcationPoint(xVar, yVar, mode))

                # Extrac the stability boundary
                # for mode in searchModes:
                #     if signs[mode] == 0 and np.sign(oldReals[mode]) != 0 and newImags[mode] > 0:
                #         self.bif.append(bifurcationPoint(t.pi, newImags[mode], mode))
                #         #searchModes = [mode] # If I found one unstable, follow this
                oldReals = newReals

        # Get the data from the bifurcation point
        self.xAxis = []
        self.yAxis = []
        for p in self.bif:
            self.xAxis.append(p.xParameter)
            self.yAxis.append(p.yParameter)

    def plot(self):
        self.axe.scatter(self.xAxis,
                         self.yAxis,
                         edgecolors='k',
                         facecolors='r',
                         s=10)
        plt.show()


class bifurcationPoint():
    def __init__(self, xPar, yPar, neval):
        self.xParameter = xPar
        self.yParameter = yPar
        self.modeNumber = neval


