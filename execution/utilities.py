import sys, pathlib
import numpy as np
from colorama import Fore

def banner(s, width=84):
    stars = '*' * width
    pad = (width + len(s)) // 2
    print(Fore.LIGHTBLUE_EX + f'{stars}\n{""}\n{s:>{pad}}\n{""}\n{stars}')
    
def section(s, width=84):
    delimiter = '-' * width
    pad = (width + len(s)) // 2
    print(Fore.LIGHTYELLOW_EX + f'{delimiter}\n{s:>{pad}}\n{delimiter}')
    
def subSection(s, width=84, where='low'):
    delimiter = '-' * width
    pad = (width + len(s)) // 2
    if where == 'low':
        print(Fore.LIGHTYELLOW_EX + f'\n{s:>{pad}}{delimiter}\n')
    else:
        print(Fore.LIGHTYELLOW_EX + f'\n{delimiter}\n{s:>{pad}}\n')

def writeVersionHeader(file):
    version = "pyFSI Output File - Martin Saravia - 2021"
    head = "#" + (67 * "-") + "#\n"
    pad = ((67 - len(version)) // 2) * " "
    message = "#" + pad + version + pad + "#"
    file.write(head)
    file.write(message + "\n")
    file.write(head)
    
def file(name):
    filePath = pathlib.Path(__file__).parent.absolute()/name
    return filePath

def setParameter(control, time):
    command = control['parameter']['command']
    values = control['parameter']['value']
    setValue = time.pf * (values[1] - values[0]) + values[0]   
    exec("control" + command + " = setValue")
    
def loadFiles(files):
    if isinstance(files, list):
        n = len(files)
        dataList = []
        for file in files:
            try:
                data = np.loadtxt(file)
            except:
                data = np.loadtxt(file, dtype=np.complex)
            dataList.append(data)
        return dataList
    else:
        try:
            data = np.loadtxt(files)
        except:
            data = np.loadtxt(files, dtype=np.complex)
    return data
    
    
def getObjectbyName(name, list):
    for obj in list: 
        if obj == name:
            return list[obj]
        
    

# def getBifurcations(eigenValuesFile, eigenVectorsFile=None):
#     values = loadFiles(eigenValuesFile)
#     #vectors = loadFiles(eigenVectorsFile)
#     size = len(values[0])
#     realValues = np.real(values)
#     realValuesSigns = np.zeros_like(realValues)
#     imagValues = np.imag(values)
    
#     realValuesSigns = np.sign(realValues)
#     for i in range(size):
#         row = realValues[:,i]
#         realValuesSigns[:,i] = ((np.roll(row, 1) - row) != 0).astype(int)
#         realValuesSigns[0,i] = 0

#     print(realValues)
#     print(realValuesSigns[:,i])
    

