import pandas as pd
from tabulate import tabulate
import numpy as np
import pprint, sys
from colorama import Fore, Back, Style

def printEigen(values, vectors):
    pd.set_option('display.max_columns', 12)
    pd.options.display.float_format = "{:,.2f}".format
    for i, val in enumerate(values):
        print("+--------------+------------+")
        print("+         Eigenpair " + str(i+1) + "       +")
        # dataval = {"Re": [np.real(values[i])],
        #            "Im": [np.imag(values[i])]}
        # dfval = pd.DataFrame(data=dataval)
        # print(tabulate(dfval, headers='keys', tablefmt='psql'))
        #"Second argument: {1:3d}, first one: {0:7.2f}".format(47.42,11)
        ReStr = "Real\n{0:6.2e}".format(np.real(values[i]))
        ImStr = "Imaginary\n{0:5.2e}".format(np.imag(values[i]))

        datavec = {ReStr: np.real(vectors[i]),
                   ImStr: np.imag(vectors[i])}
        dfvec = pd.DataFrame(data=datavec)
        print(tabulate(dfvec,
                       headers='keys',
                       tablefmt='psql',
                       floatfmt=".5f",
                       showindex="never",
                       numalign="decimal",
                       stralign="center"))


def printStateMatrix(solution):
    pd.set_option('display.max_columns', 12)
    pd.options.display.float_format = "{:,.6f}".format
    print("+--------------+------------+")
    print("+      State Matrix         +")
    dfvec = pd.DataFrame(data=solution.S)
    print(tabulate(dfvec,
                   headers='keys',
                   tablefmt='psql',
                   floatfmt=".6f",
                   showindex=True,
                   numalign="decimal",
                   stralign="center"))


def printNumbers(solution):
    flow = solution.flow()
    nflow = list(flow.dimNumbers.values())
    for n in nflow:
        n.info()

    solid = solution.solid()
    nsolid = list(solid.dimNumbers.values())
    for n in nsolid:
        n.info()

    fsi = solution
    nfsi = list(fsi.dimNumbers.values())
    for n in nfsi:
        n.info()


def printArray(array, name=False, complex=False, file=False):
    dfvec = pd.DataFrame(array)
    if name:
        header = "+----+-----------+------------+------------+"
        spaces = int((len(header) - len(name)) / 2)
        print(header, "\n", " " * spaces, name)
        
    if complex:     
        np.set_printoptions(precision=3)
        print(array)
    else:
        table = tabulate(dfvec,
                        headers='keys',
                        tablefmt='psql',
                        floatfmt=".3f",
                        showindex=True,
                        numalign="decimal",
                        stralign="center")
        if file is False:
            print(table) # Just print in the current stdout
        else:
            print(table) # Print in the current stdout
            # Save to file
            with open(file, 'w') as f:
                original_stdout = sys.stdout # Save a reference to the original standard output
                sys.stdout = f # Change the standard output to the file we created.
                print(table)
                sys.stdout = original_stdout
    

def printNumbers(obj):
    for key, val in obj.dimNumbers.items():
        pprint.pprint(val.info())

def printError(message, width=69):
    stars = '*' * width
    pad = (width + len(message)) // 2
    print(Fore.RED + f'{stars}\n{message:>{pad}}\n{stars}')

def printInfo(object, message):
    if type == 'solid':  
        print(Fore.BLUE + object.name  + "Info ---> " + message)