from pyFSI.post.prints import printError
import sys


def error(message):
    printError(message)

def errorExit(message):
    printError(message)
    sys.exit()
    
