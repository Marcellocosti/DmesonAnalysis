import os
import sys
import ROOT
from ROOT import TDirectoryFile, TList

def Load(container, path):
    '''
    Function to extract an object inside a root file.
    Supports nested containers with the following Data Types:
     - TDirectoryFile (TFile)
     - TList

    Parameters
    -----------
    container: TFile of the input file
    path: path of the object inside the root file

    Returns:
    -----------
    obj: target root object
    '''

    # Check that the input file is OK
    path = os.path.normpath(path)
    if container == None:  # pylint: disable=singleton-comparison
        print('The container %s is NULL', container.GetName())

    # Start to extract
    for name in path.split(os.sep):
        print('Trying to load %s:%s. Available keys:', container.GetName(), name)

        for key in GetKeyNames(container):
            print('    %s', key)

        if isinstance(container, TDirectoryFile):
            obj = container.Get(name)
        elif isinstance(container, TList):
            obj = container.FindObject(name)
        else:
            print('The container %s of type %s is not valid', container.GetName(), type(container))

        if obj == None:  # pylint: disable=singleton-comparison
            print('The container %s does not contain an object named %s', container.GetName(), name)
            raise NameError()
        container = obj

    print('The object %s:%s was succesfully loaded', container.GetName(), path)
    return obj

def GetKeyNames(container):  # pylint: disable=inconsistent-return-statements
    if isinstance(container, TDirectoryFile):
        return [key.GetName() for key in list(container.GetListOfKeys())]

    if isinstance(container, TList):
        it = container.MakeIterator()
        names = []
        while True:
            obj = it.Next()
            if obj == None: # pylint: disable=singleton-comparison
                break
            names.append(obj.GetName())

        return names
    print('Unknown container type %s', type(container))
    
def GetNPanels(n):
    if n<=3:
        return (n, 1)
    elif n == 4:
        return (2, 2)
    elif n <= 8:
        return ((n+1)//2, 2)
    elif n <= 12:
        return ((n+1)//3, 3)
    elif n <= 20:
        return ((n+1)//5, 4)
    else:
        sys.exit()