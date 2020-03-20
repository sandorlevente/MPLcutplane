#   paraview cut triangulator script
#   run it with pvpython from the case/scripts directory

import paraview.simple as pvs
import os

def getListOfFiles(dirName):
    # create a list of file and sub directories 
    # names in the given directory 
    listOfFile = os.listdir(dirName)
    allFiles = list()
    # Iterate over all the entries
    for entry in listOfFile:
        # Create full path
        fullPath = os.path.join(dirName, entry)
        # If entry is a directory then get the list of files in this directory 
        if os.path.isdir(fullPath):
            allFiles = allFiles + getListOfFiles(fullPath)
        else:
            allFiles.append(fullPath)
                
    return allFiles



fileList = getListOfFiles('../postProcessing/cuttingPlanes/')
fileList = [el for el in fileList if '.vtk' in el]

'''
directory='../postProcessing/triangulatedCuttingPlanes'
if not os.path.exists(directory):
    os.makedirs(directory)
'''

for i in range(0, len(fileList)):
    pvdata=pvs.LegacyVTKReader(FileNames=fileList[i])
    clean1=pvs.Clean(Input=pvdata)
    clean1.PieceInvariant = 0
    clean1.AbsoluteTolerance = 0.0000001
    clean1.ToleranceIsAbsolute = 1
    clean1.ConvertLinesToPoints = 0
    clean1.ConvertPolysToLines = 0
    clean1.ConvertStripsToPolys = 0
    triang = pvs.Triangulate(Input=clean1)
    pvs.SaveData(filename=os.path.join('../postProcessing/triangulatedCuttingPlanes/', '{}_{}{}'.format((os.path.splitext(os.path.basename(fileList[i]))[0]),os.path.basename(os.path.dirname(fileList[i])),(os.path.splitext(os.path.basename(fileList[i]))[1]))), proxy=triang, FileType='Ascii')
    print(fileList[i],'triangulated')