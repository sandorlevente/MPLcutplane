#usual import
#you will need need matplotlib, numpy, and VTK
import sys;
import numpy as np
import vtk
import os
from vtk.util.numpy_support import vtk_to_numpy
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import matplotlib.tri as triangulator

######################## utility definitions ############################

#Getting the file list from directory
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

# Reading vtk files
def vtkRead(fileName):
    '''
    Read the data from the file. 
    '''

    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(fileName)
    reader.Update()
    data = reader.GetOutput()
    return (data)

#Reading the triangulated vtk scalar data
def triangulateScalarData(data, val):
    '''
    Triangulates the data file and return the grid points in (x, y, z), 
    the connectivity (tri) and the SCALAR values (scalar):
        data:       VTK data output
        x, y, z:    coordinates,
        tri:        triangulated connections,
        val:        scalar value. 
    '''

    cells = data.GetPolys()
    triangles = cells.GetData()
    points = data.GetPoints()
    point_data = data.GetPointData().GetArray(val)

    ntri = int(triangles.GetNumberOfTuples() / 4)
    npts = points.GetNumberOfPoints()
    nvls = point_data.GetNumberOfTuples()

    tri = np.zeros((ntri, 3))
    x = np.zeros(npts)
    y = np.zeros(npts)
    z = np.zeros(npts)
    scalar = np.zeros(nvls)

    for i in range(0, ntri):
        tri[i, 0] = triangles.GetTuple(4 * i + 1)[0]
        tri[i, 1] = triangles.GetTuple(4 * i + 2)[0]
        tri[i, 2] = triangles.GetTuple(4 * i + 3)[0]

    for i in range(points.GetNumberOfPoints()):
        pt = points.GetPoint(i)
        x[i] = pt[0]
        y[i] = pt[1]
        z[i] = pt[2]

    scalar = vtk_to_numpy(data.GetPointData().GetArray(val))
    return (x, y, z, tri, scalar)

#Reading triangulated vector data
def triangulateVectorData(data, val):
    '''
    Triangulates the data file and return the grid points in (x, y, z), 
    the connectivity (tri) and the VECTOR values (ux, uy, uz):
        data:       VTK data output, 
        val:        name of vector. 
    Return: 
        x, y, z:    coordinates,
        tri:        triangulated connections,
        ux, uy, uz: vector values. 
    '''

    cells = data.GetPolys()
    triangles = cells.GetData()
    points = data.GetPoints()
    point_data = data.GetPointData().GetArray(val)

    ntri = int(triangles.GetNumberOfTuples() / 4)
    npts = points.GetNumberOfPoints()
    nvls = point_data.GetNumberOfTuples()

    tri = np.zeros((ntri, 3))
    x = np.zeros(npts)
    y = np.zeros(npts)
    z = np.zeros(npts)
    ux = np.zeros(nvls)
    uy = np.zeros(nvls)
    uz = np.zeros(nvls)

    for i in range(0, int(ntri)):
        tri[i, 0] = triangles.GetTuple(4 * i + 1)[0]
        tri[i, 1] = triangles.GetTuple(4 * i + 2)[0]
        tri[i, 2] = triangles.GetTuple(4 * i + 3)[0]

    for i in range(points.GetNumberOfPoints()):
        pt = points.GetPoint(i)
        x[i] = pt[0]
        y[i] = pt[1]
        z[i] = pt[2]

    for i in range(0, point_data.GetNumberOfTuples()):
        U = point_data.GetTuple(i)
        ux[i] = U[0]
        uy[i] = U[1]
        uz[i] = U[2]

    return (x, y, z, tri, ux, uy, uz)

#Getting the vector magnitude
def velocityMag(data, nameOfVec):
    '''
    Calculate the magnitude of a vector:
        data:       VTK data output, 
        nameOfVec:  name of the vector field (for instance: 'U').
    Return: 
        x, y, z:    grid coordinates, 
        tri:        connectivity triangles, 
        magnitude:  magnitude of the vector at given grid points. 
    '''

    (x, y, z, tri, u, v, w) = triangulateVectorData(data, nameOfVec)
    magnitude = np.sqrt(u ** 2.0 + v ** 2.0 + w ** 2.0)
    return (x, y, z, tri, magnitude)

#Getting the vector magnitude
def planarMag(data, nameOfVec, vector1, vector2):
    '''
    Calculate the magnitude of a vector:
        data:       VTK data output, 
        nameOfVec:  name of the vector field (for instance: 'U').
    Return: 
        x, y, z:    grid coordinates, 
        tri:        connectivity triangles, 
        magnitude:  magnitude of the vector at given grid points. 
    '''        
    (x, y, z, tri, u, v, w) = triangulateVectorData(data, nameOfVec)
    coords = {}
    coords['X'] = u
    coords['Y'] = v
    coords['Z'] = w
    magnitude = np.sqrt(coords[vector1] ** 2.0 + coords[vector2] ** 2.0)
    return (x, y, z, tri, magnitude)



############### Plotting definitions #############################


def  scalarPlot(fileString, quantity, xVectorDir, yVectorDir, xRange, yRange,ticks, minValue, maxValue, noColors, colorMap, cbOrientation, quantityName):
    fileName = getListOfFiles('../postProcessing/triangulatedCuttingPlanes/')
    fileName = [counter for counter in fileName if fileString in counter]
    for i in range(0, len(fileName)):    
        data = vtkRead(fileName[i])
        print(fileName[i], 'loaded')
        
        # select quantity
        x,y,z,tri,scalar = triangulateScalarData(data, quantity)
        print(f'Min and max {quantity} for {fileName[i]}: {np.min(scalar)}, {np.max(scalar)}')
        # transform data into arrays
        xNew = []
        yNew = []
        zNew = []
        #convert to mm
        for j in range(0, len(x)):
            xNew.append(x[j] * 1e+3)
            yNew.append(y[j] * 1e+3)
            zNew.append(z[j] * 1e+3)
        x = xNew
        y = yNew
        z = zNew

        coords = {}
        coords['X'] = x
        coords['Y'] = y
        coords['Z'] = z

        # ranges and labels
        xMinMax = [min(x), max(x)]
        yMinMax = [min(y), max(y)]
        zMinMax = [min(z), max(z)]
        xLabel = 'X [mm]'
        yLabel = 'Y [mm]'
        zLabel = 'Z [mm]'
        #path assembly
        figName=os.path.join('../postProcessing/figures/', '{}{}'.format((os.path.splitext(os.path.basename(fileName[i]))[0]),'.png'))
        
        fig1, ax1 = plt.subplots()
        fig_levels = np.linspace(minValue, maxValue, noColors)
        plt.xlim(xRange)
        plt.ylim(yRange)
        plt.rc('font', family='serif', size=5)
        plt.title(os.path.basename(fileName[i]))

        ax1.set_aspect('equal')
        ax1.xaxis.set_major_locator(plt.MultipleLocator(ticks))
        ax1.yaxis.set_major_locator(plt.MultipleLocator(ticks))
        #contour filled
        tcf = ax1.tricontourf(coords[xVectorDir],coords[yVectorDir],tri,scalar,fig_levels,cmap=colorMap, extend='both')
        #contour line
        ax1.tricontour(coords[xVectorDir],coords[yVectorDir],tri, scalar,fig_levels, colors='black', linewidths=0.1)
        
        clb=fig1.colorbar(tcf,orientation=cbOrientation)
        clb.set_label(quantityName)

        plt.tight_layout()

        fig1.savefig(figName, format='png', dpi=300)
        
        print('Figure created for', fileName[i])

def  vectorPlot(fileString, quantity, xVectorDir, yVectorDir, xRange, yRange,ticks, minValue, maxValue, noColors, colorMap,cbOrientation, quantityName):
    # getting file list
    fileName = getListOfFiles('../postProcessing/triangulatedCuttingPlanes/')
    fileName = [counter for counter in fileName if fileString in counter]
    for i in range(0, len(fileName)):    
        data = vtkRead(fileName[i])
        print(fileName[i], 'loaded')
        
        # select quantity
        x, y, z, tri, ux, uy, uz = triangulateVectorData(data, quantity)
        x,y,z,tri,velMag=planarMag(data, quantity, xVectorDir, yVectorDir)
        print(f'Min and max {quantity} for {fileName[i]}: {np.min(velMag)}, {np.max(velMag)}')
        # transform data into arrays
        xNew = []
        yNew = []
        zNew = []
        #convert to mm
        for j in range(0, len(x)):
            xNew.append(x[j] * 1e+3)
            yNew.append(y[j] * 1e+3)
            zNew.append(z[j] * 1e+3)
        x = xNew
        y = yNew
        z = zNew

        coords = {}
        coords['X'] = x
        coords['Y'] = y
        coords['Z'] = z

        # ranges and labels
        xMinMax = [min(x), max(x)]
        yMinMax = [min(y), max(y)]
        zMinMax = [min(z), max(z)]
        xLabel = 'X [mm]'
        yLabel = 'Y [mm]'
        zLabel = 'Z [mm]'
        #path assembly
        figName=os.path.join('../postProcessing/figures/', '{}{}'.format((os.path.splitext(os.path.basename(fileName[i]))[0]),'.png'))
        
        fig1, ax1 = plt.subplots()
        fig_levels = np.linspace(minValue, maxValue, noColors)
        plt.xlim(xRange)
        plt.ylim(yRange)
        plt.rc('font', family='serif', size=5)
        plt.title(os.path.basename(fileName[i]))

        ax1.set_aspect('equal')
        ax1.xaxis.set_major_locator(plt.MultipleLocator(ticks))
        ax1.yaxis.set_major_locator(plt.MultipleLocator(ticks))
        #contour filled
        tcf = ax1.tricontourf(coords[xVectorDir],coords[yVectorDir],tri,velMag,fig_levels,cmap=colorMap, extend='both')
        #contour line
        ax1.tricontour(coords[xVectorDir],coords[yVectorDir],tri, velMag,fig_levels, colors='black', linewidths=0.1)
        
        clb=fig1.colorbar(tcf,orientation=cbOrientation)
        clb.set_label(quantityName)

        plt.tight_layout()

        fig1.savefig(figName, format='png', dpi=300)
        
        print('Figure created for', fileName[i])


############################## DUMMY PLOT (because the fonts are not working on the first one....) ############################
scalarPlot('p_V', 'p','X','Z', [50,250],[0, 150],20, 99992, 100003, 20, 'jet', 'horizontal', 'Static pressure [Pa]')




############################## Defining plots ############################



#scalarPlot(fileString, quantity, xVectorDir (for plot), yVectorDir (for plot), xRange (for plot), yRange(for plot) ,ticks, minValue, maxValue, noColors, colorMap, cbOrientation, quantityName)
scalarPlot('p_V', 'p','X','Z', [50,250],[0, 150],20, 99992, 100003, 20, 'jet', 'horizontal', 'Static pressure [Pa]')


#vectorPlot(fileString, quantity,  xVectorDir (for plot), yVectorDir (for plot), xRange (for plot), yRange(for plot) ,ticks, minValue, maxValue, noColors, colorMap,cbOrientation, quantityName)
vectorPlot('U_H', 'U', 'X', 'Y', [50,250], [-200, -50],50, 0, 10, 20, 'jet', 'horizontal', 'vectorplot')