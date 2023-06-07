# trace generated using paraview version 5.10.1
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 10

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

timefrequency = 0
filenames = []

for i in range(6):
    filename = 'C:\\Users\\s338090\\OneDrive - Cranfield University\\LBM\\taichi-LBM\\taichi_LBM3D\\2phase\\bubble-case\\structured' + str(timefrequency) + '.vtr'
    filenames.append(filename)
    timefrequency += 1

structured=XMLRectilinearGridReader(registrationName='structured*', FileName=filenames)