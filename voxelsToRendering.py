

import vtk

import sys

from scipy.ndimage import morphology

sys.path.append('/Users/paulaljabar/Python/nibabel-1.3.0-py2.7.egg')
import nibabel as nib

import numpy


def usage(scriptName):
  print '''Usage: '

  vtkpython ''',
  print  scriptName,
  print ''' [image-input] [polydata-output]

Identify the boundary voxels between 0 and postive voxels in the input image.

Render a cube at each of the voxels and save the result in a polydata output file.
'''

  exit()


def main(*args):
  if len(args) < 3:
     usage(args[0])

  inputFileName = args[1]
  outputFileName = args[2]
  
  
  img = nib.load(inputFileName)
  imgData = img.get_data()
  
  aff = img.get_affine()
  
  pixdims = img.get_header()['pixdim'][1:4]

  pd = vtk.vtkPolyData()
  pts = vtk.vtkPoints()

  # Identify boundary voxels:

  voxelInds = numpy.where(imgData > 0)
  
  boundaryVoxels = numpy.zeros(imgData.shape)

  temp = numpy.zeros(imgData.shape)
  temp[:-1, :, :] = imgData[1:, :, :]
  inds = numpy.where(numpy.bitwise_and(imgData > 0 , temp == 0))
  boundaryVoxels[inds] = 1

  temp[:] = 0
  temp[1:, :, :] = imgData[:-1, :, :]
  inds = numpy.where(numpy.bitwise_and(imgData > 0 , temp == 0))
  boundaryVoxels[inds] = 1
  
  temp[:] = 0
  temp[:, :-1, :] = imgData[:, 1:, :]
  inds = numpy.where(numpy.bitwise_and(imgData > 0 , temp == 0))
  boundaryVoxels[inds] = 1

  temp[:] = 0
  temp[:, 1:, :] = imgData[:, :-1, :]
  inds = numpy.where(numpy.bitwise_and(imgData > 0 , temp == 0))
  boundaryVoxels[inds] = 1
  
  temp[:] = 0
  temp[:, :, :-1] = imgData[:, :, 1:]
  inds = numpy.where(numpy.bitwise_and(imgData > 0 , temp == 0))
  boundaryVoxels[inds] = 1

  temp[:] = 0
  temp[:, :, 1:] = imgData[:, :, :-1]
  inds = numpy.where(numpy.bitwise_and(imgData > 0 , temp == 0))
  boundaryVoxels[inds] = 1
  
  bdryInds = numpy.where(boundaryVoxels > 0) 

  # Create structuring element for erosion. 6-connectivity from central voxel in a 3x3x3 cube.  
  strel = numpy.zeros((3,3,3))
  inds = numpy.asarray( [[0,1,1], [1,0,1], [1,1,0], [2,1,1], [1,2,1], [1,1,2], [1,1,1]] )
  inds = (inds[:,0], inds[:,1], inds[:,2])
  strel[inds] = 1
  
  temp = morphology.binary_erosion(imgData[:,:,:,0], structure=strel)
  erodedData = numpy.zeros(imgData.shape)
  erodedData[numpy.where(temp)] = 1
  
  boundaryVoxels = imgData - erodedData
  bdryInds = numpy.where(boundaryVoxels > 0) 


  ijk1 = numpy.vstack(bdryInds[0:3])
  # Convert to homogeneous coordinates
  ijk1 = numpy.vstack( (ijk1, numpy.ones(bdryInds[3].shape)) )

  for n in range(ijk1.shape[1]):
    x,y,z,w = numpy.dot(aff, ijk1[:, n])
    pts.InsertNextPoint(x,y,z)


  pd.SetPoints(pts)
  pd.Update()
  
  
  cube = vtk.vtkCubeSource()
  cube.SetXLength(pixdims[0])
  cube.SetYLength(pixdims[1])
  cube.SetZLength(pixdims[2])
  
  glyph = vtk.vtkGlyph3D()
  glyph.SetInput(pd)
  glyph.SetSourceConnection(cube.GetOutputPort())
  glyph.ScalingOff()
  glyph.Update()
  
  pdWriter = vtk.vtkPolyDataWriter()
  pdWriter.SetFileName(outputFileName)
  pdWriter.SetInput(glyph.GetOutput())
  pdWriter.Update()
  
 
## Set up the glyph filter
#glyph = vtk.vtkGlyph3D()
#glyph.SetInputConnection(elev.GetOutputPort())
#glyph.SetSourceConnection(sph.GetOutputPort())
#glyph.ScalingOn()
#  

if __name__ == '__main__':
  sys.exit(main(*sys.argv))
