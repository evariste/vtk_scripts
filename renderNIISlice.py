import os
import sys


sys.path.append('/Users/paulaljabar/Python/nibabel-1.3.0-py2.7.egg')

import nibabel as nib
import numpy
import vtk

# Extracts a x-y slice for a given z-slice index, saves the result in a
# rectangular polydata object with scalar values containing the pixel values.

def usage(scriptName):
  print "Usage: "
  print "   python " + scriptName + " [input.nii.gz] [output.vtk] <slice>"
  print ""
  print "  Generate a rendering of a slice of a nifti image."
  print "  Output written in a vtk file."
  print "  Can set optional slice to render, default is first (index 0)"
  print " "
  exit()


# TODO: extend so that a x-z or a y-z slice can also be extracted.
def main(*args):
  if len(args) < 3:
     usage(args[0])
  
  inputFilename = args[1]
  outputFilename = args[2]

  print len(args)
  print args[0]

  if len(args) < 4:
    zSlice = 0
  else:
    zSlice = int(args[3])
    
  print 'renderNIISlice: input file  %s', (inputFilename)
  print '                output file %s', (outputFilename)

  img = nib.load(inputFilename)
  
  
  imgData = img.get_data()
  # (xdim, ydim, zdim, tdim) = img.shape
  dims = img.shape
  xdim = dims[0]
  ydim = dims[1]
  zdim = dims[2]
  
  if zSlice > zdim - 1:
    print 'Warning: chosen slice exceeds number available'
    print '         setting to slice index %d' % (zdim - 1)       
    zSlice = zdim - 1
    
  if zSlice < 0:
    print 'Warning: chosen slice negative, setting to slice index 0'
    zSlice = 0
  
  # Storage for pixel values.
  sliceScalars = vtk.vtkFloatArray()
  sliceScalars.SetNumberOfComponents(1)
  sliceScalars.SetNumberOfValues(xdim*ydim)
  
  
  count = 0
  maxVal = 0
  minVal = 1000000
  
  for j in range(ydim):
    for i in range(xdim):
      val = imgData[i, j, zSlice]    
      sliceScalars.InsertTuple1(count, val)    
      if maxVal < val:
        maxVal = val
      if minVal > val:
        minVal = val
      
      count = count + 1
  
  print '%d voxels, min and max vals %d to %d' % (count, minVal, maxVal)    

  ## Needed if writing to a jpeg image  
  #for n in range(xdim*ydim):
  #  val = sliceScalars.GetValue(n)
  #  val = 127 * (val - minVal) / maxVal
  #  sliceScalars.SetValue(n, val)
  
  
  # Make a tiled rectangle to fit the geometry of the chosen slice.
  
  plane = vtk.vtkPlaneSource()
  i2wMatrix = img.get_affine()
  # Set plane corners to match world coordinate corners of slice
  v = [ [0], [0], [zSlice], [1] ]
  [ [x], [y], [z], [w] ] = numpy.dot(i2wMatrix, v)
  plane.SetOrigin(x, y, z)
  
  v = [ [xdim-1], [0], [zSlice], [1] ]
  [ [x], [y], [z], [w] ] = numpy.dot(i2wMatrix, v)
  plane.SetPoint1(x, y, z)
  
  v = [ [0], [ydim-1], [zSlice], [1] ]
  [ [x], [y], [z], [w] ] = numpy.dot(i2wMatrix, v)
  plane.SetPoint2(x, y, z)
  
  # Number of tiles based on the number of intervals between the voxel centres
  # along each dimension.
  plane.SetXResolution(xdim-1)
  plane.SetYResolution(ydim-1)
  
  planePolyData = plane.GetOutput()
  planePolyData.Update()
  
  planePolyData.GetPointData().SetScalars(sliceScalars)
  
  triFilter = vtk.vtkTriangleFilter()
  triFilter.SetInput(planePolyData)
  
  polydataWriter = vtk.vtkPolyDataWriter()
  polydataWriter.SetInput(triFilter.GetOutput())
  polydataWriter.SetFileName(outputFilename)
  polydataWriter.Write()
  
  

if __name__ == '__main__':
  sys.exit(main(*sys.argv))


  ## Visualisation within vtkpython
  
#  # Image data for texture
#  texImgData = vtk.vtkImageData()
#  texImgData.SetExtent(0, xdim-1, 0, ydim-1, 0, 0)
#  texImgData.AllocateScalars()
#  texImgData.GetPointData().SetScalars(sliceScalars)
#  texImgData.Update()
#  
#  
#  # Texture object
#  texture = vtk.vtkTexture()
#  texture.SetInput(texImgData)

  #texturePlane = vtk.vtkTextureMapToPlane()
  #texturePlane.SetInput(plane.GetOutput())
  #
  #planeMapper = vtk.vtkPolyDataMapper()
  #planeMapper.SetInput(texturePlane.GetOutput())
  # 
  #texturedPlane = vtk.vtkActor()
  #texturedPlane.SetMapper(planeMapper)
  #texturedPlane.SetTexture(texture)
  # 
  ## Visualize the textured plane
  #renderer = vtk.vtkRenderer()
  #renderer.AddActor(texturedPlane)
  #renderer.SetBackground(.9, .9, .9)
  #renderer.ResetCamera()
  # 
  #renderWindow = vtk.vtkRenderWindow()
  #renderWindow.AddRenderer(renderer)
  #
  #renderWindowInteractor = vtk.vtkRenderWindowInteractor()
  #renderWindowInteractor.SetRenderWindow(renderWindow)
  # 
  #renderWindow.Render()
  # 
  #renderWindowInteractor.Start()
  
  
  ### Writing to an image
  #castFilter = vtk.vtkImageCast()
  #castFilter.SetOutputScalarTypeToUnsignedChar()
  #castFilter.SetInput(texImgData)
  #castFilter.Update()
  #
  #jpgWriter = vtk.vtkJPEGWriter()
  #jpgWriter.SetFileName('bla.jpg')
  #jpgWriter.SetInput(castFilter.GetOutput())
  #jpgWriter.Write()
  

