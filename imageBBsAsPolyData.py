import sys, os

# for debugging
sys.path.append('/Users/paulaljabar/Python/nibabel-1.3.0-py2.7.egg')

import nibabel as nib
import argparse
import numpy as np

import vtk
from vtk.util import numpy_support

from vtkUtilsPythonPA import vtkMatrixFromNumpyArray4x4


scriptDir = '/Users/paulaljabar/work/scripts/python'
irtkDir = '/Users/paulaljabar/work/packages/irtk/build/bin'


def main(*args):
  
  helpText = "    Get the bounding box(es) of one or more images fields of view and store in a vtk polydata set"
  
  parser = argparse.ArgumentParser(description=helpText)
  
  helpText = "out: output.vtk"
  parser.add_argument("outputBoxes", type=str, help=helpText)
  
    
  helpText = "images to find ROIs for (nifti format)"
  parser.add_argument("images", type=str, nargs='+', help=helpText, metavar='label')

  #############################################################
      
  args = parser.parse_args()
  
  output = args.outputBoxes
  input = args.images

  if len(input) == 0:
    print 'No images given, exiting'
    sys.exit(1)
    
  # Do all the input files exist?
  fileCheck = map(lambda f: os.path.isfile(f), input)
  allExist = reduce(lambda x,y: x and y, fileCheck)
  if not allExist:
    print 'Not all input files exist, exiting'
    sys.exit(1)
    

    
  append = vtk.vtkAppendPolyData()

  for f in input:
    img = nib.load(f)
    
    dimx, dimy, dimz = img.get_header()['dim'][1:4]

    outline = vtk.vtkOutlineSource()
    outline.SetGenerateFaces(1)
    outline.SetBounds(0,dimx-1,0,dimy-1,0,dimz-1)
    outline.Update()
    pdBox = outline.GetOutput()
    pdBox.Update()
    
    
    m = img.get_affine()
    m = vtkMatrixFromNumpyArray4x4(m)
    print m

    T = vtk.vtkTransform()
    T.SetMatrix(m)
    T.Update()
    
    tFilt = vtk.vtkTransformFilter()
    tFilt.SetInput(pdBox)
    tFilt.SetTransform(T)
    tFilt.Update()

    pdTr = tFilt.GetOutput()
    pdTr.Update()
    
    append.AddInput(pdTr)
    append.Update()


    print 'processed file ', f
    
  print 'writing output to ', output
  
  writer = vtk.vtkPolyDataWriter()
  writer.SetFileName(output)
  writer.SetInput(append.GetOutput())
  writer.Write()

if __name__ == '__main__':
  sys.exit(main(*sys.argv))


