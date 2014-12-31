

# README: THIS FUNCTION IS MORE OF A PRACTICE REALLY, SAREG AND THE POINT
# REGISTRATION IT USES APPEARS TO DO CENTERING ANYWAY.

scriptDir = '/Users/paulaljabar/work/scripts/python'
irtkDir = '/Users/paulaljabar/work/packages/irtk/build/bin'


import sys, os

# for debugging
sys.path.append('/Users/paulaljabar/Python/nibabel-1.3.0-py2.7.egg')
sys.path.append('/Users/paulaljabar/Python')
sys.path.append(scriptDir)

from generalUtilsPythonPA import runCommand

import argparse
import numpy
import vtk 
import math
import tempfile







def main(*args):
  
  helpText = "\
\n\
\n\
\n\
\n\
"
  
  parser = argparse.ArgumentParser(description=helpText, formatter_class=argparse.RawTextHelpFormatter)
  
  helpText = "target surface, e.g. tgt.vtk\n\
"
  parser.add_argument("tgtFile", type=str, help=helpText)

  helpText = "source surface, e.g. src.vtk\n\
"
  parser.add_argument("srcFile", type=str, help=helpText)
  
  helpText = "output transformation, e.g. src-tgt.dof\n\
A dof format transformation for mapping target locations \n\
to source locations"
  parser.add_argument("outputDof", type=str, help=helpText)
  

  #  helpText = "Optional argument"
#  parser.add_argument("-opt", type=int, nargs='+', help=helpText, metavar='label')
#  helpText = "voxel size for resampled output (default 1mm)"
#  parser.add_argument("-Xsize", type=float, nargs='?', help=helpText, default=1.0)


  args = parser.parse_args()

  #############################################################
      
  numpy.set_printoptions(2, suppress=True)


  
  print 'Reading files\n'
  
  tgtFile = args.tgtFile
  srcFile = args.srcFile
  
  if not ( tgtFile.endswith('.vtk') ):
    print "Input files must be vtk files containing landmarks."
    return 1 
  
  tgtReader = vtk.vtkPolyDataReader()
  tgtReader.SetFileName(tgtFile)
  target = vtk.vtkPolyData()
  target = tgtReader.GetOutput()
  target.Update()
    
  srcReader = vtk.vtkPolyDataReader()
  srcReader.SetFileName(srcFile)
  source = vtk.vtkPolyData()
  source = srcReader.GetOutput()
  source.Update()
  
  
  # Target centre of gravity
  nPtsTgt = target.GetNumberOfPoints()
  
  if nPtsTgt < 1:
    print 'No target points!'
    return 1
  
  nPtsSrc = source.GetNumberOfPoints()
  
  if nPtsSrc < 1:
    print 'No source points!'
    return 1
  
  cogTgt = numpy.zeros( (3,) )
  for i in range(nPtsTgt):
    pt = numpy.asarray(target.GetPoint(i))
    cogTgt = cogTgt + pt
    
  cogTgt = cogTgt / nPtsTgt
  print 'Target centre of gravity: ', cogTgt


  print 'Centering data'

  ptsTgtCentred = vtk.vtkPoints()
  ptsSrcCentred = vtk.vtkPoints()
  
  for i in range(nPtsTgt):
    pt = numpy.asarray(target.GetPoint(i))
    pt = pt - cogTgt
    ptsTgtCentred.InsertNextPoint(pt[0], pt[1], pt[2])

  for i in range(nPtsSrc):
    pt = numpy.asarray(source.GetPoint(i))
    pt = pt - cogTgt
    ptsSrcCentred.InsertNextPoint(pt[0], pt[1], pt[2])
    
  target.SetPoints(ptsTgtCentred)
  source.SetPoints(ptsSrcCentred)
  
  tempDir = tempfile.gettempdir()
  pid = str(os.getpid())
  tgtTempFile = os.path.join(tempDir , 'tmp-tgt-' + pid + '.vtk' )
  srcTempFile = os.path.join(tempDir , 'tmp-src-' + pid + '.vtk' )
  dofTempFile = os.path.join(tempDir , 'tmp-src-tgt-' + pid + '.dof' )
  
  
  tgtWriter = vtk.vtkPolyDataWriter()
  tgtWriter.SetInput(target)
  tgtWriter.SetFileName(tgtTempFile)
  tgtWriter.Update()
  
  srcWriter = vtk.vtkPolyDataWriter()
  srcWriter.SetInput(source)
  srcWriter.SetFileName(srcTempFile)
  srcWriter.Update()
  
  cmd = irtkDir + '/sareg'
  cmd = cmd + ' ' + tgtTempFile
  cmd = cmd + ' ' + srcTempFile
  cmd = cmd + ' ' + '-dofout ' + dofTempFile
  cmd = cmd + ' ' + '-epsilon 0.0001 -symmetric'
  
  _,_ = runCommand(cmd)

#  cmd = irtkDir + '/stransformation'
#  cmd = cmd + ' ' + workingDir + '/cortex-mid-40-weeks.vtk'
#  cmd = cmd + ' ' + workingDir + '/temp.vtk'
#  cmd = cmd + ' ' + '-dofin ' + dofFileT2toTemplate
#  _,_ = runCommand(cmd)


  # TODO: Remove temp files
  
  ##########################################################


  print 'done'

  
  

if __name__ == '__main__':
  sys.exit(main(*sys.argv))


