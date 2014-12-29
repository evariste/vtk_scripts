


import sys, os

#from subprocess import STDOUT
# for debugging
sys.path.append('/Users/paulaljabar/Python/')
#import dicom
#import nibabel as nib
import argparse
import numpy as np
import subprocess
#import glob
import vtk 
import math
import shutil
import tempfile



scriptDir = '/Users/paulaljabar/work/scripts/python'
irtkDir = '/Users/paulaljabar/work/packages/irtk/build/bin'


def main(*args):
  
  helpText = "Run point-based rigid registration between a pair of \
landmark sets. These are assumed to contain in order: \
   top of stem, \
   bottom  of stem, \
   one eye, \
   and the other eye. RUN WITH VTKPYTHON."
  
  parser = argparse.ArgumentParser(description=helpText)
  
  helpText = "in: landmarksA.vtk  : vtk file with landmarks A"
  parser.add_argument("landmarksFileA", type=str, help=helpText)
  
  helpText = "in: landmarksB.vtk  : vtk file with landmarks B"
  parser.add_argument("landmarksFileB", type=str, help=helpText)

  helpText = "out: prreg.dof : Transformation from landmarks A  to landmarks B "
  parser.add_argument("dofout", type=str, help=helpText)
    
#  helpText = "Optional argument"
#  parser.add_argument("-opt", type=int, nargs='+', help=helpText, metavar='label')

  #############################################################
      
  args = parser.parse_args()
  
  fileA = args.landmarksFileA
  fileB = args.landmarksFileB
  
  if not ( fileA.endswith('.vtk')  and fileB.endswith('.vtk') ):
    print "Input files must be vtk files containing landmarks."
    return 1 
  
  
  
  readerA = vtk.vtkPolyDataReader()
  readerA.SetFileName(fileA)
  pdA = vtk.vtkPolyData()
  pdA = readerA.GetOutput()
  pdA.Update()
  
  readerB = vtk.vtkPolyDataReader()
  readerB.SetFileName(fileB)
  pdB = vtk.vtkPolyData()
  pdB = readerB.GetOutput()
  pdB.Update()
  
  stemATop = np.asarray(pdA.GetPoint(0))
  stemABottom = np.asarray(pdA.GetPoint(1))
  eyeA1 = np.asarray(pdA.GetPoint(2))
  eyeA2 = np.asarray(pdA.GetPoint(3))

  stemBTop = np.asarray(pdB.GetPoint(0))
  stemBBottom = np.asarray(pdB.GetPoint(1))
  eyeB1 = np.asarray(pdB.GetPoint(2))
  eyeB2 = np.asarray(pdB.GetPoint(3))
  
  eyeAMid = 0.5 * (eyeA1 + eyeA2)
  eyeBMid = 0.5 * (eyeB1 + eyeB2)
  
  # Vector from bottom to top of stem
  stemAVec = stemATop - stemABottom
  stemBVec = stemBTop - stemBBottom
  
  # Vector from top of stem to point mid-way between eyes.
  stem2frontVecA = eyeAMid - stemATop 
  stem2frontVecB = eyeBMid - stemBTop
  
  # Vector inter-eyes
  eyeAVec = eyeA1 - eyeA2 
  eyeBVec = eyeB1 - eyeB2
  
  # Vector orthogonal to stem vector and stem-mid-eye vector.
  vA = np.cross(stemAVec, stem2frontVecA)
  vB = np.cross(stemBVec, stem2frontVecB)

  # The angle between the last vector and the inter-eye vector should be close to 0 or 180.
  # I.e. the cosine shoud be close to 1 or -1 
  cosA = np.dot(vA, eyeAVec) / np.sqrt ( np.dot(vA,vA) * np.dot(eyeAVec, eyeAVec) )
  cosB = np.dot(vB, eyeBVec) / np.sqrt ( np.dot(vB,vB) * np.dot(eyeBVec, eyeBVec) ) 

  angleA = math.acos( cosA )
  angleB = math.acos( cosB )
  
  if (angleA > math.pi/4) and (angleA < 3*math.pi/4):
    print "Warning vector between eyes is far from orthogonal to stem vector for landmarks ", fileA
    
  if (angleB > math.pi/4) and (angleB < 3*math.pi/4):
    print "Warning vector between eyes is far from orthogonal to stem vector for landmarks ", fileB
    
  targetFileName = fileA
  sourceFileName = fileB
  
  if cosA * cosB < 0:
    print "Detected opposite senses for eyes, will flip points for eyes in first landmarks"
    pdA.GetPoints().SetPoint(2, eyeA2[0], eyeA2[1], eyeA2[2])
    pdA.GetPoints().SetPoint(3, eyeA1[0], eyeA1[1], eyeA1[2])
    # Strip '.vtk' and put '-flip.vtk'
    targetFileName = targetFileName[:-4] + '-flip.vtk'
    writer = vtk.vtkPolyDataWriter()
    writer.SetInput(pdA)
    writer.SetFileName(targetFileName)
    writer.Write()
    




  
  prregExe = os.path.join(irtkDir, 'prreg')
  
   
  cmd = [prregExe, targetFileName, sourceFileName, '-dofout', args.dofout]

  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  (stdOut, stdErr) = p.communicate()
  if not p.returncode == 0:
    print "Error running prreg. Command was"
    print cmd
    print "Error provided:"
    print stdErr
    return 1
    
  
      
#  prreg dyn_tra-01.vtk  dyn_tra2-01.vtk  -dofout dyn_tra2-01-dyn_tra-01.dof
# 1033  07:39  rview dyn_tra-01.nii.gz  dyn_tra2-01.nii.gz -dofin dyn_tra2-01-dyn_tra-01.dof
  
  print stdOut

  
  
  

if __name__ == '__main__':
  sys.exit(main(*sys.argv))




#import vtk
#
#
#filename='temp.wrl'
#
#rend = vtk.vtkRenderer()
#renWin = vtk.vtkRenderWindow()
#renWin.AddRenderer(rend)
#
#renWinInt = vtk.vtkRenderWindowInteractor()
#renWinInt.SetRenderWindow(renWin)
#
#importer = vtk.vtkVRMLImporter()
#importer.SetFileName( filename )
#importer.Read()
#
#importer.SetRenderWindow(renWin)
#
#print dir(importer)
#importer.Update()
#
#renWin.Render()
#renWinInt.Start()
#
#
