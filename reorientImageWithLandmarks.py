


import sys, os
from numpy.ma.core import fabs

#from subprocess import STDOUT
# for debugging

sys.path.append('/Users/paulaljabar/Python/nibabel-1.3.0-py2.7.egg')
sys.path.append('/Users/paulaljabar/Python')

import geometryUtilsPA as geompa
import imageUtilsPythonPA as impa

import nibabel as nib
import argparse
import numpy as np
import vtk 
import math
import tempfile
import subprocess
#import dicom
#import glob
#import shutil



scriptDir = '/Users/paulaljabar/work/scripts/python'
irtkDir = '/Users/paulaljabar/work/packages/irtk/build/bin'


def normalise(vec):
  m = magnitude(vec)
  if m < 0.000001:
    print "Warning: normalise: vector magnitude very small"
  return vec / m


def magnitude(vec):
  return math.sqrt( np.dot(vec, vec) )


def runCommand(cmd):
  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  (stdOut, stdErr) = p.communicate()
  if not p.returncode == 0:
    print "Error running command"
    print cmd
    print "Error provided:"
    print stdErr
    return 1




def main(*args):
  
  helpText = "\
Re-orient then resample an MR image. The re-orientation relies \n\
on  using a set of landmark points given in a vtk file. \n\
\n\
The landmark file is assumed to contain the following points in order: \n\
   top of stem, \n\
   bottom  of stem, \n\
   one eye, \n\
   and the other eye. \n\
\n\
This script needs to be run with the vtk version of python (vtkpython)."
  
  parser = argparse.ArgumentParser(description=helpText, formatter_class=argparse.RawTextHelpFormatter)
  
  helpText = "input image, e.g. image-in.nii.gz\n\
MR image to re-orient (full path)"
  parser.add_argument("inputImage", type=str, help=helpText)

  helpText = "input landmark points, e.g. landmarks.vtk\n\
A vtk file with landmarks (full path)"
  parser.add_argument("landmarksFile", type=str, help=helpText)
  
  helpText = "output image, e.g. image-out.nii.gz\n\
The re-oriented and resampled MR image"
  parser.add_argument("outputImage", type=str, help=helpText)
  
  helpText = "output transformation, e.g. original-new.dof\n\
A dof format transformation for which the new image is \n\
the target and the original image is the source"
  parser.add_argument("outputDof", type=str, help=helpText)
  

  helpText = "voxel size for resampled output (default 1mm)"
  parser.add_argument("-Xsize", type=float, nargs='?', help=helpText, default=1.0)
  helpText = "voxel size for resampled output (default 1mm)"
  parser.add_argument("-Ysize", type=float, nargs='?', help=helpText, default=1.0)
  helpText = "voxel size for resampled output (default 1mm)"
  parser.add_argument("-Zsize", type=float, nargs='?', help=helpText, default=1.0)
  
  helpText = "Tilt angle (degrees): Controls direction of stem (default 25 degrees)"
  parser.add_argument("-tiltAngle", type=float, nargs='?', help=helpText, default=25.0)
  
#  helpText = "Optional argument"
#  parser.add_argument("-opt", type=int, nargs='+', help=helpText, metavar='label')


  args = parser.parse_args()

  #############################################################
      
  np.set_printoptions(2, suppress=True)


  
  # READ LANDMARKS
  print 'Reading landmarks\n'
  
  landmarksFile = args.landmarksFile
  
  if not ( landmarksFile.endswith('.vtk') ):
    print "Input files must be vtk files containing landmarks."
    return 1 
  
  pd_reader = vtk.vtkPolyDataReader()
  pd_reader.SetFileName(landmarksFile)
  pd = vtk.vtkPolyData()
  pd = pd_reader.GetOutput()
  pd.Update()
    

  
  # START CHECKING CONFIGURATION (basic checks)
  print 'doing basic checks on landmark configuration\n'
  
  stemTop = np.asarray(pd.GetPoint(0))
  stemBottom = np.asarray(pd.GetPoint(1))
  eyeA = np.asarray(pd.GetPoint(2))
  eyeB = np.asarray(pd.GetPoint(3))
  
  eyeMid = 0.5 * (eyeA + eyeB)
  
  # Vector from bottom to top of stem
  stemVec = stemTop - stemBottom
  
  # Vector from top of stem to point mid-way between eyes.
  stemToEyeMidVec = eyeMid - stemTop 
  
  # Vector inter-eyes
  eyeVec = eyeA - eyeB 
  
  # Vector orthogonal to stem vector and stem-mid-eye vector.
  estdLateralVec = np.cross(stemVec, stemToEyeMidVec)

  # The angle between the estimated lateral vector and the inter-eye vector
  # should be close to 0 or 180. I.e. the cosine shoud be close to 1 or -1
  cosAng = np.dot(estdLateralVec, eyeVec)
  cosAng = cosAng / magnitude(estdLateralVec) / magnitude(eyeVec)
  
  angleLatAndEyeVec = math.acos( cosAng )
  
  if (angleLatAndEyeVec > math.pi/4) and (angleLatAndEyeVec < 3*math.pi/4):
    print "Warning vector between eyes is far from orthogonal to stem vector for landmarks ", landmarksFile


  # END CHECKING .....
  
  
  
  
#  headertool  ../misc_controls/5674470N/recon-5674470N.nii.gz img-reset.nii.gz  -origin 0 0 0 -orientation 1 0 0  0 1 0 0 0 1

  ##########################################################
  print 'Resetting the header of input image to canonical coordinates\n'
  
  imgIn = nib.load(args.inputImage)
  i2w = imgIn.get_affine()    
  w2i = np.linalg.inv(i2w)
  
  hdr = imgIn.get_header()
  pixdims = hdr['pixdim']
    
  S = np.eye(4,4)
  S[0,0] = pixdims[1]
  S[1,1] = pixdims[2]
  S[2,2] = pixdims[3]
  
  dims = hdr['dim']
  T = np.eye(4)
  T[0,3] = -1*pixdims[1]*(dims[1] - 1)/2.0
  T[1,3] = -1*pixdims[2]*(dims[2] - 1)/2.0
  T[2,3] = -1*pixdims[3]*(dims[3] - 1)/2.0
  
  i2w_reset = T.dot(S)

  imgOut = nib.Nifti1Image(imgIn.get_data(), i2w_reset)
  resetImage = os.path.join( tempfile.gettempdir() , 'tmp-reset' + str(os.getpid()) + '.nii.gz' )
  nib.save(imgOut, resetImage)


  # Affine transformation from original image to canonical one  
  matOriginal2ResetImage = i2w_reset.dot(w2i)


  
  # Store all the points in an array. This set will be iteratively modified, so
  # we currently have iteration 1
  allPts1 = np.ones( (4,4) )
  allPts1[0:3, 0] = stemTop
  allPts1[0:3,1]  = stemBottom
  allPts1[0:3,2]  = eyeA
  allPts1[0:3,3]  = eyeB


  # Transform all points to frame of canonical image
  allPts1 = matOriginal2ResetImage.dot(allPts1)
  

  # All following operations take place in frame of canonical ('reset') image.
  
  ##########################################################
  print 'Estimating transformation to make eye to eye vector along the x axis\n'
  
  # Get eye to eye vector to rotate into xy plane
  eyeVec = allPts1[:,2] - allPts1[:,3]
  v = np.copy(eyeVec)
  v[2] = 0
  theta = geompa.angleBetweenVectors(eyeVec, v)
  ax = np.cross(eyeVec[0:3], v[0:3])
  
  if np.array_equal(ax, [0,0,0]):
    R = np.eye(3)
  else:
    R = geompa.rotationGivenAxisAndAngle(ax, theta)
  
  R1 = np.eye(4)
  R1[0:3,0:3] = R
    

  allPts2 = np.dot(R1, allPts1)
  
  # Get current eye to eye vector to rotate onto x axis
  eyeVec2 = allPts2[:,2] - allPts2[:,3]

  v = np.zeros( (4,) )
  v[0] = 1
  theta = geompa.angleBetweenVectors(eyeVec2, v)
  
  if theta > np.pi / 2.0:
    v[0] = -1
    theta = np.pi - theta

  ax = np.cross(eyeVec2[0:3], v[0:3])
  
  if np.array_equal(ax, [0,0,0]):
    R = np.eye(3)
  else:
    R = geompa.rotationGivenAxisAndAngle(ax, theta)
  
  R2 = np.eye(4)
  R2[0:3,0:3] = R
  
  allPts3 = np.dot(R2, allPts2)


  
  ##########################################################
  print 'Now trying to make stem vector as close to vertical as possible (incorporating tilt angle in direction of +y).\n'
  
  # Now carry out a rotation about the current vector between the eyes. The
  # rotation should make the stem vector from bottom to top point as much
  # possible in the 'up' z direction and a little 'forward', in the direction of
  # the y-axis

  tiltAngle = -1*args.tiltAngle # Defaults to -25.0 # degrees
  tiltAngle = tiltAngle * np.pi / 180
  v = np.zeros( (4,) )
  v[1] = math.sin(tiltAngle)
  v[2] = math.cos(tiltAngle)
  
  eyeVec3  = allPts3[:,2] - allPts3[:,3]
  stemVec3 = allPts3[:,0] - allPts3[:,1]
  
  nThetas = 120
  thetas = np.linspace(0, np.pi * 2, nThetas, endpoint=False)
  scores = np.zeros( thetas.shape )
  
  R3 = np.eye(4)

  for j in range(len(thetas)):
    theta = thetas[j]
    R = geompa.rotationGivenAxisAndAngle(eyeVec3[0:3], theta)
    R3[0:3,0:3] = R
    temp = np.dot(R3, stemVec3)
    scores[j] = np.dot(temp, v)
    

  j = np.argmax(scores)
  theta = thetas[j]
  R = geompa.rotationGivenAxisAndAngle(eyeVec3[0:3], theta)
  R3[0:3,0:3] = R
  
  allPts4 = np.dot(R3, allPts3)
  
  
  ##########################################################

  print 'Checking that the eyes are anterior to the stem, perform 180 rotation if not\n'

  # We would like the P-A direction to represent increasing y coordinates. So,
  # if the y-coords for the eyes have ended up less than those for the stem
  # apply a rotation of 180 around the z axis.

  # Eye y-coordinate, should be the same for both eyes.
  y1 = allPts4[1,2]
  # Mean y coordinate of bottom and top of stem.
  y2 = 0.5*(allPts4[1,0] + allPts4[1,1])

  R4 = np.eye(4)

  if y1 < y2:
    R4[0,0] = -1.0
    R4[1,1] = -1.0


  allPts5 = np.dot(R4, allPts4)
  

  
  ##########################################################
  print 'Done finding estimated transformation for points\n'
  
  # The summary rotation that would map the original set of points to the
  # current position
 
  R = R4.dot(R3).dot(R2).dot(R1)

  # Recalculate total rotation in one step to reduce numeric inaccuracy.
  allPts5 = R.dot(allPts1)
  
  
  
  
  
  
  ##########################################################
  print 'Generating dof to map image to new grid'
  tempMatFile = os.path.join( tempfile.gettempdir() , 'tmp-' + str(os.getpid()) + '.mat' )
  tempDofFile = os.path.join( tempfile.gettempdir() , 'tmp-' + str(os.getpid()) + '.dof' )
  
  impa.writeIRTKMatrix(tempMatFile, np.linalg.inv(R))

  mat2dofExe = os.path.join(irtkDir, 'mat2dof')  
  cmd = [mat2dofExe, tempMatFile, tempDofFile]

  runCommand(cmd)
  print




  ##########################################################
  xSz = args.Xsize
  ySz = args.Ysize
  zSz = args.Zsize
  print 'Upsample the reset image to a resolution in x-y-z of {0}, {1}, {2}'.format(xSz, ySz, zSz)
  
  resampledImage = os.path.join( tempfile.gettempdir() , 'tmp-upsample-' + str(os.getpid()) + '.nii.gz' )
  
  resampleExe = os.path.join(irtkDir, 'resample')
  cmd = [resampleExe, resetImage, resampledImage]
  cmd = cmd + ['-size']
  cmd = cmd + map(str, [xSz, ySz, zSz])

  runCommand(cmd)
  print


  ##########################################################
  print 'Save reoriented image data onto upsampled reset grid'
  transformExe = os.path.join(irtkDir, 'transformation')
  cmd = [transformExe, resetImage, args.outputImage]
  cmd = cmd + ['-target', resampledImage]
  cmd = cmd + ['-dofin', tempDofFile]
  cmd = cmd + ['-bspline']

  runCommand(cmd)
  print
  
  
  ##########################################################
  print 'Save concatenated transformation from original data to reset image to reoriented image'


  outputDofFile = args.outputDof
  
  matReset2Original = np.linalg.inv(matOriginal2ResetImage)
  invR = np.linalg.inv(R)
  myMat = matReset2Original.dot(invR)
  impa.writeIRTKMatrix(tempMatFile, myMat)
  mat2dofExe = os.path.join(irtkDir, 'mat2dof')

  
  cmd = [mat2dofExe, tempMatFile, outputDofFile]
  runCommand(cmd)
  print


  # TODO: Save landmarks in space of reoriented upsampled image data?


  ##########################################################
  print 'Remove temporary files.'
  fileList = [resetImage,  resampledImage, tempMatFile, tempDofFile]
  for f in fileList:
    os.remove(f)


  print 'done'


  
  
  
  
#  temp = np.zeros( (4,) )
#  temp[0:3] = eyeVec
#  temp = np.dot(w2i, temp)
#  
#  u1 = normalise(temp)
#  u1 = u1[0:3]
#
#
#  temp = np.zeros( (4,) )
#  temp[0:3] = stemToEyeMidVec
#  temp = np.dot(w2i, temp)  
#  temp = temp[0:3]
#  
#  u3 = np.cross(u1, temp)
#  u3 = normalise(u3)
#  
#  u2 = np.cross(u3, u1)
#  u2 = normalise(u2)
#
#
#  
#  R = np.eye(4,4)
#  R[0:3,0] = u1
#  R[0:3,1] = u2
#  R[0:3,2] = u3
#  
#  R = R.T
#  
#  hdr = imgIn.get_header()
#  pixdims = hdr['pixdim']
#  
#  
#  S = np.eye(4,4)
#  S[0,0] = pixdims[1]
#  S[1,1] = pixdims[2]
#  S[2,2] = pixdims[3]
#
#  
#  T = np.eye(4,4)
##  T[0:3,3] = stemTop
#  
#  M = np.dot(T, R)
#  M = np.dot(M, S)
#  imgOut = nib.Nifti1Image(imgIn.get_data(), M)
#  nib.save(imgOut, args.outputImage)



#
#
#
#  
#  prregExe = os.path.join(irtkDir, 'prreg')
#  
#   
#  cmd = [prregExe, targetFileName, sourceFileName, '-dofout', args.dofout]
#
#  p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#  (stdOut, stdErr) = p.communicate()
#  if not p.returncode == 0:
#    print "Error running prreg. Command was"
#    print cmd
#    print "Error provided:"
#    print stdErr
#    return 1
#    
#  
      
#  prreg dyn_tra-01.vtk  dyn_tra2-01.vtk  -dofout dyn_tra2-01-dyn_tra-01.dof
# 1033  07:39  rview dyn_tra-01.nii.gz  dyn_tra2-01.nii.gz -dofin dyn_tra2-01-dyn_tra-01.dof
  
#  print stdOut

  
  
  

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
