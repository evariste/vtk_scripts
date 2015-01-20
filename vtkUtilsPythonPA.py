
import vtk

def vtkMatrixFromNumpyArray4x4(npArray):
  """
  Return a VTK matrix from a numpy 4x4 array.
  """
  if (npArray.shape != (4,4)):
    print "Error: numpy array is not 4x4, cannot convert to VTK matrix"
    return
       
  vtkMat = vtk.vtkMatrix4x4()
       
  for row in range(npArray.shape[0]):
    for col in range(npArray.shape[1]):
      vtkMat.SetElement(row,col,npArray[row][col])
    
  return vtkMat