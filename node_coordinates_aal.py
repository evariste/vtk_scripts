
import sys
sys.path.append('/Users/paulaljabar/work/scripts/python/pythonEnvs/surfaceAnalysis/lib/python2.7/site-packages')

import vtk
import nibabel as nib
import numpy as np




atlasDir = '/Users/paulaljabar/work/atlases/aal-gareth'

atlasFile = atlasDir + '/tr-aal-to-neo-nr-40-all.nii.gz'

outputFile = atlasDir + '/tr-aal-to-neo-nr-40-all-label-centres.vtk'

img = nib.load(atlasFile)


data = img.get_data()
aff = img.get_affine()

labels = np.unique(data[data > 0])

pts = vtk.vtkPoints()
  
pts.SetNumberOfPoints(labels.size)
  
n = 0


for label in labels:
  i, j, k, _ = np.where(data == label)
  m_i = np.mean(i)
  m_j = np.mean(j)
  m_k = np.mean(k)
  v = np.asarray([ m_i, m_j, m_k, 1])
  x,y,z,_ = aff.dot(v)
  pts.InsertPoint(n, x, y, z)
  n = n + 1


pd = vtk.vtkPolyData()
pd.SetPoints(pts)


writer = vtk.vtkPolyDataWriter()
writer.SetInput(pd)
writer.SetFileName(outputFile)
writer.Write()



print 'hello'







