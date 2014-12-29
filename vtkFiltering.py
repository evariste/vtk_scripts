
# fslstats wm.nii.gz -w
# region wm.nii.gz  wm-crop.nii.gz -Rx1 75 -Rx2 181 -Ry1 59 -Ry2 199 -Rz1 11 -Rz2 93 
# resample_real wm-dilate.nii.gz wm-dilate-resamp.nii.gz  -size 6.875 6.875 8.0808 

# makesequence wm-crop-resamp.nii.gz wm-crop-resamp.nii.gz  wm-crop-resamp.nii.gz wm-crop-resamp.nii.gz wm-crop-resamp.nii.gz wm-crop-resamp.nii.gz wm-crop-resamp-4D.nii.gz
# padding data_tensor_nl-crop-resamp.nii.gz wm-crop-resamp-4D.nii.gz data_tensor_nl-crop-resamp-pad.nii.gz 0 0 

import vtk
import nibabel as nib
import numpy






segImg = nib.load('t2_seg.nii.gz')

segData = segImg.get_data()


tensorImg = nib.load('data_tensor_nl-crop-resamp-pad.nii.gz')

tensorData = tensorImg.get_data()

print segData.shape

print tensorData.shape

sg = vtk.vtkStructuredGrid()
sg.SetDimensions(tensorData.shape[:3])

pts = vtk.vtkPoints()
nPts = numpy.prod(tensorData.shape[:3])
print 'nPts: ' , nPts

aff = tensorImg.get_affine()

xdim = tensorData.shape[0]
ydim = tensorData.shape[1]
zdim = tensorData.shape[2]

ijk = [[i,j,k,1] for k in range(zdim)  for j in range(ydim) for i in range(xdim)]

tDataOut = vtk.vtkFloatArray()
tDataOut.SetNumberOfComponents(9)



for imCoord in ijk:
  wCoord = numpy.dot(aff, numpy.array(imCoord))
  pts.InsertNextPoint(wCoord[0], wCoord[1], wCoord[2])
  tDataIn = tensorData[imCoord[0], imCoord[1], imCoord[2], :]
  tDataIn = (tDataIn[0], tDataIn[1], tDataIn[2], \
             tDataIn[1], tDataIn[3], tDataIn[4], \
             tDataIn[2], tDataIn[4], tDataIn[5])
  tDataOut.InsertNextTupleValue(tDataIn)
  
tDataOut.Squeeze()
  
sg.SetPoints(pts)

sg.GetPointData().SetTensors(tDataOut)

sgWriter = vtk.vtkStructuredGridWriter()
sgWriter.SetInput(sg)
sgWriter.SetFileName('bla.vtk')
sgWriter.Update()



## run with vtkpython after setting python path etc using source script.
#
#from vtk import *
#
#filename = "temp.vtk"
#
#reader = vtkStructuredGridReader()
#reader.SetFileName(filename)
#reader.Update()
#
#grid = vtkStructuredGrid()
#grid = reader.GetOutput()
#grid.Update()
#
#nPoints = grid.GetNumberOfPoints()
#nCells  = grid.GetNumberOfCells()
#dims    = grid.GetDimensions()
#
#ni = dims[0]
#nj = dims[1]
#nk = dims[2]
#
#pd = vtkPolyData();
#pd.Allocate(nPoints, nPoints);
#
#pts = vtkPoints()
#pts = grid.GetPoints()
#
#pd.SetPoints(pts)
#pd.Update()
#
## help(grid.__class__)
#
#
#
#ptOn = True
#pt = [0,0,0]
#a = 0
#b = 0
#c = 0
#
#
#
#
#cells = vtkCellArray()
#cells.Initialize()
#cells.Allocate(cells.EstimateSize(3*nPoints, 2), 0)
#cells.InitTraversal()
#
#for n in range(nPoints):
#    pt = grid.GetPoint(n)
#    i = n % ni
#    j = (n / ni) %  nj
#    k = n / (ni * nj)
#    if i < ni-1:
#        cells.InsertNextCell(2)
#        cells.InsertCellPoint(n)
#        cells.InsertCellPoint(n+1)
#    if j < nj-1:
#        cells.InsertNextCell(2)
#        cells.InsertCellPoint(n)
#        cells.InsertCellPoint(n+ni)
#    if k < nk-1:
#        cells.InsertNextCell(2)
#        cells.InsertCellPoint(n)
#        cells.InsertCellPoint(n+ni*nj)
#
#
#
#pd.SetLines(cells)
#pd.Update()
#
#
#
#tubes = vtkTubeFilter()
##help(tubes.__class__)
#edges = vtkExtractEdges()
#edges.SetInput(grid)
#
#tubes.SetInputConnection(edges.GetOutputPort())
#tubes.Update()
#
#
#sphere = vtkSphereSource()
#sphere.SetRadius(3)
#sphere.SetPhiResolution(8)
#sphere.SetThetaResolution(8)
#
#
#glyphs = vtkGlyph3D()
#glyphs.SetInput(grid)
#glyphs.SetSourceConnection(sphere.GetOutputPort())
#
#
#
#appended = vtkAppendPolyData()
#appended.AddInputConnection(tubes.GetOutputPort())
#appended.AddInputConnection(glyphs.GetOutputPort())
#appended.Update()
#
#pdWriter = vtkPolyDataWriter()
#pdWriter.SetFileName('bla.vtk')
##pdWriter.SetInput(pd)
##pdWriter.SetInputConnection(tubes.GetOutputPort())
#pdWriter.SetInputConnection(appended.GetOutputPort())
#pdWriter.Update()
#
#
#
#
#
## nCells = grid.GetNumberOfCells()
## for n in range(nCells):
##     currCell = grid.GetCell(n)
##     # grid.GetCellType(n)
##     # print currCell
#
## celltypes = vtkCellTypes()
## help(celltypes.__class__)
#
## print grid.GetCellTypes(celltypes)
#
#
## n = 0
## for k in range(nk):
##     for j in range(nj):
##         for i in range(ni):
##             pt = grid.GetPoint(n)
##             print n, ' - ',  pt
##             n = n + 1
#
#
#
#
##help (grid.__class__)
