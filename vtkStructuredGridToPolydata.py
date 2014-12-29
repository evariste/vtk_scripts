# run with vtkpython after setting python path etc using source script.

import sys
import vtk

def main(*args):

  # Free form deformation that has been converted to vtkStructuredGrid
  # format using dof2vtk:
  inputName  = args[1]

  # Polydata file:
  outputName = args[2]


  reader = vtk.vtkStructuredGridReader()
  reader.SetFileName(inputName)
  reader.Update()

  grid = vtk.vtkStructuredGrid()
  grid = reader.GetOutput()
  grid.Update()
  

  nPoints = grid.GetNumberOfPoints()

  pd = vtk.vtkPolyData();
  pd.Allocate(nPoints, nPoints);

  # Assign points to polydata
  pts = vtk.vtkPoints()
  pts = grid.GetPoints()
  
  
  pd.SetPoints(pts)
  pd.Update()
  pd.Squeeze()

  # Prepare for assigning edges
  dims    = grid.GetDimensions()
  ni = dims[0]
  nj = dims[1]
  nk = dims[2]

  cells = vtk.vtkCellArray()
  cells.Initialize()
  cells.Allocate(cells.EstimateSize(3*nPoints, 2), 0)
  cells.InitTraversal()

  # Insert an edge between each pair of grid-adjacent points.
  for n in range(nPoints):
    pt = grid.GetPoint(n)
    i = n % ni
    j = (n / ni) %  nj
    k = n / (ni * nj)
    if i < ni-1:
      cells.InsertNextCell(2)
      cells.InsertCellPoint(n)
      cells.InsertCellPoint(n+1)
    if j < nj-1:
      cells.InsertNextCell(2)
      cells.InsertCellPoint(n)
      cells.InsertCellPoint(n+ni)
    if k < nk-1:
      cells.InsertNextCell(2)
      cells.InsertCellPoint(n)
      cells.InsertCellPoint(n+ni*nj)

  ptData = vtk.vtkFloatArray()
  ptData = grid.GetPointData().GetVectors()

  pd.GetPointData().AddArray(ptData)
  pd.SetLines(cells)
  pd.Update()

  # Finish off.
  pdWriter = vtk.vtkPolyDataWriter()
  pdWriter.SetFileName(outputName)
  pdWriter.SetInput(pd)
  pdWriter.Update()



if __name__ == '__main__':
  if len(sys.argv) < 3:
    print "Usage $0 input output"
    sys.exit(1)
  sys.exit(main(*sys.argv))

