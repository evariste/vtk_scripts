# run with vtkpython after setting python path etc using source script.

import sys
import vtk

def main(*args):

  # Free form deformation that has been converted to vtkStructuredGrid
  # format using dof2vtk:
  inputName  = args[1]

  # Polydata file:
  outputName = args[2]


  reader = vtk.vtkUnstructuredGridReader()
  reader.SetFileName(inputName)
  reader.Update()

  dataIn = vtk.vtkUnstructuredGrid()
  dataIn = reader.GetOutput()
  dataIn.Update()
  
  geomFilter = vtk.vtkGeometryFilter()
  geomFilter.SetInput(dataIn)
  geomFilter.Update()
  
  pd = vtk.vtkPolyData();
  pd = geomFilter.GetOutput()
  pd.Update()
  
  
  # Finish off.
  pdWriter = vtk.vtkPolyDataWriter()
  pdWriter.SetFileName(outputName)
  pdWriter.SetInput(pd)
  pdWriter.Update()



if __name__ == '__main__':
  if len(sys.argv) < 3:
    print 'Usage :   vtkpython', sys.argv[0] ,' input output'
    sys.exit(1)
  sys.exit(main(*sys.argv))

