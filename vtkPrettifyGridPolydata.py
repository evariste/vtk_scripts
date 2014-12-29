# run with vtkpython after setting python path etc using source script.

import sys
import vtk

def main(*args):
  inputName  = args[1]
  outputName = args[2]

  reader = vtk.vtkPolyDataReader()
  reader.SetFileName(inputName)
  reader.Update()

  pd = vtk.vtkPolyData()
  pd = reader.GetOutput()
  pd.Update()

  tubes = vtk.vtkTubeFilter()
  tubes.SetInput(pd)
  tubes.SetRadius(1)
  tubes.Update()
  tubeOutput = tubes.GetOutput()

  tubeLabels = vtk.vtkFloatArray()
  tubeLabels.SetNumberOfComponents(1)
  tubeLabels.SetNumberOfTuples(tubeOutput.GetNumberOfPoints())
  tubeLabels.SetName('label')

  for n in range(tubeOutput.GetNumberOfPoints()):
    tubeLabels.SetTuple1(n, 0)

  tubeOutput.GetPointData().AddArray(tubeLabels)
  tubeOutput.Update()

  sphere = vtk.vtkSphereSource()
  sphere.SetRadius(3)
  sphere.SetPhiResolution(8)
  sphere.SetThetaResolution(8)

  glyphs = vtk.vtkGlyph3D()
  glyphs.SetInput(pd)
  glyphs.SetSourceConnection(sphere.GetOutputPort())
  glyphs.Update()


  glyphOutput = glyphs.GetOutput()

  glyphLabels = vtk.vtkFloatArray()
  glyphLabels.SetNumberOfComponents(1)
  glyphLabels.SetNumberOfTuples(glyphOutput.GetNumberOfPoints())
  glyphLabels.SetName('label')
  
  
  for n in range(glyphOutput.GetNumberOfPoints()):
    glyphLabels.SetTuple1(n, 1)

  glyphOutput.GetPointData().AddArray(glyphLabels)
  glyphOutput.Update()

  appended = vtk.vtkAppendPolyData()
  appended.AddInputConnection(tubes.GetOutputPort())
  appended.AddInput(glyphOutput)
  appended.Update()

  cleaner = vtk.vtkCleanPolyData()
  cleaner.SetInputConnection(appended.GetOutputPort())
  
  pdWriter = vtk.vtkPolyDataWriter()
  pdWriter.SetFileName(outputName)
  pdWriter.SetInputConnection(cleaner.GetOutputPort())
  pdWriter.Update()



if __name__ == '__main__':
  if len(sys.argv) < 3:
    sys.exit(1)
  sys.exit(main(*sys.argv))

