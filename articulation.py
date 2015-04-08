import vtk

# Reads in 19 landmarks for a fetal skeleton and adds links for where they should be joined.

# The landmarks from the top (head) as follows:
#
#               1
#               |
# 2 -- 3 -- 4 - 5 - 6 -- 7 -- 8 
#               |
#               |
#               9
#               |
#               |
#              10
#               |
#               |
#          12 - 11 - 16
#          /          \
#         /            \
#       13              17
#       |               |
#       |               |
#      14               18
#      /                  \
#    15                    19
reader = vtk.vtkPolyDataReader()

reader.SetFileName('landmarks-reset.vtk')
reader.Update()


bodyPts = vtk.vtkPolyData()
bodyPts = reader.GetOutput()
bodyPts.Update()

cells = vtk.vtkCellArray()
cells.Initialize()
cells.Allocate(100, 0)
cells.InitTraversal()

# One-indexed:
edges = ([[1, 5], # head 
          [2, 3], [3, 4], [4, 5], # arm
          [5, 6], [6, 7], [7, 8], # arm
          [5, 9], [9, 10], [10,11], # spine
          [11,12], [12,13], [13,14], [14,15], # leg
          [11,16], [16,17], [17,18], [18,19] ]) # leg

for i in range(len(edges)):
  cells.InsertNextCell(2)
  cells.InsertCellPoint(edges[i][0]-1)
  cells.InsertCellPoint(edges[i][1]-1)

cells.Squeeze()

bodyPts.SetLines(cells)

writer = vtk.vtkPolyDataWriter()
writer.SetFileName('articulation.vtk')
writer.SetInput(bodyPts)
writer.Write()



print bodyPts.GetNumberOfPoints()
