#!python

import vtk

# Create 6 points.
p0 = [0.0, 0.0, 0.0]
p1 = [1.0, 0.0, 0.0]
p2 = [1.0, 1.0, 0.0]
p3 = [1.0, 1.0, 1.0]
p4 = [2.0, 2.0, 2.0]
p5 = [1.0, 0.0, 0.0]
# Create a vtkPoints object and store the points in it
pts = vtk.vtkPoints()
pts.InsertNextPoint(p0)
pts.InsertNextPoint(p1)
pts.InsertNextPoint(p2)
pts.InsertNextPoint(p3)
pts.InsertNextPoint(p4)
pts.InsertNextPoint(p5)
# Setup four colors - one for each line
red = [255, 0, 0]
green = [0, 255, 0]
white = [255, 255, 255]
sandybrown = [244,164,96]
blue = [0, 0, 255]
# Setup the colors array
colors = vtk.vtkUnsignedCharArray()
#R,G,B colors
colors.SetNumberOfComponents(3)
colors.SetName("MyColors")

#Setup radius for each point
widths = vtk.vtkDoubleArray()
widths.SetName('TubeRadius')
#widths.SetNumberOfTuples(5)
widths.InsertNextTuple1(0.10)
widths.InsertNextTuple1(0.05)
widths.InsertNextTuple1(0.03)
widths.InsertNextTuple1(0.01)
widths.InsertNextTuple1(0.10)
widths.InsertNextTuple1(0.01)
#widths.SetTuple1(1,0.05)
#widths.SetTuple1(2,0.03)
#widths.SetTuple1(3,0.01)
#widths.SetTuple1(4,0.08)
# Add the colors we created to the colors array
colors.InsertNextTypedTuple(red)
colors.InsertNextTypedTuple(green)
colors.InsertNextTypedTuple(sandybrown)
colors.InsertNextTypedTuple(blue)
# Create the first line (between Origin and P0)
line0 = vtk.vtkLine()
line0.GetPointIds().SetId(0,0) # the second 0 is the index of the P0 in the vtkPoints
line0.GetPointIds().SetId(1,5) # the second 1 is the index of P5 in the vtkPoints
# Create the second line (between Origin and P1)
line1 = vtk.vtkLine()
line1.GetPointIds().SetId(0,1) # the 1 is the index of the P1 in the vtkPoints
line1.GetPointIds().SetId(1,2) # the 2 is the index of P2 in the vtkPoints
#Third line
line2 = vtk.vtkLine()
line2.GetPointIds().SetId(0,2) #The index of P2 
line2.GetPointIds().SetId(1,3) #The index of P3
#Fourth Line
line3 = vtk.vtkLine()
line3.GetPointIds().SetId(0,0) #The index of P0
line3.GetPointIds().SetId(0,4) #The index of P4 
# Create a cell array to store the lines in and add the lines to it
lines = vtk.vtkCellArray()
lines.InsertNextCell(line0)
lines.InsertNextCell(line1)
lines.InsertNextCell(line2)
lines.InsertNextCell(line3)
# Create a polydata to store everything in it
linesPolyData = vtk.vtkPolyData()
# Add the points to the dataset
linesPolyData.SetPoints(pts)
# Add the lines to the dataset
linesPolyData.SetLines(lines)
# Color the lines - associate the first component (red) of the
# colors array with the first component of the cell array (line 0)
# and the second component (green) of the colors array with the
# second component of the cell array (line 1)
linesPolyData.GetCellData().AddArray(colors)
#linesPolyData.GetCellData().SetScalars(colors)
#linesPolyData.GetCellData().SetActiveScalars("MyColors");
#Add widths
linesPolyData.GetPointData().AddArray(widths)
linesPolyData.GetPointData().SetActiveScalars("TubeRadius");

#Display tubes instead of lines
tubeFilter = vtk.vtkTubeFilter()
tubeFilter.SetInputData(linesPolyData)
#tubeFilter.SetRadius(0.03)
tubeFilter.SetNumberOfSides(10)
tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()
tubeFilter.CappingOn()
tubeFilter.Update()
writer = vtk.vtkXMLPolyDataWriter()
writer.SetFileName('ColorLines.vtp')
#writer.SetInputDataObject(linesPolyData)
#writer.SetInputDataObject(tubeFilter.GetOutput())
writer.AddInputDataObject(tubeFilter.GetOutput())
writer.Write()
# Visualize
mapper = vtk.vtkPolyDataMapper()
if vtk.VTK_MAJOR_VERSION <= 5:
    mapper.SetInput(tubeFilter.GetOutput())
else:
    #mapper.SetInputData(tubeFilter.GetOutput())
    mapper.SetInputConnection(tubeFilter.GetOutputPort())
mapper.ScalarVisibilityOn()
mapper.SetScalarModeToUseCellFieldData()
mapper.SelectColorArray("MyColors")
actor = vtk.vtkActor()
actor.SetMapper(mapper)

renderer = vtk.vtkRenderer()
renderer.AddActor(actor)
renderWindow = vtk.vtkRenderWindow()
renderWindow.AddRenderer(renderer)
renderWindowInteractor = vtk.vtkRenderWindowInteractor()
renderWindowInteractor.SetRenderWindow(renderWindow)


renderWindow.Render()
renderWindowInteractor.Start()
