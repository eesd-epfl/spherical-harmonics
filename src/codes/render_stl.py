import vtk, random

def render(filename=""):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(filename)
    reader.SetFileName(filename)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(reader.GetOutputPort())
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    # Create a rendering window and renderer
    ren = vtk.vtkRenderer()
    renWin = vtk.vtkRenderWindow()
    renWin.AddRenderer(ren)
    # Create a renderwindowinteractor
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)
    # Assign actor to the renderer
    ren.AddActor(actor)

    # The object color
    # Set properties from here: https://vtk.org/doc/nightly/html/classvtkProperty.html
    actor.GetProperty().SetColor(0.9, 0.9, 0.9)
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetEdgeColor(0.1,0.1,0.1)
    ren.SetBackground(1., 1., 1.)
    
    # Enable user interface interactor
    iren.Initialize()
    renWin.Render()
    iren.Start()
    # delete the screen
    # End the render screen from VTK
    iren.GetRenderWindow().Finalize()
    iren.TerminateApp()
    del renWin, iren, ren, actor, mapper, reader
    
if __name__ == "__main__":
    filename = "test_stone.stl"
    render(filename)
    
