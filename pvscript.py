from paraview.simple import *
sm= servermanager
sm.Connect()
sphere = Sphere()
sphere.ThetaResolution= 100
sphere.PhiResolution= 100
view = sm.CreateRenderView()
view.ResetCamera()


shrinkFilter = Shrink(sphere)
Show(shrinkFilter)
GetActiveCamera().Elevation(105)
Render()

Render()

#pidscal= ProcessIdScalars(sphere)
#view = sm.CreateRenderView()
#rep = Show(pidscal, view)
#nbprocs= sm.ActiveConnection.GetNumberOfDataPartitions
#lt= MakeBlueToRedLT(0, nbprocs-1)
#lt.NumberOfTableValues= nbprocs
#rep.LookupTable= lt
#rep.ColorAttributeType= 'POINT_DATA'
#rep.ColorArrayName= "ProcessId"
#bar = CreateScalarBar(LookupTable=lt, Title='PID')
#view.Representations.append(bar)
#view.ResetCamera()
#view.StillRender()
#view.WriteImage("Sphere.png","vtkPNGWriter",1)
