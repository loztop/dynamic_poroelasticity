try: paraview.simple
except: from paraview.simple import *
paraview.simple._DisableFirstRenderCameraReset()

SIM_SquashSponge01_es_002 = ExodusIIReader( FileName=['/users/lorenzb/libmesh/chapelle_press/SIM_SquashSponge01.e-s.002'] )

AnimationScene1 = GetAnimationScene()
SIM_SquashSponge01_es_002.NodeMapArrayStatus = []
SIM_SquashSponge01_es_002.FaceVariables = []
SIM_SquashSponge01_es_002.ElementVariables = []
SIM_SquashSponge01_es_002.XMLFileName = 'Invalid result'
SIM_SquashSponge01_es_002.FaceSetResultArrayStatus = []
SIM_SquashSponge01_es_002.PointVariables = []
SIM_SquashSponge01_es_002.FaceSetArrayStatus = []
SIM_SquashSponge01_es_002.FaceMapArrayStatus = []
SIM_SquashSponge01_es_002.FileRange = [0, 0]
SIM_SquashSponge01_es_002.SideSetResultArrayStatus = []
SIM_SquashSponge01_es_002.ElementSetArrayStatus = []
SIM_SquashSponge01_es_002.EdgeVariables = []
SIM_SquashSponge01_es_002.FilePrefix = '/users/lorenzb/libmesh/chapelle_press/SIM_SquashSponge01.e-s.002'
SIM_SquashSponge01_es_002.FilePattern = '%s'
SIM_SquashSponge01_es_002.EdgeSetArrayStatus = []
SIM_SquashSponge01_es_002.SideSetArrayStatus = []
SIM_SquashSponge01_es_002.GlobalVariables = []
SIM_SquashSponge01_es_002.NodeSetArrayStatus = []
SIM_SquashSponge01_es_002.NodeSetResultArrayStatus = []
SIM_SquashSponge01_es_002.ElementMapArrayStatus = []
SIM_SquashSponge01_es_002.EdgeSetResultArrayStatus = []
SIM_SquashSponge01_es_002.EdgeMapArrayStatus = []
SIM_SquashSponge01_es_002.ElementSetResultArrayStatus = []

AnimationScene1.EndTime = 1.4142857142857135
AnimationScene1.PlayMode = 'Snap To TimeSteps'

SIM_SquashSponge01_es_002.PointVariables = ['u', 'v', 'w', 'p', 'a', 'b', 'c', 'u_nu', 'v_nu', 'w_nu', 'p_nu', 'a_nu', 'b_nu', 'c_nu', 'u_ref', 'v_ref', 'w_ref', 'p_ref', 'a_ref', 'b_ref', 'c_ref', 'fluid_U_vel', 'fluid_V_vel', 'fluid_W_vel', 'fluid_P', 'fluid_m', 'Jacobian', 'res_aux1', 'res_aux2']
SIM_SquashSponge01_es_002.ElementBlocks = ['Unnamed block ID: 0 Type: HEX27 Size: 63']

RenderView1 = GetRenderView()
DataRepresentation12 = Show()
DataRepresentation12.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation12.SelectionPointFieldDataArrayName = 'a'
DataRepresentation12.SelectionCellFieldDataArrayName = 'GlobalElementId'
DataRepresentation12.ScalarOpacityUnitDistance = 1.136497396100493
DataRepresentation12.ExtractedBlockIndex = 2
DataRepresentation12.ScaleFactor = 0.39932123213075105

defe5 = defe()

RenderView1.CameraClippingRange = [9.149961855585381, 14.455555965471646]

DataRepresentation13 = Show()
DataRepresentation13.EdgeColor = [0.0, 0.0, 0.5000076295109483]
DataRepresentation13.SelectionPointFieldDataArrayName = 'a'
DataRepresentation13.SelectionCellFieldDataArrayName = 'GlobalElementId'
DataRepresentation13.ScalarOpacityUnitDistance = 0.0
DataRepresentation13.ExtractedBlockIndex = 2
DataRepresentation13.ScaleFactor = 0.1

RenderView1.CameraClippingRange = [12.47675066499613, 12.791820126233407]

DataRepresentation12.Visibility = 0

Render()
