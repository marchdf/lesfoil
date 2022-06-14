# Pointwise V18.4R2 Journal file - Tue Sep 21 20:17:02 2021

package require PWI_Glyph 4.18.4

pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}

pw::Application clearModified

set dir /Users/mhenryde/exawind/cases/lesfoil/meshes
pw::Application setCAESolver {EXODUS II} 3
pw::Application markUndoLevel {Set Dimension 3D}

# Number of points on pressure or suction side

set ref coarse

if {$ref == "coarse"} {
   set np_pressure 101
   set np_suction 101
}

# Growth factor for initial wall normal extrusion
set gf 1.10

# Height for wall normal extrusion
set bl_dist 120.0

# First cell height
set ds 0.000011813977015662547

# Trailing and leading edge spacing
set te_dx 0.01
set le_dx 0.001

# Dimension of mesh
set dim 3

# Span length
set spanlength {0 0 0.05}

# Span spacing
set spandx 0.01
set spansteps [expr {int([lindex $spanlength end] / $spandx)}]

# length of branch connector
set branch_offset {50 0 0}

# Angle of attack
set angle -13.3

set _TMP(mode_1) [pw::Application begin DatabaseImport]
  $_TMP(mode_1) initialize -strict -type Segment $dir/lesfoil_lower.dat
  $_TMP(mode_1) read
  $_TMP(mode_1) convert
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Import Database}

set _TMP(mode_1) [pw::Application begin DatabaseImport]
  $_TMP(mode_1) initialize -strict -type Segment $dir/lesfoil_upper.dat
  $_TMP(mode_1) read
  $_TMP(mode_1) convert
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Import Database}

set _DB(1) [pw::DatabaseEntity getByName curve-2]
set _DB(2) [pw::DatabaseEntity getByName curve-1]
set _TMP(mode_1) [pw::Application begin Modify [list $_DB(1) $_DB(2)]]
  pw::Entity transform [pwu::Transform rotation -anchor {0 0 0} {0 0 1} $angle] [$_TMP(mode_1) getEntities]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Rotate

set _DB(1) [pw::DatabaseEntity getByName curve-2]
set _DB(2) [pw::DatabaseEntity getByName curve-1]
set _TMP(PW_1) [pw::Connector createOnDatabase -parametricConnectors Aligned -merge 0 -reject _TMP(unused) [list $_DB(1) $_DB(2)]]
unset _TMP(unused)
unset _TMP(PW_1)
pw::Application markUndoLevel {Connectors On DB Entities}

set _CN(1) [pw::GridEntity getByName con-2]
set _CN(2) [pw::GridEntity getByName con-1]
set _TMP(PW_1) [pw::Collection create]
$_TMP(PW_1) set [list $_CN(1)]
$_TMP(PW_1) do setDimension $np_suction
$_TMP(PW_1) delete
unset _TMP(PW_1)
pw::CutPlane refresh
pw::Application markUndoLevel Dimension

set _TMP(PW_1) [pw::Collection create]
$_TMP(PW_1) set [list $_CN(2)]
$_TMP(PW_1) do setDimension $np_pressure
$_TMP(PW_1) delete
unset _TMP(PW_1)
pw::CutPlane refresh
pw::Application markUndoLevel Dimension

pw::Display resetView -Z
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(1) $_CN(2)]]
  set _TMP(PW_1) [$_CN(2) getDistribution 1]
  $_TMP(PW_1) setBeginSpacing $te_dx
  unset _TMP(PW_1)
  set _TMP(PW_1) [$_CN(2) getDistribution 1]
  $_TMP(PW_1) setEndSpacing $le_dx
  unset _TMP(PW_1)
  set _TMP(PW_1) [$_CN(1) getDistribution 1]
  $_TMP(PW_1) setBeginSpacing $le_dx
  unset _TMP(PW_1)
  set _TMP(PW_1) [$_CN(1) getDistribution 1]
  $_TMP(PW_1) setEndSpacing $te_dx
  unset _TMP(PW_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Change Spacings}


# Branch connector

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint [$_CN(1) getPosition -arc 1]
  $_TMP(PW_1) addPoint [pwu::Vector3 add [pw::Application getXYZ [$_CN(1) getPosition -arc 1]] $branch_offset]
  set _CN(3) [pw::Connector create]
  $_CN(3) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(3) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
$_CN(3) setDimension 51
pw::CutPlane refresh
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(3)]]
  set _TMP(PW_1) [$_CN(3) getDistribution 1]
  $_TMP(PW_1) setBeginSpacing $te_dx
  unset _TMP(PW_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Change Spacing}

# Extrude normal to airfoil
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::Edge createFromConnectors [list $_CN(2) $_CN(1) $_CN(3)]]
  set _TMP(edge_1) [lindex $_TMP(PW_1) 0]
  unset _TMP(PW_1)
  $_TMP(edge_1) delete
  unset _TMP(edge_1)
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(3)
  $_TMP(edge_1) addConnector $_CN(1)
  $_TMP(edge_1) addConnector $_CN(2)
  $_TMP(edge_1) addConnector $_CN(3)
  set _DM(1) [pw::DomainStructured create]
  $_DM(1) addEdge $_TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_DM(1)]]
  $_TMP(mode_1) setKeepFailingStep true
  $_DM(1) setExtrusionSolverAttribute NormalMarchingVector {-0 -0 -1}
  $_DM(1) setExtrusionBoundaryCondition Begin ConstantX
  $_DM(1) setExtrusionBoundaryCondition End ConstantX
  $_DM(1) setExtrusionBoundaryConditionStepSuppression Begin 0
  $_DM(1) setExtrusionSolverAttribute NormalInitialStepSize $ds
  $_DM(1) setExtrusionSolverAttribute StopAtHeight Off
  $_DM(1) setExtrusionSolverAttribute NormalVolumeSmoothing 0.01
  $_DM(1) setExtrusionSolverAttribute SpacingGrowthFactor $gf
  $_DM(1) setExtrusionSolverAttribute StopAtHeight 1.0
  $_TMP(mode_1) run 500
  $_DM(1) setExtrusionSolverAttribute StopAtHeight $bl_dist
  $_DM(1) setExtrusionSolverAttribute NormalVolumeSmoothing 0.5
  $_DM(1) setExtrusionSolverAttribute SpacingGrowthFactor 1.2
  $_TMP(mode_1) run 100
$_TMP(mode_1) end
unset _TMP(mode_1)
#$_DM(1) delete
unset _TMP(edge_1)
pw::Display resetView -Z

# Wake block

set _TMP(mode_1) [pw::Application begin Create]
  set _CN(4) [pw::GridEntity getByName con-6]
  set _CN(5) [pw::GridEntity getByName con-4]
  set _TMP(PW_1) [pw::Edge createFromConnectors [list $_CN(4) $_CN(5)]]
  set _TMP(edge_1) [lindex $_TMP(PW_1) 0]
  unset _TMP(PW_1)
  set _DM(2) [pw::DomainStructured create]
  $_DM(2) addEdge $_TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_DM(2)]]
  $_TMP(mode_1) setKeepFailingStep true
  $_DM(2) setExtrusionSolverAttribute Mode Translate
  $_DM(2) setExtrusionSolverAttribute TranslateDirection {1 0 0}
  $_DM(2) setExtrusionSolverAttribute TranslateDistance 100
  $_TMP(mode_1) run 20
$_TMP(mode_1) end
unset _TMP(mode_1)
unset _TMP(edge_1)
pw::Application markUndoLevel {Extrude, Translate}

set _CN(6) [pw::GridEntity getByName con-7]
set _CN(7) [pw::GridEntity getByName con-9]
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(6) $_CN(7)]]
  pw::Connector swapDistribution Tanh [list [list $_CN(6) 1] [list $_CN(7) 1]]
  [[$_CN(7) getDistribution 1] getBeginSpacing] setValue 5
  [[$_CN(6) getDistribution 1] getBeginSpacing] setValue 5
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute


# Create the span connector

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint [$_CN(1) getPosition -arc 0]
  $_TMP(PW_1) addPoint $spanlength
  set _CN(8) [pw::Connector create]
  $_CN(8) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(8) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

$_CN(8) setDimensionFromSpacing -resetDistribution $spandx
pw::CutPlane refresh
pw::Application markUndoLevel Dimension


# Extrude in span

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::FaceStructured createFromDomains [list $_DM(2) $_DM(1)]]
  set _TMP(face_1) [lindex $_TMP(PW_1) 0]
  set _TMP(face_2) [lindex $_TMP(PW_1) 1]
  unset _TMP(PW_1)
  set _BL(1) [pw::BlockStructured create]
  $_BL(1) addFace $_TMP(face_1)
  set _BL(2) [pw::BlockStructured create]
  $_BL(2) addFace $_TMP(face_2)
$_TMP(mode_1) end
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_BL(1) $_BL(2)]]
  $_TMP(mode_1) setKeepFailingStep true
  $_BL(1) setExtrusionSolverAttribute Mode Path
  $_BL(2) setExtrusionSolverAttribute Mode Path
  $_BL(1) setExtrusionSolverAttribute PathConnectors [list $_CN(8)]
  $_BL(1) setExtrusionSolverAttribute PathUseTangent 1
  $_BL(2) setExtrusionSolverAttribute PathConnectors [list $_CN(8)]
  $_BL(2) setExtrusionSolverAttribute PathUseTangent 1
  $_TMP(mode_1) run $spansteps
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Extrude, Path}

unset _TMP(face_2)
unset _TMP(face_1)
pw::Display resetView -Z

# Volume condition

set _TMP(PW_1) [pw::VolumeCondition create]
pw::Application markUndoLevel {Create VC}

$_TMP(PW_1) setName Flow
pw::Application markUndoLevel {Name VC}

$_TMP(PW_1) apply [list $_BL(1) $_BL(2)]
pw::Application markUndoLevel {Set VC}

unset _TMP(PW_1)

# Boundary conditions

set _DM(3) [pw::GridEntity getByName dom-4]
set _DM(4) [pw::GridEntity getByName dom-5]
set _DM(5) [pw::GridEntity getByName dom-8]
set _DM(6) [pw::GridEntity getByName dom-10]
set _DM(7) [pw::GridEntity getByName dom-13]
set _DM(8) [pw::GridEntity getByName dom-14]
set _DM(9) [pw::GridEntity getByName dom-15]
set _DM(10) [pw::GridEntity getByName dom-16]
set _DM(11) [pw::GridEntity getByName dom-3]
set _DM(12) [pw::GridEntity getByName dom-7]
set _DM(13) [pw::GridEntity getByName dom-9]
set _TMP(PW_1) [pw::BoundaryCondition getByName Unspecified]
set _TMP(PW_2) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_3) [pw::BoundaryCondition getByName bc-2]
unset _TMP(PW_2)
$_TMP(PW_3) setName inflow
pw::Application markUndoLevel {Name BC}

set _TMP(PW_4) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_5) [pw::BoundaryCondition getByName bc-3]
unset _TMP(PW_4)
$_TMP(PW_5) setName outflow
pw::Application markUndoLevel {Name BC}

set _TMP(PW_6) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_7) [pw::BoundaryCondition getByName bc-4]
unset _TMP(PW_6)
$_TMP(PW_7) setName airfoil
pw::Application markUndoLevel {Name BC}

set _TMP(PW_8) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_9) [pw::BoundaryCondition getByName bc-5]
unset _TMP(PW_8)
$_TMP(PW_9) setName front
pw::Application markUndoLevel {Name BC}

set _TMP(PW_10) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

set _TMP(PW_11) [pw::BoundaryCondition getByName bc-6]
unset _TMP(PW_10)
$_TMP(PW_11) setName back
pw::Application markUndoLevel {Name BC}

$_TMP(PW_3) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_5) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_7) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_9) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_11) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_3) apply [list [list $_BL(1) $_DM(5)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_5) apply [list [list $_BL(2) $_DM(8)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_9) apply [list [list $_BL(1) $_DM(6)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_9) apply [list [list $_BL(2) $_DM(10)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_11) apply [list [list $_BL(1) $_DM(1)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_11) apply [list [list $_BL(2) $_DM(2)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_7) apply [list [list $_BL(1) $_DM(4)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_7) apply [list [list $_BL(1) $_DM(3)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_5) apply [list [list $_BL(2) $_DM(9)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_5) apply [list [list $_BL(2) $_DM(7)]]
pw::Application markUndoLevel {Set BC}

unset _TMP(PW_1)
unset _TMP(PW_3)
unset _TMP(PW_5)
unset _TMP(PW_7)
unset _TMP(PW_9)
unset _TMP(PW_11)

# Save the mesh

set _TMP(mode_1) [pw::Application begin CaeExport [pw::Entity sort [list $_BL(1) $_BL(2)]]]
  $_TMP(mode_1) initialize -strict -type CAE $dir/lesfoil-c-$ref.exo
  $_TMP(mode_1) verify
  $_TMP(mode_1) write
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application save $dir/lesfoil-c-$ref.pw
