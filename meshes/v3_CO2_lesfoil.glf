# Pointwise V18.4R2 Journal file - Tue Sep 21 20:17:02 2021
# (Nov 2022) Modification by S. Haering of original script by
# Marc de Frahan to convert full C-mesh into C-mesh embedded in an O-mesh


#package require PWI_Glyph 4.18.4
package require PWI_Glyph 6.22.1

pw::Application setUndoMaximumLevels 5
pw::Application reset
pw::Application markUndoLevel {Journal Reset}
pw::Application clearModified

set dir /Users/mhenryde/exawind/cases/lesfoil/meshes
pw::Application setCAESolver {EXODUS II} 3
pw::Application markUndoLevel {Set Dimension 3D}


# Number of points on pressure or suction side
set ref fine

if {$ref == "coarse"} {
   set np_pressure 101
   set np_suction 101
   set spandx 0.01

   # Trailing and leading edge spacing
   #set te_dx 0.01
   set te_dx 0.001
   set le_dx 0.005
   set branch_number 3
}
if {$ref == "fine"} {
   #set np_pressure 101
   ###set np_pressure 201
   set np_pressure 401
   #set np_suction 801
   ###set np_suction 901
   set np_suction 1501
   set spandx 0.0025

   # Trailing and leading edge spacing
   #set te_dx 0.0025
   ###set te_dx 0.00025
   set te_dx 0.0002
   #set le_dx 0.00025
   ###set le_dx 0.0005
   set le_dx 0.0001
   set branch_number 5
}
if {$ref == "finer"} {
   #set np_pressure 101
   set np_pressure 301
   set np_suction 1601
   set spandx 0.0025

   # Trailing and leading edge spacing
   #set te_dx 0.0025
   set te_dx 0.00025
   #set le_dx 0.00025
   set le_dx 0.0003
   set branch_number 5
}

# Growth factor for initial wall normal extrusion
set gf 1.08
set gf2 1.05

# Height for wall normal extrusion
set bl_dist 120.0

# First cell height (~1.2E-5)
set ds 0.000011813977015662547
#set ds 2e-5

# Dimension of mesh
set dim 3

# Span length
#set spanlength {0 0 0.05}
set spanlength 0.05
#set spanlength 0.1
#set spanlength 10.0

# Span spacing
#set spansteps [expr {int([lindex $spanlength end] / $spandx)}]
set spansteps [expr {int($spanlength / $spandx)}]
#set spansteps 4

# length of branch connector
#set branch_offset {50 0 0}
set branch_offset {0.001 0 0}

# discretization of "branch" => moved to mesh specs
#set branch_number 51

# Angle of attack
set angle -13.3

# outer params
set N_wake 100
set d_wake1 0.00025
set d_wake2 0.0004
set d_in 0.0008

# separtion pt spacing
set d_sep 0.0002

# and begin...

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
$_CN(3) setDimension $branch_number
pw::CutPlane refresh
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(3)]]
  set _TMP(PW_1) [$_CN(3) getDistribution 1]
  $_TMP(PW_1) setBeginSpacing $te_dx
  unset _TMP(PW_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Change Spacing}


# ADDITION FOR SEAPARATION POINT SPACING
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(1)]]
  $_CN(1) addBreakPoint -X 0.13
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(1)]]
  [[$_CN(1) getDistribution 2] getBeginSpacing] setValue $d_sep
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(1)]]
  [[$_CN(1) getDistribution 1] getEndSpacing] setValue $d_sep
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute



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
  $_DM(1) setExtrusionSolverAttribute StopAtHeight 0.005
  $_TMP(mode_1) run 500
  $_DM(1) setExtrusionSolverAttribute StopAtHeight 0.02
  $_DM(1) setExtrusionSolverAttribute NormalVolumeSmoothing 0.01
  $_DM(1) setExtrusionSolverAttribute SpacingGrowthFactor 1.0
  $_TMP(mode_1) run 500
#  $_DM(1) setExtrusionSolverAttribute StopAtHeight 0.1
#  $_DM(1) setExtrusionSolverAttribute StopAtHeight 0.05
#  $_DM(1) setExtrusionSolverAttribute NormalVolumeSmoothing 0.1
#  $_DM(1) setExtrusionSolverAttribute SpacingGrowthFactor $gf2
#  $_TMP(mode_1) run 500
#  $_DM(1) setExtrusionSolverAttribute StopAtHeight $bl_dist
#  $_DM(1) setExtrusionSolverAttribute NormalVolumeSmoothing 0.5
#  $_DM(1) setExtrusionSolverAttribute SpacingGrowthFactor 1.2
#  $_TMP(mode_1) run 500
$_TMP(mode_1) end
unset _TMP(mode_1)
#$_DM(1) delete
unset _TMP(edge_1)
pw::Display resetView -Z


# second extrusion to get break points
set _TMP(mode_1) [pw::Application begin Create]
  set _CN(4) [pw::GridEntity getByName con-5]
  set _TMP(PW_1) [pw::Edge createFromConnectors [list $_CN(4)]]
  set _TMP(edge_1) [lindex $_TMP(PW_1) 0]
  unset _TMP(PW_1)
  set _DM(2) [pw::DomainStructured create]
  $_DM(2) addEdge $_TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_DM(2)]]
  $_TMP(mode_1) setKeepFailingStep true
  $_DM(2) setExtrusionSolverAttribute SpacingGrowthFactor 1.08
  $_DM(2) setExtrusionSolverAttribute NormalInitialStepSize 0.0004
  $_DM(2) setExtrusionSolverAttribute StopAtHeight Off
  $_DM(2) setExtrusionSolverAttribute StopAtHeight 0.05
  $_DM(2) setExtrusionBoundaryCondition End Splay 0.005
  $_DM(2) setExtrusionBoundaryConditionStepSuppression End 0
  $_DM(2) setExtrusionBoundaryCondition Begin Splay 0.0075
  $_DM(2) setExtrusionBoundaryConditionStepSuppression Begin 0
  $_TMP(mode_1) run 32
  $_TMP(mode_1) run -1
  pw::Display setShowDomains 1
$_TMP(mode_1) end
unset _TMP(mode_1)
unset _TMP(edge_1)
pw::Application markUndoLevel {Extrude, Normal}


# build trailing edge transition box
pw::Display setShowDomains 0
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentCircle create]
  set _CN(5) [pw::GridEntity getByName con-8]
  set _CN(6) [pw::GridEntity getByName con-9]
  set _CN(7) [pw::GridEntity getByName con-7]
  $_TMP(PW_1) addPoint [$_CN(5) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(7) getPosition -arc 1]
  $_TMP(PW_1) setCenterPoint {0.959400110697872 -0.210849159470519 0} {0 0 1}
  set _CN(8) [pw::Connector create]
  $_CN(8) addSegment $_TMP(PW_1)
  $_CN(8) calculateDimension
  unset _TMP(PW_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentCircle create]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  set _CN(9) [pw::GridEntity getByName con-6]
  $_TMP(PW_1) addPoint [$_CN(4) getPosition -arc 1]
  $_TMP(PW_1) addPoint {1.005 -0.20853257 0}
  set _CN(10) [pw::Connector create]
  $_CN(10) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(10) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  set _CN(11) [pw::GridEntity getByName con-4]
  $_TMP(PW_1) addPoint [$_CN(11) getPosition -arc 1]
  $_TMP(PW_1) addPoint {1.005 -0.25955164 0}
  set _CN(12) [pw::Connector create]
  $_CN(12) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(12) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint [$_CN(12) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(10) getPosition -arc 1]
  set _CN(13) [pw::Connector create]
  $_CN(13) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(13) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(split_params) [list]
lappend _TMP(split_params) [$_CN(8) getParameter -closest [pw::Application getXYZ [list 1.04335055313673 -0.30458951398644574 0.0]]]
lappend _TMP(split_params) [$_CN(8) getParameter -closest [pw::Application getXYZ [list 1.0795818989902435 -0.17354944615661771 0.0]]]
set _TMP(PW_1) [$_CN(8) split $_TMP(split_params)]
unset _TMP(PW_1)
unset _TMP(split_params)
pw::Application markUndoLevel Split

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  set _CN(14) [pw::GridEntity getByName con-10-split-1]
  set _CN(15) [pw::GridEntity getByName con-10-split-2]
  $_TMP(PW_1) addPoint [$_CN(14) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(12) getPosition -arc 1]
  set _CN(16) [pw::Connector create]
  $_CN(16) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(16) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}


set _CN(17) [pw::GridEntity getByName con-10-split-3]
set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(15) $_CN(17)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _CN(18) [pw::GridEntity getByName con-10-split-2]
set _TMP(split_params) [list]
lappend _TMP(split_params) [$_CN(18) getParameter -closest [pw::Application getXYZ [list 1.0618453805444221 -0.13779856170429922 0.0]]]
set _TMP(PW_1) [$_CN(18) split $_TMP(split_params)]
unset _TMP(PW_1)
unset _TMP(split_params)
pw::Application markUndoLevel Split

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  set _CN(19) [pw::GridEntity getByName con-10-split-2-split-1]
  set _CN(20) [pw::GridEntity getByName con-10-split-2-split-2]
  $_TMP(PW_1) addPoint [$_CN(10) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(19) getPosition -arc 1]
  set _CN(21) [pw::Connector create]
  $_CN(21) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(21) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(13) $_CN(19)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(11) 1] [list $_CN(9) 1]]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(16) $_CN(21)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(7) 1]]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

# transition box grid size
set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(10) $_CN(8) $_CN(12) $_CN(20)]
  #$_TMP(PW_1) do setDimension -resetDistribution 60
  $_TMP(PW_1) do setDimension -resetDistribution 80
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Dimension]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(21)]]
  set _TMP(dist_1) [pw::DistributionGeneral create [list [list $_CN(7) 1]]]
  # Clear spacings so the distribution will scale properly
  $_TMP(dist_1) setBeginSpacing 0
  $_TMP(dist_1) setEndSpacing 0
  $_TMP(dist_1) setVariable [[$_CN(21) getDistribution 1] getVariable]
  $_CN(21) setDistribution -lockEnds 1 $_TMP(dist_1)
  unset _TMP(dist_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(21)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(16)]]
  set _TMP(dist_1) [pw::DistributionGeneral create [list [list $_CN(7) 1]]]
  $_TMP(dist_1) reverse
  # Clear spacings so the distribution will scale properly
  $_TMP(dist_1) setBeginSpacing 0
  $_TMP(dist_1) setEndSpacing 0
  $_TMP(dist_1) setVariable [[$_CN(16) getDistribution 1] getVariable]
  $_CN(16) setDistribution -lockEnds 1 $_TMP(dist_1)
  unset _TMP(dist_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(16)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(19)]]
  $_CN(19) removeAllBreakPoints
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(19)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(12)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(13)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(10)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(9)
  $_TMP(edge_4) addConnector $_CN(11)
  set _DM(3) [pw::DomainStructured create]
  $_DM(3) addEdge $_TMP(edge_1)
  $_DM(3) addEdge $_TMP(edge_2)
  $_DM(3) addEdge $_TMP(edge_3)
  $_DM(3) addEdge $_TMP(edge_4)
  unset _TMP(edge_4)
  unset _TMP(edge_3)
  unset _TMP(edge_2)
  unset _TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(8)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(16)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(12)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(6)
  set _DM(4) [pw::DomainStructured create]
  $_DM(4) addEdge $_TMP(edge_1)
  $_DM(4) addEdge $_TMP(edge_2)
  $_DM(4) addEdge $_TMP(edge_3)
  $_DM(4) addEdge $_TMP(edge_4)
  unset _TMP(edge_4)
  unset _TMP(edge_3)
  unset _TMP(edge_2)
  unset _TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(16)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(19)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(21)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(13)
  set _DM(5) [pw::DomainStructured create]
  $_DM(5) addEdge $_TMP(edge_1)
  $_DM(5) addEdge $_TMP(edge_2)
  $_DM(5) addEdge $_TMP(edge_3)
  $_DM(5) addEdge $_TMP(edge_4)
  unset _TMP(edge_4)
  unset _TMP(edge_3)
  unset _TMP(edge_2)
  unset _TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(21)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(20)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(7)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(10)
  set _DM(6) [pw::DomainStructured create]
  $_DM(6) addEdge $_TMP(edge_1)
  $_DM(6) addEdge $_TMP(edge_2)
  $_DM(6) addEdge $_TMP(edge_3)
  $_DM(6) addEdge $_TMP(edge_4)
  unset _TMP(edge_4)
  unset _TMP(edge_3)
  unset _TMP(edge_2)
  unset _TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_1) [pw::Application begin Create]
$_TMP(mode_1) abort
unset _TMP(mode_1)


# modify some spacing in transition box
pw::Display setShowDomains 0
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(10)]]
  [[$_CN(10) getDistribution 1] getBeginSpacing] setValue 0.0002
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(10)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(12)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(12)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(12)]]
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(10)]]
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(12)]]
  [[$_CN(12) getDistribution 1] getEndSpacing] setValue 0.002
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(12)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

pw::Display setShowDomains 1
pw::Display setShowDomains 0
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(21)]]
  [[$_CN(21) getDistribution 1] getBeginSpacing] setValue 0.0002
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(21)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(16)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(16)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

pw::Display setShowDomains 1
pw::Display setShowDomains 0
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(21)]]
  [[$_CN(21) getDistribution 1] getEndSpacing] setValue 0.003
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(21)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(16)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(16)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

pw::Display setShowDomains 1
pw::Display setShowDomains 0
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(21)]]
  $_CN(21) replaceDistribution 1 [pw::DistributionTanh create]
  [$_CN(21) getDistribution 1] setBeginSpacing 0.0
  [$_CN(21) getDistribution 1] setEndSpacing 0.0
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(21)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(16)]]
  $_CN(16) replaceDistribution 1 [pw::DistributionTanh create]
  [$_CN(16) getDistribution 1] setBeginSpacing 0.0
  [$_CN(16) getDistribution 1] setEndSpacing 0.0
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(16)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

pw::Display setShowDomains 1



# make outer O
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentCircle create]
  $_TMP(PW_1) addPoint {0 50 0}
  $_TMP(PW_1) addPoint {0 -50 0}
  $_TMP(PW_1) setCenterPoint {0 0 0} {0 0 1}
  set _CN(22) [pw::Connector create]
  $_CN(22) addSegment $_TMP(PW_1)
  $_CN(22) calculateDimension
  unset _TMP(PW_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentCircle create]
  $_TMP(PW_1) addPoint [$_CN(22) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(22) getPosition -arc 0]
  $_TMP(PW_1) setCenterPoint {0 0 0} {0 0 1}
  set _CN(23) [pw::Connector create]
  $_CN(23) addSegment $_TMP(PW_1)
  $_CN(23) calculateDimension
  unset _TMP(PW_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentCircle create]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(split_params) [list]
lappend _TMP(split_params) [$_CN(22) getParameter -Y 0]
set _TMP(PW_1) [$_CN(22) split $_TMP(split_params)]
unset _TMP(PW_1)
unset _TMP(split_params)
pw::Application markUndoLevel Split

set _TMP(split_params) [list]
lappend _TMP(split_params) [$_CN(23) getParameter -Y 0]
set _TMP(PW_1) [$_CN(23) split $_TMP(split_params)]
unset _TMP(PW_1)
unset _TMP(split_params)
pw::Application markUndoLevel Split


# extrude to more of a circle
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::Edge createFromConnectors [list $_CN(5) $_CN(8) $_CN(18) $_CN(20)]]
  set _TMP(edge_1) [lindex $_TMP(PW_1) 0]
  unset _TMP(PW_1)
  set _DM(7) [pw::DomainStructured create]
  $_DM(7) addEdge $_TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_DM(7)]]
  $_TMP(mode_1) setKeepFailingStep true
  $_DM(7) setExtrusionSolverAttribute NormalMarchingVector {-0 -0 -1}
  $_DM(7) setExtrusionSolverAttribute NormalInitialStepSize 0.004
  $_DM(7) setExtrusionSolverAttribute StopAtHeight Off
  $_DM(7) setExtrusionSolverAttribute StopAtHeight 50.0
  $_DM(7) setExtrusionSolverAttribute NormalExplicitSmoothing 0.5
  $_DM(7) setExtrusionSolverAttribute NormalImplicitSmoothing 1.0
  $_DM(7) setExtrusionSolverAttribute NormalKinseyBarthSmoothing 0.1
  $_DM(7) setExtrusionSolverAttribute NormalVolumeSmoothing 0.5
  $_TMP(mode_1) run 1000
  pw::Display setShowDomains 1
$_TMP(mode_1) end
unset _TMP(mode_1)
unset _TMP(edge_1)
pw::Application markUndoLevel {Extrude, Normal}


# merge and smooth
set _TMP(face_1) [pw::FaceStructured create]
$_TMP(face_1) delete
unset _TMP(face_1)
set _TMP(PW_2) [pw::DomainStructured join -reject _TMP(ignored) [list $_DM(4) $_DM(2)]]
unset _TMP(ignored)
unset _TMP(PW_2)
pw::Application markUndoLevel Join

set _DM(8) [pw::GridEntity getByName dom-2]
set _TMP(face_1) [pw::FaceStructured create]
$_TMP(face_1) delete
unset _TMP(face_1)
set _TMP(PW_2) [pw::DomainStructured join -reject _TMP(ignored) [list $_DM(8) $_DM(5)]]
unset _TMP(ignored)
unset _TMP(PW_2)
pw::Application markUndoLevel Join



# more box spacing
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(12)]]
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(10)]]
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(12)]]
  $_CN(12) replaceDistribution 1 [pw::DistributionTanh create]
  [$_CN(12) getDistribution 1] setBeginSpacing 0.0
  [$_CN(12) getDistribution 1] setEndSpacing 0.0
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(12)]]
  [[$_CN(12) getDistribution 1] getBeginSpacing] setValue 0.0003
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(12)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(12)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(10)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(10)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(10)]]
  [[$_CN(10) getDistribution 1] getBeginSpacing] setValue 0.0002
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(10)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

# smooth
set _DM(9) [pw::GridEntity getByName dom-5]
set _TMP(face_1) [pw::FaceStructured create]
$_TMP(face_1) addDomain $_DM(9)
$_TMP(face_1) addDomain -linkage [list 1 1 1 1 3 1 0] $_DM(6)
set _TMP(PW_2) [$_TMP(face_1) joinDomains]
unset _TMP(PW_2)
$_TMP(face_1) delete
unset _TMP(face_1)
pw::Application markUndoLevel Join

set _DM(10) [pw::GridEntity getByName dom-5]
set _TMP(face_1) [pw::FaceStructured create]
$_TMP(face_1) delete
unset _TMP(face_1)
set _TMP(PW_2) [pw::DomainStructured join -reject _TMP(ignored) [list $_DM(10) $_DM(7)]]
unset _TMP(ignored)
unset _TMP(PW_2)
pw::Application markUndoLevel Join

set _DM(11) [pw::GridEntity getByName dom-7]
set _TMP(mode_1) [pw::Application begin EllipticSolver [list $_DM(11)]]
  set _TMP(SUB_1) [$_DM(11) getSubGrid 1]
  set _TMP(SUB_2) [$_DM(11) getSubGrid 2]
  set _TMP(SUB_3) [$_DM(11) getSubGrid 3]
  set _TMP(SUB_4) [$_DM(11) getSubGrid 4]
  set _TMP(SUB_5) [$_DM(11) getSubGrid 5]
  set _TMP(SUB_6) [$_DM(11) getSubGrid 6]
  set _TMP(SUB_7) [$_DM(11) getSubGrid 7]
  set _TMP(SUB_8) [$_DM(11) getSubGrid 8]
  foreach ent [list $_DM(11)] bc [list 2] {
  $ent setEllipticSolverAttribute -edge $bc EdgeControl StegerSorenson
}
  foreach ent [list $_DM(11)] bc [list 2] {
#  $ent setEllipticSolverAttribute -edge $bc EdgeConstraint Floating
  $ent setEllipticSolverAttribute -edge $bc EdgeConstraint Fixed
}
#  $_DM(11) setEllipticSolverAttribute InteriorControl Laplace
#  $_DM(11) setEllipticSolverAttribute InteriorControl ThomasMiddlecoff
  $_DM(11) setEllipticSolverAttribute InteriorControl Fixed
  $_TMP(mode_1) setActiveSubGrids $_DM(11) [list]
  $_TMP(mode_1) run 100
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Solve

set _TMP(mode_1) [pw::Application begin EllipticSolver [list $_DM(11)]]
$_TMP(mode_1) abort
unset _TMP(mode_1)
unset _TMP(SUB_2)
unset _TMP(SUB_3)
unset _TMP(SUB_4)
unset _TMP(SUB_5)
unset _TMP(SUB_6)
unset _TMP(SUB_7)
unset _TMP(SUB_8)
unset _TMP(SUB_1)

# split for bc's
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::GridShape create]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _CN(24) [pw::GridEntity getByName con-16-split-2]
set _CN(25) [pw::GridEntity getByName con-17-split-2]
pw::Entity delete [list $_CN(22) $_CN(23) $_CN(24) $_CN(25)]
pw::Application markUndoLevel Delete

# fine tuning may be required here...
set _CN(26) [pw::GridEntity getByName con-17]
set _TMP(split_params) [list]
#lappend _TMP(split_params) [$_CN(26) getParameter -closest [pw::Application getXYZ [$_CN(26) getXYZ -parameter 0.215]]]
#lappend _TMP(split_params) [$_CN(26) getParameter -closest [pw::Application getXYZ [$_CN(26) getXYZ -parameter 0.715]]]
lappend _TMP(split_params) [$_CN(26) getParameter -closest [pw::Application getXYZ [$_CN(26) getXYZ -parameter 0.1825]]]
lappend _TMP(split_params) [$_CN(26) getParameter -closest [pw::Application getXYZ [$_CN(26) getXYZ -parameter 0.6825]]]
set _TMP(PW_1) [$_CN(26) split $_TMP(split_params)]
unset _TMP(PW_1)
unset _TMP(split_params)
pw::Application markUndoLevel Split


# Extrude, jacobian evals seem off
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::FaceStructured createFromDomains [list $_DM(6) $_DM(1) $_DM(3)]]
  set _TMP(face_1) [lindex $_TMP(PW_1) 0]
  set _TMP(face_2) [lindex $_TMP(PW_1) 1]
  set _TMP(face_3) [lindex $_TMP(PW_1) 2]
  unset _TMP(PW_1)
  set _BL(1) [pw::BlockStructured create]
  $_BL(1) addFace $_TMP(face_1)
  set _BL(2) [pw::BlockStructured create]
  $_BL(2) addFace $_TMP(face_2)
  set _BL(3) [pw::BlockStructured create]
  $_BL(3) addFace $_TMP(face_3)
$_TMP(mode_1) end
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_BL(1) $_BL(2) $_BL(3)]]
  $_TMP(mode_1) setKeepFailingStep true
  $_BL(1) setExtrusionSolverAttribute Mode Translate
  $_BL(2) setExtrusionSolverAttribute Mode Translate
  $_BL(3) setExtrusionSolverAttribute Mode Translate
  $_BL(1) setExtrusionSolverAttribute TranslateDirection {1 0 0}
  $_BL(2) setExtrusionSolverAttribute TranslateDirection {1 0 0}
  $_BL(3) setExtrusionSolverAttribute TranslateDirection {1 0 0}
  $_BL(1) setExtrusionSolverAttribute StopAtPositiveSkewJacobian false
  $_BL(2) setExtrusionSolverAttribute StopAtPositiveSkewJacobian false
  $_BL(3) setExtrusionSolverAttribute StopAtPositiveSkewJacobian false
  $_BL(1) setExtrusionSolverAttribute StopAtZeroJacobian false
  $_BL(2) setExtrusionSolverAttribute StopAtZeroJacobian false
  $_BL(3) setExtrusionSolverAttribute StopAtZeroJacobian false
  $_BL(1) setExtrusionSolverAttribute StopAtNegativeJacobian false
  $_BL(2) setExtrusionSolverAttribute StopAtNegativeJacobian false
  $_BL(3) setExtrusionSolverAttribute StopAtNegativeJacobian false
  $_BL(1) setExtrusionSolverAttribute StopAtNegativeSkewJacobian false
  $_BL(2) setExtrusionSolverAttribute StopAtNegativeSkewJacobian false
  $_BL(3) setExtrusionSolverAttribute StopAtNegativeSkewJacobian false
  $_BL(1) setExtrusionSolverAttribute TranslateDirection {0 0 1}
  $_BL(2) setExtrusionSolverAttribute TranslateDirection {0 0 1}
  $_BL(3) setExtrusionSolverAttribute TranslateDirection {0 0 1}
  $_BL(1) setExtrusionSolverAttribute TranslateDistance $spanlength
  $_BL(2) setExtrusionSolverAttribute TranslateDistance $spanlength
  $_BL(3) setExtrusionSolverAttribute TranslateDistance $spanlength
  $_TMP(mode_1) run $spansteps
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Extrude, Translate}

unset _TMP(face_3)
unset _TMP(face_2)
unset _TMP(face_1)


# tag bc's
set _TMP(PW_1) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

unset _TMP(PW_1)
set _TMP(PW_1) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

unset _TMP(PW_1)
set _TMP(PW_1) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

unset _TMP(PW_1)
set _TMP(PW_1) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

unset _TMP(PW_1)
set _TMP(PW_1) [pw::BoundaryCondition getByName bc-2]
$_TMP(PW_1) setName inlet
pw::Application markUndoLevel {Name BC}

set _TMP(PW_2) [pw::BoundaryCondition getByName bc-3]
$_TMP(PW_2) setName outlet
pw::Application markUndoLevel {Name BC}

set _TMP(PW_3) [pw::BoundaryCondition getByName bc-4]
$_TMP(PW_3) setName wall
pw::Application markUndoLevel {Name BC}

set _TMP(PW_4) [pw::BoundaryCondition getByName bc-5]
$_TMP(PW_4) setName periodic
pw::Application markUndoLevel {Name BC}

set _DM(12) [pw::GridEntity getByName dom-9]
set _DM(13) [pw::GridEntity getByName dom-23]
set _DM(14) [pw::GridEntity getByName dom-22]
set _DM(15) [pw::GridEntity getByName dom-15]
set _DM(16) [pw::GridEntity getByName dom-10]
set _DM(17) [pw::GridEntity getByName dom-24]
set _DM(18) [pw::GridEntity getByName dom-27]
set _DM(19) [pw::GridEntity getByName dom-33]
$_TMP(PW_3) apply [list [list $_BL(2) $_DM(6)] [list $_BL(1) $_DM(12)] [list $_BL(2) $_DM(13)] [list $_BL(1) $_DM(1)] [list $_BL(3) $_DM(3)] [list $_BL(2) $_DM(14)] [list $_BL(1) $_DM(15)] [list $_BL(1) $_DM(16)] [list $_BL(2) $_DM(17)] [list $_BL(2) $_DM(18)] [list $_BL(3) $_DM(19)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_4) apply [list [list $_BL(2) $_DM(6)] [list $_BL(2) $_DM(18)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_1) apply [list [list $_BL(2) $_DM(13)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_2) apply [list [list $_BL(2) $_DM(14)] [list $_BL(2) $_DM(17)]]
pw::Application markUndoLevel {Set BC}

unset _TMP(PW_1)
unset _TMP(PW_2)
unset _TMP(PW_3)
unset _TMP(PW_4)

# fini...?

# Set volume conditions

set _TMP(PW_1) [pw::VolumeCondition create]
pw::Application markUndoLevel {Create VC}

$_TMP(PW_1) setName Flow
pw::Application markUndoLevel {Name VC}

$_TMP(PW_1) apply [list $_BL(2)]
pw::Application markUndoLevel {Set VC}

set _DM(20) [pw::GridEntity getByName dom-21]
set _DM(21) [pw::GridEntity getByName dom-20]
set _DM(22) [pw::GridEntity getByName dom-16]
set _DM(23) [pw::GridEntity getByName dom-17]
set _DM(24) [pw::GridEntity getByName dom-18]
set _DM(25) [pw::GridEntity getByName dom-13]
$_TMP(PW_1) apply [list $_BL(1)]
pw::Application markUndoLevel {Set VC}

set _DM(26) [pw::GridEntity getByName dom-14]
set _DM(27) [pw::GridEntity getByName dom-8]
set _DM(28) [pw::GridEntity getByName dom-12]
$_TMP(PW_1) apply [list $_BL(3)]
pw::Application markUndoLevel {Set VC}

unset _TMP(PW_1)

# Boundary conditions

set _TMP(PW_1) [pw::BoundaryCondition getByName inlet]
$_TMP(PW_1) setName inflow
pw::Application markUndoLevel {Name BC}

set _TMP(PW_2) [pw::BoundaryCondition getByName outlet]
$_TMP(PW_2) setName outflow
pw::Application markUndoLevel {Name BC}

set _TMP(PW_3) [pw::BoundaryCondition getByName wall]
$_TMP(PW_3) setName airfoil
pw::Application markUndoLevel {Name BC}

set _TMP(PW_4) [pw::BoundaryCondition getByName periodic]
$_TMP(PW_4) setName front
pw::Application markUndoLevel {Name BC}

set _TMP(PW_5) [pw::BoundaryCondition create]
pw::Application markUndoLevel {Create BC}

unset _TMP(PW_5)
set _TMP(PW_5) [pw::BoundaryCondition getByName bc-6]
$_TMP(PW_5) setName back
pw::Application markUndoLevel {Name BC}

set _TMP(PW_6) [pw::BoundaryCondition getByName Unspecified]
$_TMP(PW_6) apply [list [list $_BL(1) $_DM(12)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_5) apply [list [list $_BL(2) $_DM(7)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_4) apply [list [list $_BL(1) $_DM(15)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_5) apply [list [list $_BL(1) $_DM(1)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_4) apply [list [list $_BL(3) $_DM(19)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_5) apply [list [list $_BL(3) $_DM(3)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_3) apply [list [list $_BL(1) $_DM(12)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_1) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_2) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_3) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_4) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

$_TMP(PW_5) setPhysicalType -usage CAE {Side Set}
pw::Application markUndoLevel {Change BC Type}

unset _TMP(PW_6)
unset _TMP(PW_1)
unset _TMP(PW_2)
unset _TMP(PW_3)
unset _TMP(PW_4)
unset _TMP(PW_5)

# Save the mesh

set _TMP(mode_1) [pw::Application begin CaeExport [pw::Entity sort [list $_BL(1) $_BL(2) $_BL(3)]]]
  $_TMP(mode_1) initialize -strict -type CAE $dir/CO2_lesfoil-$ref.exo
  $_TMP(mode_1) verify
  $_TMP(mode_1) write
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application save $dir/CO2_lesfoil-$ref.pw