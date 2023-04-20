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


# Number of points on pressure or suction side, np_suction over-written
   set np_pressure 301   
   set np_suction 801 
   set np_s1 100
   set np_s2 100
   set np_s3 550
   set np_te 50
   set spandx 0.0025

   # Trailing and leading edge spacing, te_dx over-written
   set te_dx 0.0004
   set te_dx1 0.001
   set te_dx2 0.0004
   #set le_dx 0.0003
   set le_dx 0.0009

# separtion pt spacing
#set d_sep 0.0002
set d_sep 0.0006

# 1st to 1nd layer af mesh size
set d_junction 0.0005

# Growth factor for initial wall normal extrusion
set gf 1.08
set gf2 1.05

# Height for wall normal extrusion
set bl_dist 120.0

# First cell height (~1.2E-5)
#set ds 0.000011813977015662547
set ds 0.5e-5

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
set branch_offset {0.0002 0 0}

# discretization of "branch" => moved to mesh specs
#set branch_number 51

# Angle of attack
set angle -13.3

# outer params
set N_wake 100
set d_wake1 0.00025
set d_wake2 0.0004
set d_in 0.0008



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

#set _DB(1) [pw::DatabaseEntity getByName curve-2]
#set _DB(2) [pw::DatabaseEntity getByName curve-1]
#set _TMP(mode_1) [pw::Application begin Modify [list $_DB(1) $_DB(2)]]
#  pw::Entity transform [pwu::Transform rotation -anchor {0 0 0} {0 0 1} $angle] [$_TMP(mode_1) getEntities]
#$_TMP(mode_1) end
#unset _TMP(mode_1)
#pw::Application markUndoLevel Rotate

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


### big redo from here on

# split sep region up for more control
set _TMP(split_params) [list]
lappend _TMP(split_params) [$_CN(1) getParameter -arc [expr {0.01 * 10}]]
lappend _TMP(split_params) [$_CN(1) getParameter -arc [expr {0.01 * 16}]]
set _TMP(PW_1) [$_CN(1) split $_TMP(split_params)]
unset _TMP(PW_1)
unset _TMP(split_params)
pw::Application markUndoLevel Split


# split for T-E control
set _CN(3) [pw::GridEntity getByName con-2-split-3]
set _TMP(split_params) [list]
lappend _TMP(split_params) [$_CN(3) getParameter -X 0.94999999999999996]
set _TMP(PW_1) [$_CN(3) split $_TMP(split_params)]
unset _TMP(PW_1)
unset _TMP(split_params)
pw::Application markUndoLevel Split

set _TMP(split_params) [list]
lappend _TMP(split_params) [$_CN(2) getParameter -X 0.94999999999999996]
set _TMP(PW_1) [$_CN(2) split $_TMP(split_params)]
unset _TMP(PW_1)
unset _TMP(split_params)
pw::Application markUndoLevel Split

# re-dim
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(1)]]
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(1)]
  $_TMP(PW_1) do setDimension -resetDistribution $np_s1
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Dimension]
  set _CN(4) [pw::GridEntity getByName con-2-split-2]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(4)]
  $_TMP(PW_1) do setDimension -resetDistribution $np_s2
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(3)]
  $_TMP(PW_1) do setDimension -resetDistribution $np_s3
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Dimension]
  set _CN(5) [pw::GridEntity getByName con-1-split-2]
  set _CN(6) [pw::GridEntity getByName con-2-split-3-split-2]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(6)]
  $_TMP(PW_1) do setDimension -resetDistribution $np_te
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(2)]
  $_TMP(PW_1) do setDimension -resetDistribution $np_te
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


# re-distribute
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(1)]]
  [[$_CN(1) getDistribution 1] getBeginSpacing] setValue $le_dx
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(1)]]
  [[$_CN(1) getDistribution 1] getEndSpacing] setValue $d_sep
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(1)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(4)]]
  [[$_CN(4) getDistribution 1] getBeginSpacing] setValue $d_sep
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(4)]]
  [[$_CN(4) getDistribution 1] getEndSpacing] setValue $d_sep
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(4)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(3)]]
  [[$_CN(3) getDistribution 1] getBeginSpacing] setValue $d_sep
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(3)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(3)]]
  [[$_CN(3) getDistribution 1] getEndSpacing] setValue $te_dx1
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(3)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(6)]]
  [[$_CN(6) getDistribution 1] getBeginSpacing] setValue $te_dx1
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(6)]]
  [[$_CN(6) getDistribution 1] getEndSpacing] setValue $te_dx2
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(6)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(2)]]
  [[$_CN(2) getDistribution 1] getBeginSpacing] setValue $te_dx2
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(2)]]
  [[$_CN(2) getDistribution 1] getEndSpacing] setValue $te_dx1
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(2)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(5)]]
  [[$_CN(5) getDistribution 1] getBeginSpacing] setValue $te_dx1
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(5)]]
  [[$_CN(5) getDistribution 1] getEndSpacing] setValue $le_dx
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(5)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute


# extrude first layer
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::Edge createFromConnectors [list $_CN(5) $_CN(3) $_CN(4) $_CN(1)]]
  set _TMP(edge_1) [lindex $_TMP(PW_1) 0]
  unset _TMP(PW_1)
  set _DM(1) [pw::DomainStructured create]
  $_DM(1) addEdge $_TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_DM(1)]]
  $_TMP(mode_1) setKeepFailingStep true
  $_DM(1) setExtrusionSolverAttribute NormalInitialStepSize $ds
  $_DM(1) setExtrusionSolverAttribute SpacingGrowthFactor 1.1
  $_DM(1) setExtrusionSolverAttribute NormalMarchingVector {-0 -0 -1}
  $_DM(1) setExtrusionSolverAttribute StopAtHeight Off
  $_DM(1) setExtrusionSolverAttribute StopAtHeight 0.005
  $_DM(1) setExtrusionSolverAttribute Mode NormalHyperbolic
  $_DM(1) setExtrusionSolverAttribute NormalExplicitSmoothing 0.51
  $_DM(1) setExtrusionSolverAttribute NormalImplicitSmoothing 0.99
  $_DM(1) setExtrusionSolverAttribute NormalKinseyBarthSmoothing 0.01
  $_DM(1) setExtrusionSolverAttribute NormalVolumeSmoothing 0.01
  $_DM(1) setExtrusionBoundaryCondition Begin Splay 0.0
  $_DM(1) setExtrusionBoundaryConditionStepSuppression Begin 0
  $_DM(1) setExtrusionBoundaryCondition End Splay 0.0
  $_DM(1) setExtrusionBoundaryConditionStepSuppression End 0
  $_TMP(mode_1) run 500
  $_DM(1) setExtrusionSolverAttribute StopAtHeight 0.02
  $_DM(1) setExtrusionSolverAttribute NormalVolumeSmoothing 0.01
  $_DM(1) setExtrusionSolverAttribute SpacingGrowthFactor 1.0
  $_TMP(mode_1) run 500  
$_TMP(mode_1) end
unset _TMP(mode_1)
unset _TMP(edge_1)
pw::Application markUndoLevel {Extrude, Normal}


# build trailing edge box
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint [$_CN(6) getPosition -arc 1]
  $_TMP(PW_1) addPoint {1.006 0.017 0}
  set _CN(7) [pw::Connector create]
  $_CN(7) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(7) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  set _CN(8) [pw::GridEntity getByName con-2]
  set _CN(9) [pw::GridEntity getByName con-3]
  $_TMP(PW_1) addPoint [$_CN(7) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(8) getPosition -arc 1]
  set _CN(10) [pw::Connector create]
  $_CN(10) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(10) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint [$_CN(6) getPosition -arc 1]
  $_TMP(PW_1) addPoint {1.006 -0.038 0}
  set _CN(11) [pw::Connector create]
  $_CN(11) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(11) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  set _CN(12) [pw::GridEntity getByName con-1]
  $_TMP(PW_1) addPoint [$_CN(11) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(12) getPosition -arc 1]
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


set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint [$_CN(7) getPosition -arc 1]
  $_TMP(PW_1) addPoint {1.06 0.017 0}
  set _CN(14) [pw::Connector create]
  $_CN(14) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(14) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint [$_CN(14) getPosition -arc 1]
  $_TMP(PW_1) addPoint {1.06 -0.038 0}
  set _CN(15) [pw::Connector create]
  $_CN(15) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(15) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint [$_CN(15) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(11) getPosition -arc 1]
  set _CN(16) [pw::Connector create]
  $_CN(16) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(16) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint [$_CN(14) getPosition -arc 1]
  $_TMP(PW_1) addPoint {1.1 0.07 0}
  set _CN(17) [pw::Connector create]
  $_CN(17) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(17) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) addPoint [$_CN(15) getPosition -arc 1]
  $_TMP(PW_1) addPoint {1.1 -0.085 0}
  set _CN(18) [pw::Connector create]
  $_CN(18) addSegment $_TMP(PW_1)
  unset _TMP(PW_1)
  $_CN(18) calculateDimension
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Create 2 Point Connector}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentSpline create]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentCircle create]
  $_TMP(PW_1) addPoint [$_CN(18) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(17) getPosition -arc 1]
  $_TMP(PW_1) setAngle 90 {0 0 1}
  set _CN(19) [pw::Connector create]
  $_CN(19) addSegment $_TMP(PW_1)
  $_CN(19) calculateDimension
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

# dim current box
# 
set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(10)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(6) 1]]
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

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(13)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(2) 1]]
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

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(7)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(9) 1]]
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

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(11)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(12) 1]]
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

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(7)]]
  $_CN(7) replaceDistribution 1 [pw::DistributionTanh create]
  [$_CN(7) getDistribution 1] setBeginSpacing 0.0
  [$_CN(7) getDistribution 1] setEndSpacing 0.0
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(7)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(11)]]
  $_CN(11) replaceDistribution 1 [pw::DistributionTanh create]
  [$_CN(11) getDistribution 1] setBeginSpacing 0.0
  [$_CN(11) getDistribution 1] setEndSpacing 0.0
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(11)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(10)]]
  set _TMP(dist_1) [pw::DistributionGeneral create [list [list $_CN(6) 1]]]
  # Clear spacings so the distribution will scale properly
  $_TMP(dist_1) setBeginSpacing 0
  $_TMP(dist_1) setEndSpacing 0
  $_TMP(dist_1) setVariable [[$_CN(10) getDistribution 1] getVariable]
  $_CN(10) setDistribution -lockEnds 1 $_TMP(dist_1)
  unset _TMP(dist_1)
  set _TMP(dist_1) [pw::DistributionGeneral create [list [list $_CN(6) 1]]]
  $_TMP(dist_1) reverse
  # Clear spacings so the distribution will scale properly
  $_TMP(dist_1) setBeginSpacing 0
  $_TMP(dist_1) setEndSpacing 0
  $_TMP(dist_1) setVariable [[$_CN(10) getDistribution 1] getVariable]
  $_CN(10) setDistribution -lockEnds 1 $_TMP(dist_1)
  unset _TMP(dist_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(10)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(13)]]
  set _TMP(dist_1) [pw::DistributionGeneral create [list [list $_CN(2) 1]]]
  # Clear spacings so the distribution will scale properly
  $_TMP(dist_1) setBeginSpacing 0
  $_TMP(dist_1) setEndSpacing 0
  $_TMP(dist_1) setVariable [[$_CN(13) getDistribution 1] getVariable]
  $_CN(13) setDistribution -lockEnds 1 $_TMP(dist_1)
  unset _TMP(dist_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(13)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(15)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(7) 1] [list $_CN(11) 1]]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(19)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(15) 1]]
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

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(15)]]
  $_CN(15) replaceDistribution 1 [pw::DistributionTanh create]
  [$_CN(15) getDistribution 1] setBeginSpacing 0.0
  [$_CN(15) getDistribution 1] setEndSpacing 0.0
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(15)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(19)]]
  $_CN(19) replaceDistribution 1 [pw::DistributionTanh create]
  [$_CN(19) getDistribution 1] setBeginSpacing 0.0
  [$_CN(19) getDistribution 1] setEndSpacing 0.0
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(19)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(16)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(15) 1]]
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


# mesh existing box

# Appended by Fidelity Pointwise V18.6R2 - Wed Mar 15 00:04:32 2023

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(6)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(7)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(10)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(9)
  set _DM(2) [pw::DomainStructured create]
  $_DM(2) addEdge $_TMP(edge_1)
  $_DM(2) addEdge $_TMP(edge_2)
  $_DM(2) addEdge $_TMP(edge_3)
  $_DM(2) addEdge $_TMP(edge_4)
  unset _TMP(edge_4)
  unset _TMP(edge_3)
  unset _TMP(edge_2)
  unset _TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(13)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(11)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(2)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(12)
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
  $_TMP(edge_1) addConnector $_CN(16)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(15)
  $_TMP(edge_1) delete
  $_TMP(edge_2) delete
  unset _TMP(edge_2)
  unset _TMP(edge_1)
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Create]
$_TMP(mode_1) abort
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(14)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(16) 1]]
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

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(14)]]
  $_CN(14) replaceDistribution 1 [pw::DistributionTanh create]
  [$_CN(14) getDistribution 1] setBeginSpacing 0.0
  [$_CN(14) getDistribution 1] setEndSpacing 0.0
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(14)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(16)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(15)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(14)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(7)
  $_TMP(edge_4) addConnector $_CN(11)
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
$_TMP(mode_1) abort
unset _TMP(mode_1)


# second layer af extrusion
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::Edge createFromConnectors [list $_CN(8) $_CN(13) $_CN(10)]]
  set _TMP(edge_1) [lindex $_TMP(PW_1) 0]
  unset _TMP(PW_1)
  set _DM(5) [pw::DomainStructured create]
  $_DM(5) addEdge $_TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_DM(5)]]
  $_TMP(mode_1) setKeepFailingStep true
  $_DM(5) setExtrusionSolverAttribute NormalMarchingVector {-0 -0 -1}
  $_DM(5) setExtrusionSolverAttribute NormalInitialStepSize $d_junction
  $_DM(5) setExtrusionSolverAttribute StopAtHeight Off
  $_DM(5) setExtrusionSolverAttribute StopAtHeight 0.07
  $_DM(5) setExtrusionSolverAttribute SpacingGrowthFactor 1.1  
  $_DM(5) setExtrusionSolverAttribute Mode NormalHyperbolic
  $_DM(5) setExtrusionSolverAttribute NormalExplicitSmoothing 0.51
  $_DM(5) setExtrusionSolverAttribute NormalImplicitSmoothing 0.99
  $_DM(5) setExtrusionSolverAttribute NormalKinseyBarthSmoothing 0.01
  $_DM(5) setExtrusionSolverAttribute NormalVolumeSmoothing 0.01
  $_DM(5) setExtrusionBoundaryCondition Begin ConstantX
  $_DM(5) setExtrusionBoundaryConditionStepSuppression Begin 0
  $_DM(5) setExtrusionBoundaryCondition End ConstantX
  $_DM(5) setExtrusionBoundaryConditionStepSuppression End 0
  $_TMP(mode_1) run 32
$_TMP(mode_1) end
unset _TMP(mode_1)
unset _TMP(edge_1)
pw::Application markUndoLevel {Extrude, Normal}

# create trailing box flair to extrude mesh connectors
# 
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::SegmentCircle create]
  set _CN(20) [pw::GridEntity getByName con-15]
  set _CN(21) [pw::GridEntity getByName con-16]
  $_TMP(PW_1) addPoint [$_CN(17) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(20) getPosition -arc 1]
  $_TMP(PW_1) setAngle 45 {0 0 1}
  set _CN(22) [pw::Connector create]
  $_CN(22) addSegment $_TMP(PW_1)
  $_CN(22) calculateDimension
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
  set _TMP(PW_1) [pw::SegmentCircle create]
  set _CN(23) [pw::GridEntity getByName con-14]
  $_TMP(PW_1) addPoint [$_CN(23) getPosition -arc 1]
  $_TMP(PW_1) addPoint [$_CN(18) getPosition -arc 1]
  $_TMP(PW_1) setAngle 45 {0 0 1}
  set _CN(24) [pw::Connector create]
  $_CN(24) addSegment $_TMP(PW_1)
  $_CN(24) calculateDimension
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

# copy connector mesh to new
# 
set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(17)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(21) 1]]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(18)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(17) 1]]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(22)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(14) 1]]
  $_TMP(PW_1) delete
  unset _TMP(PW_1)
  $_TMP(mode_1) balance -resetGeneralDistributions
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Dimension

set _TMP(mode_1) [pw::Application begin Dimension]
  set _TMP(PW_1) [pw::Collection create]
  $_TMP(PW_1) set [list $_CN(24)]
  $_TMP(PW_1) do setDimensionFromSubConnectors -resetDistribution [list  [list $_CN(16) 1]]
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


# mesh remain transition regions
# 
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(24)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(18)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(16)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(23)
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
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(19)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(17)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(15)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(18)
  set _DM(7) [pw::DomainStructured create]
  $_DM(7) addEdge $_TMP(edge_1)
  $_DM(7) addEdge $_TMP(edge_2)
  $_DM(7) addEdge $_TMP(edge_3)
  $_DM(7) addEdge $_TMP(edge_4)
  unset _TMP(edge_4)
  unset _TMP(edge_3)
  unset _TMP(edge_2)
  unset _TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel {Assemble Domain}

set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(edge_1) [pw::Edge create]
  $_TMP(edge_1) addConnector $_CN(17)
  set _TMP(edge_2) [pw::Edge create]
  $_TMP(edge_2) addConnector $_CN(22)
  set _TMP(edge_3) [pw::Edge create]
  $_TMP(edge_3) addConnector $_CN(21)
  set _TMP(edge_4) [pw::Edge create]
  $_TMP(edge_4) addConnector $_CN(14)
  set _DM(8) [pw::DomainStructured create]
  $_DM(8) addEdge $_TMP(edge_1)
  $_DM(8) addEdge $_TMP(edge_2)
  $_DM(8) addEdge $_TMP(edge_3)
  $_DM(8) addEdge $_TMP(edge_4)
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


# touchups to corners of transition box
# 
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(17)]]
  [[$_CN(17) getDistribution 1] getBeginSpacing] setValue $d_junction
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(17)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(18)]]
  [[$_CN(18) getDistribution 1] getBeginSpacing] setValue $d_junction
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute

set _TMP(mode_1) [pw::Application begin Modify [list $_CN(18)]]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Distribute


# THIS CONCLUDES MAIN AF MESG GENERATION
# joining, smoothing, and forming outer O region after this


# join conns and doms
set _TMP(face_1) [pw::FaceStructured create]
$_TMP(face_1) delete
unset _TMP(face_1)
set _TMP(PW_1) [pw::DomainStructured join -reject _TMP(ignored) [list $_DM(7) $_DM(8)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _DM(9) [pw::GridEntity getByName dom-7]
set _TMP(face_1) [pw::FaceStructured create]
$_TMP(face_1) delete
unset _TMP(face_1)
set _TMP(PW_1) [pw::DomainStructured join -reject _TMP(ignored) [list $_DM(9) $_DM(5)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _DM(10) [pw::GridEntity getByName dom-5]
set _TMP(face_1) [pw::FaceStructured create]
$_TMP(face_1) addDomain $_DM(6)
$_TMP(face_1) addDomain -linkage [list 2 1 2 1 4 1 0] $_DM(10)
set _TMP(PW_1) [$_TMP(face_1) joinDomains]
unset _TMP(PW_1)
$_TMP(face_1) delete
unset _TMP(face_1)
pw::Application markUndoLevel Join

set _TMP(face_1) [pw::FaceStructured create]
$_TMP(face_1) delete
unset _TMP(face_1)
set _TMP(PW_1) [pw::DomainStructured join -reject _TMP(ignored) [list $_DM(1) $_DM(2)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _DM(11) [pw::GridEntity getByName dom-1]
set _TMP(face_1) [pw::FaceStructured create]
$_TMP(face_1) delete
unset _TMP(face_1)
set _TMP(PW_1) [pw::DomainStructured join -reject _TMP(ignored) [list $_DM(11) $_DM(3)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(22) $_CN(19)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _CN(25) [pw::GridEntity getByName con-13]
set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(20) $_CN(25)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _CN(26) [pw::GridEntity getByName con-15]
set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(24) $_CN(26)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(6) $_CN(3)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _CN(27) [pw::GridEntity getByName con-2-split-3-split-1]
set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(4) $_CN(27)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _CN(28) [pw::GridEntity getByName con-2-split-2]
set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(1) $_CN(28)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _CN(29) [pw::GridEntity getByName con-2-split-1]
set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(29) $_CN(5)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _CN(30) [pw::GridEntity getByName con-1-split-2]
set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(30) $_CN(2)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(10) $_CN(8)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

set _CN(31) [pw::GridEntity getByName con-2]
set _TMP(PW_1) [pw::Connector join -reject _TMP(ignored) -keepDistribution [list $_CN(31) $_CN(13)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join

# small amount of smoothing in af layer 2
# 
set _DM(12) [pw::GridEntity getByName dom-6]
set _TMP(mode_1) [pw::Application begin EllipticSolver [list $_DM(12)]]
  set _TMP(SUB_1) [$_DM(12) getSubGrid 1]
  set _TMP(SUB_2) [$_DM(12) getSubGrid 2]
  set _TMP(SUB_3) [$_DM(12) getSubGrid 3]
  set _TMP(SUB_4) [$_DM(12) getSubGrid 4]
  set _TMP(SUB_5) [$_DM(12) getSubGrid 5]
  set _TMP(SUB_6) [$_DM(12) getSubGrid 6]
  $_DM(12) setEllipticSolverAttribute InteriorControl Laplace
  set _TMP(ENTS) [pw::Collection create]
  $_TMP(ENTS) set [list $_DM(12)]
  $_TMP(ENTS) do setInitializeMethod Standard
  $_DM(12) setEllipticSolverAttribute ShapeConstraint Free
  $_TMP(ENTS) delete
  unset _TMP(ENTS)
  foreach ent [list $_DM(12)] bc [list 1] {
  $ent setEllipticSolverAttribute -edge $bc EdgeConstraint Orthogonal
}
  $_TMP(mode_1) setActiveSubGrids $_DM(12) [list]
  $_TMP(mode_1) run -entities [list $_DM(12)] Initialize
  $_TMP(mode_1) setActiveSubGrids $_DM(12) [list]
  $_TMP(mode_1) run 5
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Solve

set _TMP(mode_1) [pw::Application begin EllipticSolver [list $_DM(12)]]
$_TMP(mode_1) abort
unset _TMP(mode_1)
unset _TMP(SUB_1)
unset _TMP(SUB_2)
unset _TMP(SUB_3)
unset _TMP(SUB_4)
unset _TMP(SUB_5)
unset _TMP(SUB_6)

# Extrude to full O
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::Edge createFromConnectors [list $_CN(24)]]
  set _TMP(edge_1) [lindex $_TMP(PW_1) 0]
  unset _TMP(PW_1)
  set _DM(9) [pw::DomainStructured create]
  $_DM(9) addEdge $_TMP(edge_1)
$_TMP(mode_1) end
unset _TMP(mode_1)
set _TMP(mode_1) [pw::Application begin ExtrusionSolver [list $_DM(9)]]
  $_TMP(mode_1) setKeepFailingStep true
  $_DM(9) setExtrusionSolverAttribute NormalInitialStepSize 0.005
  $_DM(9) setExtrusionSolverAttribute StopAtHeight Off
  $_DM(9) setExtrusionSolverAttribute StopAtHeight 20
  $_DM(9) setExtrusionSolverAttribute SpacingGrowthFactor 1.08  
  $_DM(9) setExtrusionSolverAttribute Mode NormalHyperbolic
  $_DM(9) setExtrusionSolverAttribute NormalVolumeSmoothing 0.01
  $_DM(9) setExtrusionSolverAttribute NormalKinseyBarthSmoothing 0.4
  $_DM(9) setExtrusionSolverAttribute NormalImplicitSmoothing 0.99
  $_DM(9) setExtrusionSolverAttribute NormalExplicitSmoothing 0.51
  $_TMP(mode_1) run 200
$_TMP(mode_1) end
unset _TMP(mode_1)
unset _TMP(edge_1)
pw::Application markUndoLevel {Extrude, Normal}

set _TMP(face_1) [pw::FaceStructured create]
$_TMP(face_1) delete
unset _TMP(face_1)
set _TMP(PW_1) [pw::DomainStructured join -reject _TMP(ignored) [list $_DM(8) $_DM(9)]]
unset _TMP(ignored)
unset _TMP(PW_1)
pw::Application markUndoLevel Join


# rotate by AoA
pw::Display setShowDomains 1
set _CN(19) [pw::GridEntity getByName con-20]
set _CN(20) [pw::GridEntity getByName con-19]
set _TMP(mode_1) [pw::Application begin Modify [list $_CN(11) $_CN(7) $_CN(14) $_CN(15) $_CN(18) $_CN(16) $_DM(9) $_DB(2) $_DB(1) $_DM(4) $_DM(3) $_CN(13) $_CN(6) $_CN(19) $_CN(20)]]
  pw::Entity transform [pwu::Transform rotation -anchor {0 0 0} {0 0 1} $angle] [$_TMP(mode_1) getEntities]
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application markUndoLevel Rotate


# split along y = max/min to tag bcs, might require futzing

# Appended by Fidelity Pointwise V18.6R2 - Wed Mar 15 01:01:09 2023

pw::Display setShowDomains 0
set _TMP(split_params) [list]
lappend _TMP(split_params) [$_CN(19) getParameter -closest [pw::Application getXYZ [$_CN(19) getXYZ -parameter 0.43783188003860563]]]
lappend _TMP(split_params) [$_CN(19) getParameter -closest [pw::Application getXYZ [$_CN(19) getXYZ -parameter 0.81242632730830555]]]
set _TMP(PW_1) [$_CN(19) split $_TMP(split_params)]
unset _TMP(PW_1)
unset _TMP(split_params)
pw::Application markUndoLevel Split

# extrude to 3d
# 
pw::Display setShowDomains 1
set _TMP(mode_1) [pw::Application begin Create]
  set _TMP(PW_1) [pw::FaceStructured createFromDomains [list $_DM(3) $_DM(4) $_DM(9)]]
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
pw::Display setShowDomains 0
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

set _DM(10) [pw::GridEntity getByName dom-15]
set _DM(11) [pw::GridEntity getByName dom-14]
set _DM(12) [pw::GridEntity getByName dom-13]
set _DM(13) [pw::GridEntity getByName dom-12]
set _DM(14) [pw::GridEntity getByName dom-24]
set _DM(15) [pw::GridEntity getByName dom-27]
set _DM(16) [pw::GridEntity getByName dom-29]
$_TMP(PW_4) apply [list [list $_BL(2) $_DM(10)] [list $_BL(2) $_DM(11)] [list $_BL(3) $_DM(3)] [list $_BL(2) $_DM(12)] [list $_BL(1) $_DM(4)] [list $_BL(1) $_DM(13)] [list $_BL(2) $_DM(9)] [list $_BL(2) $_DM(14)] [list $_BL(3) $_DM(15)] [list $_BL(3) $_DM(16)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_1) apply [list [list $_BL(2) $_DM(11)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_2) apply [list [list $_BL(2) $_DM(10)] [list $_BL(2) $_DM(12)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_3) apply [list [list $_BL(3) $_DM(15)]]
pw::Application markUndoLevel {Set BC}

unset _TMP(PW_1)
unset _TMP(PW_2)
unset _TMP(PW_3)
unset _TMP(PW_4)
# 
# fini
# 

# Appended by Fidelity Pointwise 2022.1.2 - Thu Apr 20 15:14:53 2023

pw::Display setShowDomains 1
set _TMP(PW_1) [pw::VolumeCondition create]
pw::Application markUndoLevel {Create VC}

$_TMP(PW_1) setName Flow
pw::Application markUndoLevel {Name VC}

$_TMP(PW_1) apply [list $_BL(1) $_BL(2) $_BL(3)]
pw::Application markUndoLevel {Set VC}

set _DM(17) [pw::GridEntity getByName dom-10]
set _DM(18) [pw::GridEntity getByName dom-11]
set _DM(19) [pw::GridEntity getByName dom-7]
set _DM(20) [pw::GridEntity getByName dom-8]
set _DM(21) [pw::GridEntity getByName dom-9]
set _DM(22) [pw::GridEntity getByName dom-17]
set _DM(23) [pw::GridEntity getByName dom-16]
set _DM(24) [pw::GridEntity getByName dom-19]
unset _TMP(PW_1)
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

$_TMP(PW_5) apply [list [list $_BL(2) $_DM(9)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_5) apply [list [list $_BL(3) $_DM(3)]]
pw::Application markUndoLevel {Set BC}

$_TMP(PW_5) apply [list [list $_BL(1) $_DM(4)]]
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

unset _TMP(PW_1)
unset _TMP(PW_2)
unset _TMP(PW_3)
unset _TMP(PW_4)
unset _TMP(PW_5)

set _TMP(mode_1) [pw::Application begin CaeExport [pw::Entity sort [list $_BL(1) $_BL(2) $_BL(3)]]]
  $_TMP(mode_1) initialize -strict -type CAE $dir/CO2_lesfoil-v5.exo
  $_TMP(mode_1) verify
  $_TMP(mode_1) write
$_TMP(mode_1) end
unset _TMP(mode_1)
pw::Application save $dir/CO2_lesfoil-v5.pw