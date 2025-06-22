# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

source "$model_files_path/nltha_params.tcl"

set addZerosFor 		[expr $addZerosFor*1.0];
set dtAnalysis 			[expr $dtAnalysis*1.0];
set gmSkew 				[expr $gmSkew*1.0];
set maxDimKrylov 		[expr int($maxDimKrylov)];
set maxNumIterDyn 		[expr int($maxNumIterDyn)];
set maxNumIterDynBasic 	[expr int($maxNumIterDynBasic)];
set scaleFac_L 			[expr $scaleFac_L*1.0];
set scaleFac_T 			[expr $scaleFac_T*1.0];
set scaleFac_V 			[expr $scaleFac_V*1.0];
set showTest 			[expr int($showTest)];
set showTestBasic 		[expr int($showTestBasic)];
set tolDynBasic 		[expr $tolDynBasic*1.0];
set tolDynDisp 			[expr $tolDynDisp*1.0];
set tolDynUnb 			[expr $tolDynUnb*1.0];

set gm_dirctn_L 1;			# ground-motion direction
set gm_dirctn_T 2;			# ground-motion direction
set gm_dirctn_V 3;			# ground-motion direction

# set up ground-motion-analysis parameters
set gmDataFile_L [open "gm_1.txt" "r"];
set gmDataFile_T [open "gm_2.txt" "r"];
set gmDataFile_V [open "gm_3.txt" "r"];

set gmData_L [read $gmDataFile_L];
set gmData_T [read $gmDataFile_T];
set gmData_V [read $gmDataFile_V];

close $gmDataFile_L;
close $gmDataFile_T;
close $gmDataFile_V;

set npts_L [lindex $gmData_L 0]
set npts_T [lindex $gmData_T 0]
set npts_V [lindex $gmData_V 0]

set dtRec_L [lindex $gmData_L 1];
set dtRec_T [lindex $gmData_T 1];
set dtRec_V [lindex $gmData_V 1];

set gmInput_L {};
for {set i 1} {$i <= $npts_L} {incr i} {
lappend gmInput_L [expr [expr [lindex $gmData_L [expr $i + 1]]*cos($gmSkew*$pi/180)] + [expr -[lindex $gmData_T [expr $i + 1]]*sin($gmSkew*$pi/180)]];
}
for {set i 1} {$i <= [expr int($addZerosFor/$dtRec_L)]} {incr i 1} {
lappend gmInput_L 0.;
}

set gmInput_T {};
for {set i 1} {$i <= $npts_T} {incr i} {
lappend gmInput_T [expr [expr [lindex $gmData_L [expr $i + 1]]*sin($gmSkew*$pi/180)] + [expr [lindex $gmData_T [expr $i + 1]]*cos($gmSkew*$pi/180)]];
}
for {set i 1} {$i <= [expr int($addZerosFor/$dtRec_T)]} {incr i 1} {
lappend gmInput_T 0.;
}

set gmInput_V {};
for {set i 1} {$i <= $npts_V} {incr i} {
lappend gmInput_V [lindex $gmData_V [expr $i + 1]];
}
for {set i 1} {$i <= [expr int($addZerosFor/$dtRec_V)]} {incr i 1} {
lappend gmInput_V 0.;
}

set totalAnalysisTimes	[list [expr ($npts_L + int($addZerosFor/$dtRec_L))*$dtRec_L] [expr ($npts_T + int($addZerosFor/$dtRec_T))*$dtRec_T] [expr ($npts_V + int($addZerosFor/$dtRec_V))*$dtRec_V]];

set totalAnalysisTime [lindex $totalAnalysisTimes 0]
for {set i 1} {$i < 3} {incr i 1} {
    if {[lindex $totalAnalysisTimes $i] < $totalAnalysisTime} {
        set totalAnalysisTime [lindex $totalAnalysisTimes $i]
    }
}

# ###################################################################################################################################################################################################
# ################################################################################# SET LOAD PATTERN ################################################################################################
# ###################################################################################################################################################################################################

#  perform Dynamic Ground-Motion Analysis
set gmLoadTag_L 2;	# LoadTag for Uniform ground motion excitation
set gmLoadTag_T 3;	# LoadTag for Uniform ground motion excitation
set gmLoadTag_V 4;	# LoadTag for Uniform ground motion excitation

set gmFact_L [expr $g*$scaleFac_L];		# data in input file is in g Unifts -- ACCELERATION TH
set gmFact_T [expr -$g*$scaleFac_T];		# data in input file is in g Unifts -- ACCELERATION TH
set gmFact_V [expr $g*$scaleFac_V];		# data in input file is in g Unifts -- ACCELERATION TH

set tsTag_L $gmLoadTag_L;
set tsTag_T $gmLoadTag_T;
set tsTag_V $gmLoadTag_V;

timeSeries Path $tsTag_L -dt $dtRec_L -values $gmInput_L -factor $gmFact_L;
timeSeries Path $tsTag_T -dt $dtRec_T -values $gmInput_T -factor $gmFact_T;
if {$vert == 1} {
    timeSeries Path $tsTag_V -dt $dtRec_V -values $gmInput_V -factor $gmFact_V;
}

pattern UniformExcitation  $gmLoadTag_L  $gm_dirctn_L -accel  $tsTag_L;		# create Unifform excitation
pattern UniformExcitation  $gmLoadTag_T  $gm_dirctn_T -accel  $tsTag_T;		# create Unifform excitation
if {$vert == 1} {
    pattern UniformExcitation  $gmLoadTag_V  $gm_dirctn_V -accel  $tsTag_V;		# create Unifform excitation
}

# ###################################################################################################################################################################################################
# ############################################################################### SET ANALYSIS PARAMETERS ###########################################################################################
# ###################################################################################################################################################################################################

set NewmarkGamma 0.50;	# Newmark-integrator gamma parameter (also HHT)
set NewmarkBeta  0.25;	# Newmark-integrator beta parameter
test $testBasic $tolDynBasic $maxNumIterDynBasic $showTestBasic;
algorithm {*}$algorithmBasic;
integrator Newmark $NewmarkGamma $NewmarkBeta
# integrator TRBDF2
numberer RCM;
constraints Transformation
system UmfPack;
# analysis VariableTransient -numSublevels 10 -numSubSteps 5
analysis Transient

# ###################################################################################################################################################################################################
# ############################################################################### RUN ANALYSIS ######################################################################################################
# ###################################################################################################################################################################################################

# set analysisLogFile [open "analysis_log.txt" "w"];

set nSteps [expr int($totalAnalysisTime/$dtAnalysis)];
set TIME_start [clock clicks -milliseconds]
set ok [analyze $nSteps $dtAnalysis];
set TIME_taken [expr [clock clicks -milliseconds] - $TIME_start]
if {$ok == 0} {
	set tCurrent [getTime];
	# puts $analysisLogFile "TIME: $tCurrent ($totalAnalysisTime) >> CONVERGED!"
} else {
	set tCurrent [getTime];
	set ok 0;
	while {$ok == 0 && $tCurrent <= $totalAnalysisTime && [expr $TIME_taken*2.77778e-7] <= 1.0} {
		set TIME_taken [expr [clock clicks -milliseconds] - $TIME_start]
        set ok [analyze 1 $dtAnalysis];
		# if {$ok == 0} { puts $analysisLogFile "TIME: $tCurrent ($totalAnalysisTime) >> CONVERGED!" }
		#####################################################################################################
		############################################# DISP INCR #############################################
		#####################################################################################################
		if {$ok != 0} {
			# puts $analysisLogFile "Try Krylov Newton with DispIncr"
			test NormDispIncr $tolDynDisp $maxNumIterDyn $showTest;
			algorithm KrylovNewton -maxDim $maxDimKrylov;
			set ok [analyze 1 $dtAnalysis];
			if {$ok == 0} {
				test $testBasic $tolDynBasic $maxNumIterDynBasic $showTestBasic;
				algorithm {*}$algorithmBasic;
			}
		};
		if {$ok != 0} {
			# puts $analysisLogFile "Try Newton initial with DispIncr"
			test NormDispIncr $tolDynDisp $maxNumIterDyn $showTest;
			algorithm Newton -initial
			set ok [analyze 1 $dtAnalysis];
			if {$ok == 0} {
				test $testBasic $tolDynBasic $maxNumIterDynBasic $showTestBasic;
				algorithm {*}$algorithmBasic;
			}
		};
        if {$ok != 0} {
			# puts $analysisLogFile "Try Newton -Hall 50-50 with Unbalance"
			test NormDispIncr $tolDynDisp $maxNumIterDyn $showTest;
			algorithm Newton -Hall 0.1 0.9;
			set ok [analyze 1 $dtAnalysis];
			if {$ok == 0} {
				test $testBasic $tolDynBasic $maxNumIterDynBasic $showTestBasic;
				algorithm {*}$algorithmBasic;
			}
		};
		set tCurrent [getTime];
	}
}

# close $analysisLogFile;

if {$ok == 0} {
    set tCurrent [getTime];
    if {[expr $tCurrent/$totalAnalysisTime] >= 0.80} {
        puts "NLTHA COMPLETE!"
        set NLTHA_OK 1;
    } else {
        puts "NLTHA FAILED!"
        set NLTHA_OK 0;
    }
} else {
    puts "NLTHA FAILED!"
    set NLTHA_OK 0;
}
wipeAnalysis
