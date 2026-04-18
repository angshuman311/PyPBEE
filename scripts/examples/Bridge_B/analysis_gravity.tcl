# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

if {$write_model_files == 1} {
	set GravityAnalysisResultsDirectory "GravityAnalysisResults";
	file mkdir $GravityAnalysisResultsDirectory;

	for {set i 1} {$i <= [llength $ColumnElementList]} {incr i} {
		set ColElements [lindex $ColumnElementList [expr $i - 1]];

		set EleTag [lindex $ColElements 0];
		set colElemForceFile "colElemForce_col_$i\_elem_$EleTag.txt"
		recorder Element -file $GravityAnalysisResultsDirectory/$colElemForceFile -time -ele $EleTag localForce

		if {[llength $ColElements] > 1} {
			set EleTag [lindex $ColElements [expr [llength $ColElements] - 1]];
			set colElemForceFile "colElemForce_col_$i\_elem_$EleTag.txt"
			recorder Element -file $GravityAnalysisResultsDirectory/$colElemForceFile -time -ele $EleTag localForce
		}
	}
}

set nStepsGrav 10;
test EnergyIncr 1e-20 2000 0;
algorithm KrylovNewton;
integrator LoadControl [expr 1./$nStepsGrav];
numberer RCM;
constraints Transformation
system ProfileSPD;
analysis Static;
set ok [analyze [expr $nStepsGrav]];
loadConst -time 0.0;

if {$ok == 0} {
	puts "Gravity Analysis COMPLETE!"
	set GravityAnalysisDone "YES";
} else {
	puts "Gravity Analysis FAILED!"
}

if {$write_model_files == 1} {
	remove recorders;
	for {set i 1} {$i <= [llength $ColumnElementList]} {incr i} {
		set ColElements [lindex $ColumnElementList [expr $i - 1]];

		set EleTag [lindex $ColElements 0];
		set colElemForceFile "colElemForce_col_$i\_elem_$EleTag.txt"

		set colElemForceFileID [open $GravityAnalysisResultsDirectory/$colElemForceFile "r"];
		set colElemForceCurr [read $colElemForceFileID];
		close $colElemForceFileID;
		file delete "$GravityAnalysisResultsDirectory/$colElemForceFile";

		set colElemForceCurr [split $colElemForceCurr "\n"];
		set colElemForceCurr [lindex $colElemForceCurr [expr $nStepsGrav - 1]];
		set colElemAxialForceCurr [lindex $colElemForceCurr 1];
		set colElemAxialForceCurr [expr abs($colElemAxialForceCurr)];

		set predFileID1 [open "$model_info_directory/predictor_info_col_rebar_strain_damage_col_$i\_edge_1.txt" "a"];
		set predFileID2 [open "$model_info_directory/predictor_info_col_rebar_strain_damage_col_$i\_edge_2.txt" "a"];
		puts $predFileID1 "$colElemAxialForceCurr";
		close $predFileID1;

		if {[llength $ColElements] > 1} {
			set EleTag [lindex $ColElements [expr [llength $ColElements] - 1]];
			set colElemForceFile "colElemForce_col_$i\_elem_$EleTag.txt"

			set colElemForceFileID [open $GravityAnalysisResultsDirectory/$colElemForceFile "r"];
			set colElemForceCurr [read $colElemForceFileID];
			close $colElemForceFileID;
			file delete "$GravityAnalysisResultsDirectory/$colElemForceFile";

			set colElemForceCurr [split $colElemForceCurr "\n"];
			set colElemForceCurr [lindex $colElemForceCurr [expr $nStepsGrav - 1]];
			set colElemAxialForceCurr [lindex $colElemForceCurr 1];
			set colElemAxialForceCurr [expr abs($colElemAxialForceCurr)];
		}

		puts $predFileID2 "$colElemAxialForceCurr";
		close $predFileID2;
	}
	file delete -force -- $GravityAnalysisResultsDirectory
}
wipeAnalysis
