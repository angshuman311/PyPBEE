# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

set IDGravLoadTag 1;

pattern Plain $IDGravLoadTag "Linear" {

# ##################################################################################################
if {$write_model_files == 1} {
puts $ExportID "pattern Plain $IDGravLoadTag \"Linear\" \{"
puts $ExportID "#Deck Secondary Nodes"
}
# ##################################################################################################

for {set i 0} {$i <= [expr [llength $secondaryNodesDeck] - 1]} {incr i 1} {

	set nodeID [expr [lindex $secondaryNodesDeck $i]];
	set L_trib [expr ([lindex $L [expr $nodeID/100]])/$nEleSpan];
	
	load $nodeID 0.\
				0.\
				[expr -$wconc*$Adeck*$L_trib]\
				0.\
				0\
				0;
	
	if {$write_model_files == 1} {
	puts $ExportID "load $nodeID 0.\
				0.\
				[expr -$wconc*$Adeck*$L_trib]\
				0.\
				0\
				0;"
	}
	
}




# ##################################################################################################
if {$write_model_files == 1} {
puts $ExportID "#Deck Primary Nodes"
}
# ##################################################################################################

for {set i 0} {$i <= [expr [llength $primaryNodesDeck] - 1]} {incr i 1} {

	set nodeID [expr [lindex $primaryNodesDeck $i]];
	
	if {$i == 0} {
		set L_trib [expr ([lindex $L [expr $nodeID/100]])/$nEleSpan];
		load $nodeID 0.\
					0.\
					[expr -0.5*$wconc*$Adeck*$L_trib]\
					0.\
					0\
					0;
		
		if {$write_model_files == 1} {
		puts $ExportID "load $nodeID 0.\
					0.\
					[expr -0.5*$wconc*$Adeck*$L_trib]\
					0.\
					0\
					0;"
		}
		
	} elseif {$i == $nSpans} {
		set L_trib [expr ([lindex $L end])/$nEleSpan];
		load $nodeID 0.\
					0.\
					[expr -0.5*$wconc*$Adeck*$L_trib]\
					0.\
					0.\
					0.;
		
		if {$write_model_files == 1} {		
		puts $ExportID "load $nodeID 0.\
					0.\
					[expr -0.5*$wconc*$Adeck*$L_trib]\
					0.\
					0.\
					0.;"
		}
		
	} else {
		
		set L_trib_ahead [expr ([lindex $L [expr $nodeID/100]])/$nEleSpan];
		set L_trib_behind [expr ([lindex $L [expr $nodeID/100 - 1]])/$nEleSpan];
		
		
		load $nodeID 0.\
					0.\
					[expr -(0.5*$wconc*$Adeck*$L_trib_ahead + 0.5*$wconc*$Adeck*$L_trib_behind)]\
					0.\
					0.\
					0.;
		
		if {$write_model_files == 1} {
		puts $ExportID "load $nodeID 0.\
					0.\
					[expr -(0.5*$wconc*$Adeck*$L_trib_ahead + 0.5*$wconc*$Adeck*$L_trib_behind)]\
					0.\
					0.\
					0.;";
		}
		
	}

}



set secondaryNodesColFlag [info exists secondaryNodesColumn];
if {$secondaryNodesColFlag == 1} {
	###################################################################################################
	if {$write_model_files == 1} {
	puts $ExportID "#Column Secondary Nodes"
	}
	###################################################################################################
	for {set i 0} {$i <= [expr [llength $secondaryNodesColumn] - 1]} {incr i 1} {
	
		set nodeID [expr [lindex $secondaryNodesColumn $i]];
		set Acolcurr [expr [lindex $Acol [expr $nodeID/1000 - 1]]];
		set Dcolcurr [expr [lindex $Dcol [expr $nodeID/1000 - 1]]];
		set HCol_trib [expr ([lindex $HCol [expr $nodeID/1000 - 1]])/$nEleCol];
		load $nodeID 0.\
					 0.\
					 [expr -$wconc*$Acolcurr*$HCol_trib]\
					 0.\
					 0.\
					 0.;  
		
		if {$write_model_files == 1} {
		puts $ExportID "load $nodeID 0.\
					 0.\
					 [expr -$wconc*$Acolcurr*$HCol_trib]\
					 0.\
					 0.\
					 0.;";
		}
		
	}

}

# ##################################################################################################
if {$write_model_files == 1} {
puts $ExportID "#Column Primary Nodes"
}
# ##################################################################################################

for {set i 0} {$i <= [expr [llength $primaryNodesColumn] - 1]} {incr i 1} {

	set nodeID [expr [lindex $primaryNodesColumn $i]];
	
	set Acapcurr [expr [lindex $Acap [expr $nodeID/1000 - 1]]];
	set LcapcurrR [expr [lindex $Lcap [expr 2*($nodeID/1000 - 1) - 1]]];
	set LcapcurrL [expr [lindex $Lcap [expr 2*($nodeID/1000 - 1)]]];
	set TotalLcapcurr [expr $LcapcurrL + $LcapcurrR];
	set TotalLcapcurrSkewed [expr $TotalLcapcurr/cos($skew*$pi/180.0)];
	set TotWeightCapCurr [expr 1.*$wconc*$Acapcurr*$TotalLcapcurrSkewed];
	set Acolcurr [expr [lindex $Acol [expr $nodeID/1000 - 1]]];
	set Dcolcurr [expr [lindex $Dcol [expr $nodeID/1000 - 1]]];
	set HCol_trib [expr ([lindex $HCol [expr $nodeID/1000 - 1]])/$nEleCol];
	set tempnodeID [expr $nodeID - ($nodeID/100)*100]
	if {$tempnodeID == [expr $nEleCol + 1]}  {
		
		load $nodeID 0.\
					0.\
					[expr -0.5*$wconc*$Acolcurr*$HCol_trib + (-1./$nCols)*$TotWeightCapCurr]\
					0.\
					0.\
					0.;  
		
		if {$write_model_files == 1} {
		puts $ExportID "load $nodeID 0.\
					0.\
					[expr -0.5*$wconc*$Acolcurr*$HCol_trib + (-1./$nCols)*$TotWeightCapCurr]\
					0.\
					0.\
					0.;";
		}
	} else {
		load $nodeID 0.\
					0.\
					[expr -0.5*$wconc*$Acolcurr*$HCol_trib]\
					0.\
					0.\
					0.;  
		
		if {$write_model_files == 1} {
		puts $ExportID "load $nodeID 0.\
					0.\
					[expr -0.5*$wconc*$Acolcurr*$HCol_trib]\
					0.\
					0.\
					0.;";
		}
	}
	
}
}

if {$write_model_files == 1} {
puts $ExportID "\}"
}
