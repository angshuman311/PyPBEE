# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

set rotMassSwitchTorsion 0.;
set rotMassSwitch 0.;


###################################################################################################
if {$write_model_files == 1} {
puts $ExportID "#Deck Secondary Nodes"
}
###################################################################################################

for {set i 0} {$i <= [expr [llength $secondaryNodesDeck] - 1]} {incr i 1} {

	set nodeID [expr [lindex $secondaryNodesDeck $i]];
	set L_trib [expr ([lindex $L [expr $nodeID/100]])/$nEleSpan];
	
	mass $nodeID [expr $mconc*$Adeck*$L_trib]\
				[expr $mconc*$Adeck*$L_trib]\
				[expr $mconc*$Adeck*$L_trib]\
				[expr $rotMassSwitchTorsion*$mconc*$L_trib*($Iydeck + $Izdeck)] \
				[expr $rotMassSwitch*$mconc*$L_trib*($Iydeck)] \
				[expr $rotMassSwitch*$mconc*$L_trib*($Izdeck)];
	
	if {$write_model_files == 1} {
	puts $ExportID "mass $nodeID [expr $mconc*$Adeck*$L_trib]\
				[expr $mconc*$Adeck*$L_trib]\
				[expr $mconc*$Adeck*$L_trib]\
				[expr $rotMassSwitchTorsion*$mconc*$L_trib*($Iydeck + $Izdeck)] \
				[expr $rotMassSwitch*$mconc*$L_trib*($Iydeck)] \
				[expr $rotMassSwitch*$mconc*$L_trib*($Izdeck)];"
	}
	
}




###################################################################################################
if {$write_model_files == 1} {
puts $ExportID "#Deck Primary Nodes"
}
###################################################################################################

for {set i 0} {$i <= [expr [llength $primaryNodesDeck] - 1]} {incr i 1} {

	set nodeID [expr [lindex $primaryNodesDeck $i]];
	
	if {$i == 0} {
		set L_trib [expr ([lindex $L [expr $nodeID/100]])/$nEleSpan];
		mass $nodeID [expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr $rotMassSwitchTorsion*0.5*$mconc*$L_trib*($Iydeck + $Izdeck)] \
					[expr $rotMassSwitch*0.5*$mconc*$L_trib*($Iydeck)] \
					[expr $rotMassSwitch*0.5*$mconc*$L_trib*($Izdeck)];
		
		if {$write_model_files == 1} {
		puts $ExportID "mass $nodeID [expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr $rotMassSwitchTorsion*0.5*$mconc*$L_trib*($Iydeck + $Izdeck)] \
					[expr $rotMassSwitch*0.5*$mconc*$L_trib*($Iydeck)] \
					[expr $rotMassSwitch*0.5*$mconc*$L_trib*($Izdeck)];"
		}
		
	} elseif {$i == $nSpans} {
		set L_trib [expr ([lindex $L end])/$nEleSpan];
		mass $nodeID [expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr $rotMassSwitchTorsion*0.5*$mconc*$L_trib*($Iydeck + $Izdeck)] \
					[expr $rotMassSwitch*0.5*$mconc*$L_trib*($Iydeck)] \
					[expr $rotMassSwitch*0.5*$mconc*$L_trib*($Izdeck)];
		
		if {$write_model_files == 1} {
		puts $ExportID "mass $nodeID [expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr 0.5*$mconc*$Adeck*$L_trib]\
					[expr $rotMassSwitchTorsion*0.5*$mconc*$L_trib*($Iydeck + $Izdeck)] \
					[expr $rotMassSwitch*0.5*$mconc*$L_trib*($Iydeck)] \
					[expr $rotMassSwitch*0.5*$mconc*$L_trib*($Izdeck)];"
		}
		
	} else {
		set L_trib_ahead [expr ([lindex $L [expr $nodeID/100]])/$nEleSpan];
		set L_trib_behind [expr ([lindex $L [expr $nodeID/100 - 1]])/$nEleSpan];
		mass $nodeID [expr 0.5*$mconc*$Adeck*$L_trib_ahead + 0.5*$mconc*$Adeck*$L_trib_behind]\
					 [expr 0.5*$mconc*$Adeck*$L_trib_ahead + 0.5*$mconc*$Adeck*$L_trib_behind]\
					 [expr 0.5*$mconc*$Adeck*$L_trib_ahead + 0.5*$mconc*$Adeck*$L_trib_behind]\
					 [expr $rotMassSwitchTorsion*0.5*$mconc*$L_trib_ahead*($Iydeck + $Izdeck) + $rotMassSwitchTorsion*0.5*$mconc*$L_trib_behind*($Iydeck + $Izdeck)] \
					 [expr $rotMassSwitch*0.5*$mconc*$L_trib_ahead*($Iydeck)           + $rotMassSwitch*0.5*$mconc*$L_trib_behind*($Iydeck)] \
					 [expr $rotMassSwitch*0.5*$mconc*$L_trib_ahead*($Izdeck)           + $rotMassSwitch*0.5*$mconc*$L_trib_behind*($Izdeck)];
		
		if {$write_model_files == 1} {
		puts $ExportID "mass $nodeID [expr 0.5*$mconc*$Adeck*$L_trib_ahead + 0.5*$mconc*$Adeck*$L_trib_behind]\
					 [expr 0.5*$mconc*$Adeck*$L_trib_ahead + 0.5*$mconc*$Adeck*$L_trib_behind]\
					 [expr 0.5*$mconc*$Adeck*$L_trib_ahead + 0.5*$mconc*$Adeck*$L_trib_behind]\
					 [expr $rotMassSwitchTorsion*0.5*$mconc*$L_trib_ahead*($Iydeck + $Izdeck) + $rotMassSwitchTorsion*0.5*$mconc*$L_trib_behind*($Iydeck + $Izdeck)] \
					 [expr $rotMassSwitch*0.5*$mconc*$L_trib_ahead*($Iydeck)           + $rotMassSwitch*0.5*$mconc*$L_trib_behind*($Iydeck)] \
					 [expr $rotMassSwitch*0.5*$mconc*$L_trib_ahead*($Izdeck)           + $rotMassSwitch*0.5*$mconc*$L_trib_behind*($Izdeck)];";
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
		mass $nodeID [expr $mconc*$Acolcurr*$HCol_trib]\
					 [expr $mconc*$Acolcurr*$HCol_trib]\
					 [expr $mconc*$Acolcurr*$HCol_trib]\
					 [expr $rotMassSwitch*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]])] \
					 [expr $rotMassSwitch*$mconc*$HCol_trib*([lindex $Izcol [expr $nodeID/1000 - 1]])] \
					 [expr $rotMassSwitchTorsion*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]] + [lindex $Izcol [expr $nodeID/1000 - 1]])];  
		if {$write_model_files == 1} {
		puts $ExportID "mass $nodeID [expr $mconc*$Acolcurr*$HCol_trib]\
					 [expr $mconc*$Acolcurr*$HCol_trib]\
					 [expr $mconc*$Acolcurr*$HCol_trib]\
					 [expr $rotMassSwitch*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]])] \
					 [expr $rotMassSwitch*$mconc*$HCol_trib*([lindex $Izcol [expr $nodeID/1000 - 1]])] \
					 [expr $rotMassSwitchTorsion*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]] + [lindex $Izcol [expr $nodeID/1000 - 1]])];";
		}
		
	}

}



###################################################################################################
if {$write_model_files == 1} {
puts $ExportID "#Column Primary Nodes"
}
###################################################################################################

for {set i 0} {$i <= [expr [llength $primaryNodesColumn] - 1]} {incr i 1} {
	set nodeID [expr [lindex $primaryNodesColumn $i]];
	
	set Acapcurr [expr [lindex $Acap [expr $nodeID/1000 - 1]]];
	set LcapcurrR [expr [lindex $Lcap [expr 2*($nodeID/1000 - 1) - 1]]];
	set LcapcurrL [expr [lindex $Lcap [expr 2*($nodeID/1000 - 1)]]];
	set TotalLcapcurr [expr $LcapcurrL + $LcapcurrR];
	set TotalLcapcurrSkewed [expr $TotalLcapcurr/cos($skew*$pi/180.0)];
	set TotMassCapCurr [expr 1.*$mconc*$Acapcurr*$TotalLcapcurrSkewed];
	set Acolcurr [expr [lindex $Acol [expr $nodeID/1000 - 1]]];
	set Dcolcurr [expr [lindex $Dcol [expr $nodeID/1000 - 1]]];
	set HCol_trib [expr ([lindex $HCol [expr $nodeID/1000 - 1]])/$nEleCol];
	set tempnodeID [expr $nodeID - ($nodeID/100)*100]
	if {$tempnodeID == [expr $nEleCol + 1]} {
		mass $nodeID [expr 0.5*$mconc*$Acolcurr*$HCol_trib + (1./$nCols)*$TotMassCapCurr]\
					[expr 0.5*$mconc*$Acolcurr*$HCol_trib  + (1./$nCols)*$TotMassCapCurr]\
					[expr 0.5*$mconc*$Acolcurr*$HCol_trib  + (1./$nCols)*$TotMassCapCurr]\
					[expr $rotMassSwitch*0.5*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]])] \
					[expr $rotMassSwitch*0.5*$mconc*$HCol_trib*([lindex $Izcol [expr $nodeID/1000 - 1]])] \
					[expr $rotMassSwitchTorsion*0.5*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]] + [lindex $Izcol [expr $nodeID/1000 - 1]])];  
		
		if {$write_model_files == 1} {
		puts $ExportID "mass $nodeID [expr 0.5*$mconc*$Acolcurr*$HCol_trib + (1./$nCols)*$TotMassCapCurr]\
					[expr 0.5*$mconc*$Acolcurr*$HCol_trib  + (1./$nCols)*$TotMassCapCurr]\
					[expr 0.5*$mconc*$Acolcurr*$HCol_trib  + (1./$nCols)*$TotMassCapCurr]\
					[expr $rotMassSwitch*0.5*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]])] \
					[expr $rotMassSwitch*0.5*$mconc*$HCol_trib*([lindex $Izcol [expr $nodeID/1000 - 1]])] \
					[expr $rotMassSwitchTorsion*0.5*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]] + [lindex $Izcol [expr $nodeID/1000 - 1]])];";
		}
	} else {
	mass $nodeID [expr 0.5*$mconc*$Acolcurr*$HCol_trib] \
				 [expr 0.5*$mconc*$Acolcurr*$HCol_trib] \
				 [expr 0.5*$mconc*$Acolcurr*$HCol_trib] \
				 [expr $rotMassSwitch*0.5*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]])] \
				 [expr $rotMassSwitch*0.5*$mconc*$HCol_trib*([lindex $Izcol [expr $nodeID/1000 - 1]])] \
				 [expr $rotMassSwitchTorsion*0.5*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]] + [lindex $Izcol [expr $nodeID/1000 - 1]])];
	
	if {$write_model_files == 1} {
	puts $ExportID "mass $nodeID [expr 0.5*$mconc*$Acolcurr*$HCol_trib] \
				 [expr 0.5*$mconc*$Acolcurr*$HCol_trib] \
				 [expr 0.5*$mconc*$Acolcurr*$HCol_trib] \
				 [expr $rotMassSwitch*0.5*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]])] \
				 [expr $rotMassSwitch*0.5*$mconc*$HCol_trib*([lindex $Izcol [expr $nodeID/1000 - 1]])] \
				 [expr $rotMassSwitchTorsion*0.5*$mconc*$HCol_trib*([lindex $Iycol [expr $nodeID/1000 - 1]] + [lindex $Izcol [expr $nodeID/1000 - 1]])];"
	}
	}
	
	
}
