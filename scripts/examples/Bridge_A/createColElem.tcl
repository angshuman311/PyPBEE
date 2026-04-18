# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

# createColElem $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $maxIters $tol $ColElements $ExportID $ColEleSecFileID
proc createColElem {iNodeTag jNodeTag nEleCol colEleType nIntPts locations weights secTags transfTagCol maxIters tol ColElements ExportID ColEleSecFileID} {
	
	# Create col elements from bottom to top
	# Loop over one less than number of elements 
	for {set i1 1} {$i1 <= [expr ($nEleCol - 1)]} {incr i1 1} {
		# in here => there are more than one elements. [lindex $secTags 2] will be taken care of outside loop
		# create list to hold secTags for each integration point in the element
		set colElemSecTags {};
		if {$i1 == 1} {
			# First element, use [lindex $secTags 0] for first integration point
			# and [lindex $secTags 1] for rest
			lappend colElemSecTags [expr [lindex $secTags 0]];
			for {set intPtIter 1} {$intPtIter < $nIntPts} {incr intPtIter 1} {
				lappend colElemSecTags [expr [lindex $secTags 1]];
			}
		} else {
			# use [lindex $secTags 1] for all
			for {set intPtIter 1} {$intPtIter <= $nIntPts} {incr intPtIter 1} {
				lappend colElemSecTags [expr [lindex $secTags 1]];
			}
		}
		set integration "UserDefined $nIntPts $colElemSecTags $locations $weights";
		
		# create element
		element $colEleType [expr $iNodeTag + ($i1 - 1)] [expr $iNodeTag + ($i1 - 1)] [expr $iNodeTag + $i1] $transfTagCol $integration -iter $maxIters $tol;
		if {$ExportID != -1} {
		puts $ExportID "element $colEleType [expr $iNodeTag + ($i1 - 1)] [expr $iNodeTag + ($i1 - 1)] [expr $iNodeTag + $i1] $transfTagCol $integration -iter $maxIters $tol;";
		}
		lappend ColElements [expr $iNodeTag + ($i1 - 1)];
		
		# write a file containing all eleTags and corresponding secTags for the current column
		if {$ColEleSecFileID != -1} {
			set putsVal "[expr $iNodeTag + ($i1 - 1)]";
			append putsVal " " "$colElemSecTags";
			puts $ColEleSecFileID $putsVal;
		}
	}
	
	# here => last element or the only element
	# create list to hold secTags for each integration point in the element
	set colElemSecTags {};
	if {$nEleCol == 1} {
		# Only one element, use [lindex $secTags 0] for first integration point
		# [lindex $secTags 1] for all but first and last, and [lindex $secTags 2] for last
		lappend colElemSecTags [expr [lindex $secTags 0]];
		for {set intPtIter 1} {$intPtIter < [expr $nIntPts - 1]} {incr intPtIter 1} {
			lappend colElemSecTags [expr [lindex $secTags 1]];
		}
		lappend colElemSecTags [expr [lindex $secTags 2]];
	} else {
		# last element, use [lindex $secTags 1] for all but last, and [lindex $secTags 2] for last
		for {set intPtIter 1} {$intPtIter < $nIntPts} {incr intPtIter 1} {
			lappend colElemSecTags [expr [lindex $secTags 1]];
		}
		lappend colElemSecTags [expr [lindex $secTags 2]];
	}
	set integration "UserDefined $nIntPts $colElemSecTags $locations $weights";
	
	# create element
	element $colEleType [expr $jNodeTag - 1] [expr $jNodeTag - 1] $jNodeTag $transfTagCol $integration -iter $maxIters $tol;
	if {$ExportID != -1} {
	puts $ExportID "element $colEleType [expr $jNodeTag - 1] [expr $jNodeTag - 1] $jNodeTag $transfTagCol $integration -iter $maxIters $tol;";
	}
	lappend ColElements [expr $jNodeTag - 1];
	
	# write a file containing all eleTags and corresponding secTags for the current column
	if {$ColEleSecFileID != -1} {
		set putsVal "[expr $jNodeTag - 1]";
		append putsVal " " "$colElemSecTags";
		puts $ColEleSecFileID $putsVal;
	}
	
	return $ColElements;
}
