# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

proc connectColBaseToDeck {iNodeTag jNodeTag nEleCol colEleType nIntPts locations weights secTags transfTagCol transfTagColDeckConnection maxIters tol ColElements hingeParams intType \
	existsBH BHEleType BHintPts BHSecTag BHintType BHElements deckNode Ubig bond BondSecTag skewr ExportID ColEleSecFileID} {
    
    # proc currently assumes secTags is list of length equals one of the following
	# 1. 1 : same section throughout the length of the column
	# 2. 2 : secTags[0] for col bottom through all but top, and secTags[1] for top
	# 3. 3 : secTags[0] for col bottom, secTags[1] for all but bottom and top, secTags[2] for col top
	
	# Modify secTags to accomodate the above
	# Case 1.
	if {([llength $secTags] != 3) && ([llength $secTags] == 1)} {
		set secTags [list $secTags $secTags $secTags];
	}
	# Case 2.
	if {([llength $secTags] != 3) && ([llength $secTags] == 2)} {
		set secTags [linsert $secTags 0 [lindex $secTags 0]];
	}
    
	if {$existsBH == 1} {
        set BHSecTagThis [expr $BHSecTag + ([lindex $secTags 0] - ([lindex $secTags 0]/1000)*1000)]
		set BHElements [createBaseHingeElem $iNodeTag $BHEleType $BHintPts $BHSecTagThis $transfTagCol $BHintType $BHElements $bond $BondSecTag $skewr $ExportID]
	}
    
    if {[llength $hingeParams] == 0} {
        set ColElements [createColElem $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $maxIters $tol $ColElements $ExportID $ColEleSecFileID];
	} elseif {[llength $hingeParams] == 4} {
        set integration "$intType $nIntPts [lindex $secTags 0] [lindex $hingeParams 0] [lindex $hingeParams 2] [lindex $secTags 1] [lindex $hingeParams 1] [lindex $hingeParams 3] [lindex $secTags 2]"
        element $colEleType $iNodeTag $iNodeTag $jNodeTag $transfTagCol $integration -iter $maxIters $tol;
        if {$ExportID != -1} {
            puts $ExportID "element $colEleType $iNodeTag $iNodeTag $jNodeTag $transfTagCol $integration -iter $maxIters $tol;";
        }
        lappend ColElements $iNodeTag;
        
        # write a file containing all eleTags and corresponding secTags for the current column
        if {$ColEleSecFileID != -1} {
            set putsVal "$iNodeTag";
            append putsVal " " "[lindex $secTags 0] [lindex $secTags 1] [lindex $secTags 1] [lindex $secTags 2]";
            puts $ColEleSecFileID $putsVal;
        }
    } elseif {[llength $hingeParams] == 2} {
        set integration "$intType [lindex $secTags 0] [lindex $hingeParams 0] [lindex $secTags 2] [lindex $hingeParams 1] [lindex $secTags 1]"
        element $colEleType $iNodeTag $iNodeTag $jNodeTag $transfTagCol $integration -iter $maxIters $tol;
        if {$ExportID != -1} {
            puts $ExportID "element $colEleType $iNodeTag $iNodeTag $jNodeTag $transfTagCol $integration -iter $maxIters $tol;";
        }
        lappend ColElements $iNodeTag;
        
        # write a file containing all eleTags and corresponding secTags for the current column
        if {$ColEleSecFileID != -1} {
            set putsVal "$iNodeTag";
            if {$intType == "HingeRadau"} {
                append putsVal " " "[lindex $secTags 0] [lindex $secTags 1] [lindex $secTags 1] [lindex $secTags 1] [lindex $secTags 1] [lindex $secTags 2]";
            }
            if {$intType == "HingeRadauTwo"} {
                append putsVal " " "[lindex $secTags 0] [lindex $secTags 0] [lindex $secTags 1] [lindex $secTags 1] [lindex $secTags 2] [lindex $secTags 2]";
            }
            if {$intType == "HingeMidpoint" || $intType == "HingeEndpoint"} {
                append putsVal " " "[lindex $secTags 0] [lindex $secTags 1] [lindex $secTags 1] [lindex $secTags 2]";
            }
            puts $ColEleSecFileID $putsVal;
        }
    } 
    
    
    if {$bond == 1} {
        set tempTag [lindex $secTags [expr [llength $secTags] - 1]]
        set BondSecTagThis [expr $BondSecTag + ($tempTag - ($tempTag/1000)*1000)]
        element zeroLengthSection $jNodeTag $jNodeTag [expr $jNodeTag + 1] $BondSecTagThis -orient 0 0 1 [expr -1.*cos($skewr)] [expr -1.*sin($skewr)] 0. -doRayleigh 0; 
        element elasticBeamColumn [expr $jNodeTag + 1] [expr $jNodeTag + 1] $deckNode 1. $Ubig $Ubig 1. 1. 1. $transfTagColDeckConnection;
        if {$ExportID != -1} {
            puts $ExportID "element zeroLengthSection $jNodeTag $jNodeTag [expr $jNodeTag + 1] $BondSecTagThis -orient 0 0 1 [expr -1.*cos($skewr)] [expr -1.*sin($skewr)] 0. -doRayleigh 0;"
            puts $ExportID "element elasticBeamColumn [expr $jNodeTag + 1] [expr $jNodeTag + 1] $deckNode 1. $Ubig $Ubig 1. 1. 1. $transfTagColDeckConnection;"
        }
    } else {    
        element elasticBeamColumn $jNodeTag $jNodeTag $deckNode 1. $Ubig $Ubig 1. 1. 1. $transfTagColDeckConnection;  
        if {$ExportID != -1} {
            puts $ExportID "element elasticBeamColumn $jNodeTag $jNodeTag $deckNode 1. $Ubig $Ubig 1. 1. 1. $transfTagColDeckConnection;"
        }
	}
	return [list $BHElements $ColElements];
}
