# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------


proc connectSSI {iNodeTag jNodeTag ID_SSI dirList skewr ExportID} {
    set orientVec [list [expr -1.*sin($skewr)] [expr cos($skewr)] 0. [expr -1.*cos($skewr)] [expr -1.*sin($skewr)] 0.]
    set matTagList_stiff {};
    set matTagList_damp {};
    for {set i 0} {$i < [llength $dirList]} {incr i 1} {
        lappend matTagList_stiff [expr $ID_SSI + [lindex $dirList $i]]; 
        lappend matTagList_damp [expr ($ID_SSI * 10) + [lindex $dirList $i]]; 
    }
    element zeroLength [expr $iNodeTag + 1] $iNodeTag $jNodeTag -mat {*}$matTagList_stiff -dir {*}$dirList -orient {*}$orientVec;
    element zeroLength [expr $iNodeTag + 2] $iNodeTag $jNodeTag -mat {*}$matTagList_damp -dir {*}$dirList -orient {*}$orientVec;
    if {$ExportID != -1} {
        puts $ExportID "element zeroLength [expr $iNodeTag + 1] $iNodeTag $jNodeTag -mat $matTagList_stiff -dir $dirList -orient $orientVec;"
        puts $ExportID "element zeroLength [expr $iNodeTag + 2] $iNodeTag $jNodeTag -mat $matTagList_damp -dir $dirList -orient $orientVec;"
    }
} 


# ------------------------------------------------------------------------------------------------------------
# Piers
# ------------------------------------------------------------------------------------------------------------

if {$nCols > 1} {
    set ctrCol 0;
    for {set i 0} {$i < [llength $fixNodesColumn]} {incr i 1} {
        set bentDone 0;
        set currNodeTag [lindex $fixNodesColumn [expr $i]];
        incr ctrCol 1
        set currNodeCoords [nodeCoord $currNodeTag];
        
        if {$ctrCol == $nCols} {
            set bentDone 1;
            set ctrCol 0;
        }
        
        # baseHinge yes, bond yes, ssi no
        if {$existsBH == 1 && $bond == 1 && $ssi == 0} {
            fix [expr $currNodeTag * 10] 1 1 1 1 1 1;
            fix $currNodeTag 1 1 0 0 0 1;
        }
        
        # baseHinge yes, bond no, ssi no
        if {$existsBH == 1 && $bond == 0 && $ssi == 0} {
            fix $currNodeTag 1 1 1 1 1 1;
        }
        
        # baseHinge no, bond yes, ssi no
        if {$existsBH == 0 && $bond == 1 && $ssi == 0} {
            fix $currNodeTag 1 1 1 0 0 0;
        }
        
        # baseHinge no, bond no, ssi no
        if {$existsBH == 0 && $bond == 0 && $ssi == 0} {
            fix $currNodeTag 1 1 1 0 0 0;
        }
        
        # baseHinge yes, bond yes, ssi yes
        if {$existsBH == 1 && $bond == 1 && $ssi == 1} {
            # pass
        }
        
        # baseHinge yes, bond no, ssi yes
        if {$existsBH == 1 && $bond == 0 && $ssi == 1} {
            # pass
        }
        
        # baseHinge no, bond yes, ssi yes
        if {$existsBH == 0 && $bond == 1 && $ssi == 1} {
            node [expr ($currNodeTag - 1)] {*}$currNodeCoords;
            lappend allNodes [expr ($currNodeTag - 1)];
            if {$write_model_files == 1} {
                puts $ExportID "node [expr ($currNodeTag - 1)] $currNodeCoords;"
            }
            element zeroLength [expr ($currNodeTag - 1)] [expr ($currNodeTag - 1)] $currNodeTag -mat $ID_Rigid $ID_Rigid $ID_Rigid -dir 1 2 3 -orient [expr -1.*sin($skewr)] [expr cos($skewr)] 0. [expr -1.*cos($skewr)] [expr -1.*sin($skewr)] 0.;
            if {$write_model_files == 1} {
                puts $ExportID "element zeroLength [expr ($currNodeTag - 1)] [expr ($currNodeTag - 1)] $currNodeTag -mat $ID_Rigid $ID_Rigid $ID_Rigid -dir 1 2 3 -orient [expr -1.*sin($skewr)] [expr cos($skewr)] 0. [expr -1.*cos($skewr)] [expr -1.*sin($skewr)] 0.;"
            }
        }
        
        # baseHinge no, bond no, ssi yes
        if {$existsBH == 0 && $bond == 0 && $ssi == 1} {
            node [expr ($currNodeTag - 1)] {*}$currNodeCoords;
            lappend allNodes [expr ($currNodeTag - 1)];
            if {$write_model_files == 1} {
                puts $ExportID "node [expr ($currNodeTag - 1)] $currNodeCoords;"
            }
            element zeroLength [expr ($currNodeTag - 1)] [expr ($currNodeTag - 1)] $currNodeTag -mat $ID_Rigid $ID_Rigid $ID_Rigid -dir 1 2 3 -orient [expr -1.*sin($skewr)] [expr cos($skewr)] 0. [expr -1.*cos($skewr)] [expr -1.*sin($skewr)] 0.;
            if {$write_model_files == 1} {
                puts $ExportID "element zeroLength [expr ($currNodeTag - 1)] [expr ($currNodeTag - 1)] $currNodeTag -mat $ID_Rigid $ID_Rigid $ID_Rigid -dir 1 2 3 -orient [expr -1.*sin($skewr)] [expr cos($skewr)] 0. [expr -1.*cos($skewr)] [expr -1.*sin($skewr)] 0.;"
            }
        }
        
        if {$bentDone == 1} {
            # get all nodes in current bent
            set currBentNodeTags {}
            for {set j 0} {$j < $nCols} {incr j 1} {
                lappend currBentNodeTags [lindex $fixNodesColumn [expr $i - $j]];
            }
            
            # get coords of pivotNode
            set pivotNodeCoord_x 0.
            set pivotNodeCoord_y 0.
            set pivotNodeCoord_z 0.
            for {set j 0} {$j < [llength $currBentNodeTags]} {incr j 1} {
                set pivotNodeCoord_x [expr $pivotNodeCoord_x + [lindex [nodeCoord [lindex $currBentNodeTags $j]] 0]]
                set pivotNodeCoord_y [expr $pivotNodeCoord_y + [lindex [nodeCoord [lindex $currBentNodeTags $j]] 1]]
                set pivotNodeCoord_z [expr $pivotNodeCoord_z + [lindex [nodeCoord [lindex $currBentNodeTags $j]] 2]]
            }
            set pivotNodeCoord_x [expr 1./[llength $currBentNodeTags] * $pivotNodeCoord_x]
            set pivotNodeCoord_y [expr 1./[llength $currBentNodeTags] * $pivotNodeCoord_y]
            set pivotNodeCoord_z [expr 1./[llength $currBentNodeTags] * $pivotNodeCoord_z]
            set pivotNodeCoords [list $pivotNodeCoord_x $pivotNodeCoord_y $pivotNodeCoord_z];
            
            # baseHinge yes, bond yes, ssi yes
            if {$existsBH == 1 && $bond == 1 && $ssi == 1} {
                set pivotNodeTag [expr ($currNodeTag/1000) * 10000];
                if {[expr $nCols%2] == 0} {
                    node $pivotNodeTag {*}$pivotNodeCoords;
                    lappend allNodes $pivotNodeTag;
                    if {$write_model_files == 1} {
                        puts $ExportID "node $pivotNodeTag $pivotNodeCoords;"
                    }
                } 
                for {set j 0} {$j < [llength $currBentNodeTags]} {incr j 1} {
                    if {[expr [lindex $currBentNodeTags $j] * 10] != $pivotNodeTag} {
                        element elasticBeamColumn [expr $pivotNodeTag + $j + 1] $pivotNodeTag [expr [lindex $currBentNodeTags $j] * 10] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;
                        if {$write_model_files == 1} {
                            puts $ExportID "element elasticBeamColumn [expr $pivotNodeTag + $j + 1] $pivotNodeTag [expr [lindex $currBentNodeTags $j] * 10] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;"
                        }
                    }
                }
                node [expr $pivotNodeTag * 10] {*}$pivotNodeCoords;
                lappend allNodes [expr $pivotNodeTag * 10];
                if {$write_model_files == 1} {
                    puts $ExportID "node [expr $pivotNodeTag * 10] $pivotNodeCoords;"
                }
                fix [expr $pivotNodeTag * 10] 1 1 1 1 1 1;
                connectSSI [expr $pivotNodeTag * 10] $pivotNodeTag $ID_SSI [list 1 2 3 4 5 6] $skewr $ExportID
            }
            
            # baseHinge yes, bond no, ssi yes
            if {$existsBH == 1 && $bond == 0 && $ssi == 1} {
                set pivotNodeTag [expr ($currNodeTag/1000) * 1000];
                if {[expr $nCols%2] == 0} {
                    node $pivotNodeTag {*}$pivotNodeCoords;
                    lappend allNodes $pivotNodeTag;
                    if {$write_model_files == 1} {
                        puts $ExportID "node $pivotNodeTag $pivotNodeCoords;"
                    }
                } 
                for {set j 0} {$j < [llength $currBentNodeTags]} {incr j 1} {
                    if {[lindex $currBentNodeTags $j] != $pivotNodeTag} {
                        element elasticBeamColumn [expr $pivotNodeTag*10 + $j + 1] $pivotNodeTag [lindex $currBentNodeTags $j] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;
                        if {$write_model_files == 1} {
                            puts $ExportID "element elasticBeamColumn [expr $pivotNodeTag*10 + $j + 1] $pivotNodeTag [lindex $currBentNodeTags $j] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;"
                        }
                    }
                }
                node [expr $pivotNodeTag * 100] {*}$pivotNodeCoords;
                lappend allNodes [expr $pivotNodeTag * 100];
                if {$write_model_files == 1} {
                    puts $ExportID "node [expr $pivotNodeTag * 100] $pivotNodeCoords;"
                }
                fix [expr $pivotNodeTag * 100] 1 1 1 1 1 1;
                connectSSI [expr $pivotNodeTag * 100] $pivotNodeTag $ID_SSI [list 1 2 3 4 5 6] $skewr $ExportID
            }
            
            # baseHinge no, bond yes, ssi yes
            if {$existsBH == 0 && $bond == 1 && $ssi == 1} {
                set pivotNodeTag [expr ($currNodeTag/1000) * 1000];
                if {[expr $nCols%2] == 0} {
                    node $pivotNodeTag {*}$pivotNodeCoords;
                    lappend allNodes $pivotNodeTag;
                    if {$write_model_files == 1} {
                        puts $ExportID "node $pivotNodeTag $pivotNodeCoords;"
                    }
                } 
                for {set j 0} {$j < [llength $currBentNodeTags]} {incr j 1} {
                    if {[expr [lindex $currBentNodeTags $j] - 1] != $pivotNodeTag} {
                        element elasticBeamColumn [expr $pivotNodeTag*10 + $j + 1] $pivotNodeTag [expr [lindex $currBentNodeTags $j] - 1] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;
                        if {$write_model_files == 1} {
                            puts $ExportID "element elasticBeamColumn [expr $pivotNodeTag*10 + $j + 1] $pivotNodeTag [expr [lindex $currBentNodeTags $j] - 1] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;"
                        }
                    }
                }
                node [expr $pivotNodeTag * 100] {*}$pivotNodeCoords;
                lappend allNodes [expr $pivotNodeTag * 100];
                if {$write_model_files == 1} {
                    puts $ExportID "node [expr $pivotNodeTag * 100] $pivotNodeCoords;"
                }
                fix [expr $pivotNodeTag * 100] 1 1 1 1 1 1;
                connectSSI [expr $pivotNodeTag * 100] $pivotNodeTag $ID_SSI [list 1 2 3 4 5 6] $skewr $ExportID
            }
            
            # baseHinge no, bond no, ssi yes
            if {$existsBH == 0 && $bond == 0 && $ssi == 1} {
                set pivotNodeTag [expr ($currNodeTag/1000) * 1000];
                if {[expr $nCols%2] == 0} {
                    node $pivotNodeTag {*}$pivotNodeCoords;
                    lappend allNodes $pivotNodeTag;
                    if {$write_model_files == 1} {
                        puts $ExportID "node $pivotNodeTag $pivotNodeCoords;"
                    }
                } 
                for {set j 0} {$j < [llength $currBentNodeTags]} {incr j 1} {
                    if {[expr [lindex $currBentNodeTags $j] - 1] != $pivotNodeTag} {
                        element elasticBeamColumn [expr $pivotNodeTag*10 + $j + 1] $pivotNodeTag [expr [lindex $currBentNodeTags $j] - 1] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;
                        if {$write_model_files == 1} {
                            puts $ExportID "element elasticBeamColumn [expr $pivotNodeTag*10 + $j + 1] $pivotNodeTag [expr [lindex $currBentNodeTags $j] - 1] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;"
                        }
                    }
                }
                node [expr $pivotNodeTag * 100] {*}$pivotNodeCoords;
                lappend allNodes [expr $pivotNodeTag * 100];
                if {$write_model_files == 1} {
                    puts $ExportID "node [expr $pivotNodeTag * 100] $pivotNodeCoords;"
                }
                fix [expr $pivotNodeTag * 100] 1 1 1 1 1 1;
                connectSSI [expr $pivotNodeTag * 100] $pivotNodeTag $ID_SSI [list 1 2 3 4 5 6] $skewr $ExportID
            }
        }
        
    }
}

# ------------------------------------------------------------------------------------------------------------
# Abutments
# ------------------------------------------------------------------------------------------------------------
for {set i 0} {$i <= [expr [llength $AbutmentFixedNodes] - 1]} {incr i 1} {
	fix [expr [lindex $AbutmentFixedNodes $i]] 1 1 1 1 1 1;
}