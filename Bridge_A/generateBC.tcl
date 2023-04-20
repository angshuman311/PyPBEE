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

if {$nCols == 1} {
    for {set i 0} {$i < [llength $fixNodesColumn]} {incr i 1} {
		
        set currNodeTag [lindex $fixNodesColumn [expr $i]];
        set currNodeCoords [nodeCoord $currNodeTag];
        
        # bond yes, ssi no
        if {$bond == 1 && $ssi == 0} {
            fix [expr ($currNodeTag - 1) * 10] 1 1 1 1 1 1;
            fix $currNodeTag 1 1 0 0 0 1;
			if {$write_model_files == 1} {
				puts $ExportID "fix [expr ($currNodeTag - 1) * 10] 1 1 1 1 1 1;"
				puts $ExportID "fix $currNodeTag 1 1 0 0 0 1;"
			}
        }
        
        # bond no, ssi no
        if {$bond == 0 && $ssi == 0} {
            fix $currNodeTag 1 1 1 1 1 1;
			if {$write_model_files == 1} {
				puts $ExportID "fix $currNodeTag 1 1 1 1 1 1;"
			}
        }
        
        # bond yes, ssi yes
        if {$bond == 1 && $ssi == 1} {
            node [expr ($currNodeTag - 1) * 100] {*}$currNodeCoords;
            lappend allNodes [expr ($currNodeTag - 1) * 100]
            if {$write_model_files == 1} {
                puts $ExportID "node [expr ($currNodeTag - 1) * 100] $currNodeCoords;"
            }
            fix [expr ($currNodeTag - 1) * 100] 1 1 1 1 1 1;
			if {$write_model_files == 1} {
				puts $ExportID "fix [expr ($currNodeTag - 1) * 100] 1 1 1 1 1 1;"
			}
            connectSSI [expr ($currNodeTag - 1) * 100] [expr ($currNodeTag - 1) * 10] $ID_SSI [list 1 2 3 4 5 6] $skewr $ExportID
        }
        
        # bond no, ssi yes
        if {$bond == 0 && $ssi == 1} {
            node [expr ($currNodeTag - 1) * 100] {*}$currNodeCoords;
            lappend allNodes [expr ($currNodeTag - 1) * 100];
            if {$write_model_files == 1} {
                puts $ExportID "node [expr ($currNodeTag - 1) * 100] $currNodeCoords;"
            }
            fix [expr ($currNodeTag - 1) * 100] 1 1 1 1 1 1;
			if {$write_model_files == 1} {
				puts $ExportID "fix [expr ($currNodeTag - 1) * 100] 1 1 1 1 1 1;"
			}
            connectSSI [expr ($currNodeTag - 1) * 100] $currNodeTag $ID_SSI [list 1 2 3 4 5 6] $skewr $ExportID   
        }
    }

}

if {$nCols > 1} {
    for {set i 0} {$i < [llength $fixNodesColumn]} {incr i 1} {
        
        set currNodeTag [lindex $fixNodesColumn [expr $i]];
        set currNodeCoords [nodeCoord $currNodeTag];
        
        # baseHinge yes, bond yes, ssi no
        if {$existsBH == 1 && $bond == 1 && $ssi == 0} {
            fix [expr $currNodeTag * 10] 1 1 1 1 1 1;
            fix $currNodeTag 1 1 0 0 0 1;
			if {$write_model_files == 1} {
				puts $ExportID "fix [expr $currNodeTag * 10] 1 1 1 1 1 1;"
				puts $ExportID "fix $currNodeTag 1 1 0 0 0 1;"
			}
        }
        
        # baseHinge yes, bond no, ssi no
        if {$existsBH == 1 && $bond == 0 && $ssi == 0} {
            fix $currNodeTag 1 1 1 1 1 1;
			if {$write_model_files == 1} {
				puts $ExportID "fix $currNodeTag 1 1 1 1 1 1;"
			}
        }
        
        # baseHinge no, bond yes, ssi no
        if {$existsBH == 0 && $bond == 1 && $ssi == 0} {
            fix $currNodeTag 1 1 1 0 0 0;
			if {$write_model_files == 1} {
				puts $ExportID "fix $currNodeTag 1 1 1 0 0 0;"
			}
        }
        
        # baseHinge no, bond no, ssi no
        if {$existsBH == 0 && $bond == 0 && $ssi == 0} {
            fix $currNodeTag 1 1 1 0 0 0;
			if {$write_model_files == 1} {
				puts $ExportID "fix $currNodeTag 1 1 1 0 0 0;"
			}
        }
        
        # baseHinge yes, bond yes, ssi yes
        if {$existsBH == 1 && $bond == 1 && $ssi == 1} {
            node [expr $currNodeTag * 100] {*}$currNodeCoords;
            lappend allNodes [expr $currNodeTag * 100];
            if {$write_model_files == 1} {
                puts $ExportID "node [expr $currNodeTag * 100] $currNodeCoords;"
            }
            fix [expr $currNodeTag * 100] 1 1 1 1 1 1;
			if {$write_model_files == 1} {
				puts $ExportID "fix [expr $currNodeTag * 100] 1 1 1 1 1 1;"
			}
            connectSSI [expr $currNodeTag * 100] [expr $currNodeTag * 10] $ID_SSI [list 1 2 3 4 5 6] $skewr $ExportID
        }
        
        # baseHinge yes, bond no, ssi yes
        if {$existsBH == 1 && $bond == 0 && $ssi == 1} {
            node [expr $currNodeTag * 100] {*}$currNodeCoords;
            lappend allNodes [expr $currNodeTag * 100];
            if {$write_model_files == 1} {
                puts $ExportID "node [expr $currNodeTag * 100] $currNodeCoords;"
            }
            fix [expr $currNodeTag * 100] 1 1 1 1 1 1;
			if {$write_model_files == 1} {
				puts $ExportID "fix [expr $currNodeTag * 100] 1 1 1 1 1 1;"
			}
            connectSSI [expr $currNodeTag * 100] $currNodeTag $ID_SSI [list 1 2 3 4 5 6] $skewr $ExportID
        }
        
        # baseHinge no, bond yes, ssi yes
        if {$existsBH == 0 && $bond == 1 && $ssi == 1} {
            node [expr ($currNodeTag - 1) * 100] {*}$currNodeCoords;
            lappend allNodes [expr ($currNodeTag - 1) * 100];
            if {$write_model_files == 1} {
                puts $ExportID "node [expr ($currNodeTag - 1) * 100] $currNodeCoords;"
            }
            fix [expr ($currNodeTag - 1) * 100] 1 1 1 1 1 1;
			if {$write_model_files == 1} {
				puts $ExportID "fix [expr ($currNodeTag - 1) * 100] 1 1 1 1 1 1;"
			}
            connectSSI [expr ($currNodeTag - 1) * 100] $currNodeTag $ID_SSI [list 1 2 3] $skewr $ExportID
        }
        
        # baseHinge no, bond no, ssi yes
        if {$existsBH == 0 && $bond == 0 && $ssi == 1} {
            node [expr ($currNodeTag - 1) * 100] {*}$currNodeCoords;
            lappend allNodes [expr ($currNodeTag - 1) * 100];
            if {$write_model_files == 1} {
                puts $ExportID "node [expr ($currNodeTag - 1) * 100] $currNodeCoords;"
            }
            fix [expr ($currNodeTag - 1) * 100] 1 1 1 1 1 1;
			if {$write_model_files == 1} {
				puts $ExportID "fix [expr ($currNodeTag - 1) * 100] 1 1 1 1 1 1;"
			}
            connectSSI [expr ($currNodeTag - 1) * 100] $currNodeTag $ID_SSI [list 1 2 3] $skewr $ExportID
        }
        
    }
}

# ------------------------------------------------------------------------------------------------------------
# Abutments
# ------------------------------------------------------------------------------------------------------------
for {set i 0} {$i < [llength $AbutmentFixedNodes]} {incr i 1} {
	fix [expr [lindex $AbutmentFixedNodes $i]] 1 1 1 1 1 1;
	if {$write_model_files == 1} {
		puts $ExportID "fix [expr [lindex $AbutmentFixedNodes $i]] 1 1 1 1 1 1;"
	}
}