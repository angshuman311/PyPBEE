# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

#Global X - along longitudinal direction of bridge 
#Global Y - along transverse direction of bridge
#Global Z - along height of bridge
#Local X - along height of bridge - along Global Z 
#Local Y - normal to skewed bent
#Local Z - along skewed bent 

set primaryNodesColumn {};
set secondaryNodesColumn {};
set fixNodesColumn {};

set Xi [expr [lindex $L 1]];		# abscissa of primary nodes  
for {set i 1} {$i <= $nBent} {incr i 1} {
	set colHgtCurrentBent [expr [lindex $HCol $i]];
	set BaseHingeCurrentHeight [expr [lindex $BaseHingeHeight $i]];
	if {[expr $nCols%2] == 0} {
		for {set k 1} {$k <= [expr $nCols/2]} {incr k 1} {
		
			set yOffset [expr [lindex $Lcap [expr 2*$i - 1]] - ($nCols/2 - $k)*([lindex $Lcap [expr 2*$i - 1]])/($nCols/2)];
			set xOffset [expr [lindex $dxc [expr 2*$i - 1]] - ($nCols/2 - $k)*([lindex $dxc [expr 2*$i - 1]])/($nCols/2)]
			node [expr ($i + 1)*1000 + (2*$k - 1)*100] [expr $Xi + $xOffset] [expr -$yOffset] [expr ($HColumnMax - $colHgtCurrentBent) - $BaseHingeCurrentHeight];
			if {$write_model_files == 1} {
			puts $ExportID "node [expr ($i + 1)*1000 + (2*$k - 1)*100] [expr $Xi + $xOffset] [expr -$yOffset] [expr ($HColumnMax - $colHgtCurrentBent) - $BaseHingeCurrentHeight];"
			}
			
			set yOffset [expr [lindex $Lcap [expr 2*$i]] - ($nCols/2 - $k)*([lindex $Lcap [expr 2*$i]])/($nCols/2)];
			set xOffset [expr [lindex $dxc [expr 2*$i]] - ($nCols/2 - $k)*([lindex $dxc [expr 2*$i]])/($nCols/2)]
			node [expr ($i + 1)*1000 + (2*$k)*100] [expr $Xi - $xOffset] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) - $BaseHingeCurrentHeight];
			if {$write_model_files == 1} {
			puts $ExportID "node [expr ($i + 1)*1000 + (2*$k)*100] [expr $Xi - $xOffset] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) - $BaseHingeCurrentHeight];"
			}
			
			lappend fixNodesColumn [expr ($i + 1)*1000 + (2*$k - 1)*100];
			lappend fixNodesColumn [expr ($i + 1)*1000 + (2*$k)*100];
			
			for {set j 1} {$j <= [expr $nEleCol + 1]} {incr j 1} {
						
				set yOffset [expr [lindex $Lcap [expr 2*$i - 1]] - ($nCols/2 - $k)*([lindex $Lcap [expr 2*$i - 1]])/($nCols/2)];
				set xOffset [expr [lindex $dxc [expr 2*$i - 1]] - ($nCols/2 - $k)*([lindex $dxc [expr 2*$i - 1]])/($nCols/2)]
				node [expr ($i + 1)*1000 + (2*$k - 1)*100 + $j] [expr $Xi + $xOffset] [expr -$yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];
				if {$write_model_files == 1} {
				puts $ExportID "node [expr ($i + 1)*1000 + (2*$k - 1)*100 + $j] [expr $Xi + $xOffset] [expr -$yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];"
				}
				
				set yOffset [expr [lindex $Lcap [expr 2*$i]] - ($nCols/2 - $k)*([lindex $Lcap [expr 2*$i]])/($nCols/2)];
				set xOffset [expr [lindex $dxc [expr 2*$i]] - ($nCols/2 - $k)*([lindex $dxc [expr 2*$i]])/($nCols/2)]
				node [expr ($i + 1)*1000 + (2*$k)*100 + $j] [expr $Xi - $xOffset] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];
				if {$write_model_files == 1} {
				puts $ExportID "node [expr ($i + 1)*1000 + (2*$k)*100 + $j] [expr $Xi - $xOffset] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];"
				}
				
				if {$j == 1 || $j == [expr $nEleCol + 1]} {
					lappend primaryNodesColumn [expr ($i + 1)*1000 + (2*$k - 1)*100 + $j];
					lappend primaryNodesColumn [expr ($i + 1)*1000 + (2*$k)*100 + $j];
				} else {
					lappend secondaryNodesColumn [expr ($i + 1)*1000 + (2*$k - 1)*100 + $j];
					lappend secondaryNodesColumn [expr ($i + 1)*1000 + (2*$k)*100 + $j];
				}	
			}		
		}
	} else {
		if {$nCols != 1} {
			set yOffset 0.0;
			node [expr ($i + 1)*1000] [expr $Xi] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) - $BaseHingeCurrentHeight];
			if {$write_model_files == 1} {
			puts $ExportID "node [expr ($i + 1)*1000] [expr $Xi] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) - $BaseHingeCurrentHeight];"
			}
			
			lappend fixNodesColumn [expr ($i + 1)*1000];
			
			for {set j 1} {$j <= [expr $nEleCol + 1]} {incr j 1} {
			
				set yOffset 0.0;
				node [expr ($i + 1)*1000 + $j] [expr $Xi] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];
				if {$write_model_files == 1} {
				puts $ExportID "node [expr ($i + 1)*1000 + $j] [expr $Xi] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];"
				}
				
				if {$j == 1 || $j == [expr $nEleCol + 1]} {
					lappend primaryNodesColumn [expr ($i + 1)*1000 + $j];
				} else {
					lappend secondaryNodesColumn [expr ($i + 1)*1000 + $j];
				}
			}
			
			for {set k 1} {$k <= [expr ($nCols-1)/2]} {incr k 1} {
		
				set yOffset [expr [lindex $Lcap [expr 2*$i - 1]] - (($nCols-1)/2 - $k)*([lindex $Lcap [expr 2*$i - 1]])/(($nCols-1)/2)];
				set xOffset [expr [lindex $dxc [expr 2*$i - 1]] - (($nCols-1)/2 - $k)*([lindex $dxc [expr 2*$i - 1]])/(($nCols-1)/2)]
				node [expr ($i + 1)*1000 + (2*$k - 1)*100] [expr $Xi + $xOffset] [expr -$yOffset] [expr ($HColumnMax - $colHgtCurrentBent) - $BaseHingeCurrentHeight];
				if {$write_model_files == 1} {
				puts $ExportID "node [expr ($i + 1)*1000 + (2*$k - 1)*100] [expr $Xi + $xOffset] [expr -$yOffset] [expr ($HColumnMax - $colHgtCurrentBent) - $BaseHingeCurrentHeight];"
				}
				
				set yOffset [expr [lindex $Lcap [expr 2*$i]] - (($nCols-1)/2 - $k)*([lindex $Lcap [expr 2*$i]])/(($nCols-1)/2)];
				set xOffset [expr [lindex $dxc [expr 2*$i]] - (($nCols-1)/2 - $k)*([lindex $dxc [expr 2*$i]])/(($nCols-1)/2)]
				node [expr ($i + 1)*1000 + (2*$k)*100] [expr $Xi - $xOffset] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) - $BaseHingeCurrentHeight];
				if {$write_model_files == 1} {
				puts $ExportID "node [expr ($i + 1)*1000 + (2*$k)*100] [expr $Xi - $xOffset] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) - $BaseHingeCurrentHeight];"
				}
				
				lappend fixNodesColumn [expr ($i + 1)*1000 + (2*$k - 1)*100];
				lappend fixNodesColumn [expr ($i + 1)*1000 + (2*$k)*100];
				
				for {set j 1} {$j <= [expr $nEleCol + 1]} {incr j 1} {
				
					set yOffset [expr [lindex $Lcap [expr 2*$i - 1]] - (($nCols-1)/2 - $k)*([lindex $Lcap [expr 2*$i - 1]])/(($nCols-1)/2)];
					set xOffset [expr [lindex $dxc [expr 2*$i - 1]] - (($nCols-1)/2 - $k)*([lindex $dxc [expr 2*$i - 1]])/(($nCols-1)/2)]
					node [expr ($i + 1)*1000 + (2*$k - 1)*100 + $j] [expr $Xi + $xOffset] [expr -$yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];
					if {$write_model_files == 1} {
					puts $ExportID "node [expr ($i + 1)*1000 + (2*$k - 1)*100 + $j] [expr $Xi + $xOffset] [expr -$yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];"
					}
					
					set yOffset [expr [lindex $Lcap [expr 2*$i]] - (($nCols-1)/2 - $k)*([lindex $Lcap [expr 2*$i]])/(($nCols-1)/2)];
					set xOffset [expr [lindex $dxc [expr 2*$i]] - (($nCols-1)/2 - $k)*([lindex $dxc [expr 2*$i]])/(($nCols-1)/2)]
					node [expr ($i + 1)*1000 + (2*$k)*100 + $j] [expr $Xi - $xOffset] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];
					if {$write_model_files == 1} {
					puts $ExportID "node [expr ($i + 1)*1000 + (2*$k)*100 + $j] [expr $Xi - $xOffset] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];"
					}
					
					if {$j == 1 || $j == [expr $nEleCol + 1]} {
						lappend primaryNodesColumn [expr ($i + 1)*1000 + (2*$k - 1)*100 + $j];
						lappend primaryNodesColumn [expr ($i + 1)*1000 + (2*$k)*100 + $j];
					} else {
						lappend secondaryNodesColumn [expr ($i + 1)*1000 + (2*$k - 1)*100 + $j];
						lappend secondaryNodesColumn [expr ($i + 1)*1000 + (2*$k)*100 + $j];
					}
				}	
			}
		} else {
			set yOffset 0.0;
			
			for {set j 1} {$j <= [expr $nEleCol + 1]} {incr j 1} {
			
				set yOffset 0.0;
				node [expr ($i + 1)*1000 + $j] [expr $Xi] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];
				if {$write_model_files == 1} {
				puts $ExportID "node [expr ($i + 1)*1000 + $j] [expr $Xi] [expr $yOffset] [expr ($HColumnMax - $colHgtCurrentBent) + ($j-1)*$colHgtCurrentBent/$nEleCol];"
				}
				
				if {$j == 1 || $j == [expr $nEleCol + 1]} {
					if {$j == 1} {
						lappend fixNodesColumn [expr ($i + 1)*1000 + $j];
					}
					lappend primaryNodesColumn [expr ($i + 1)*1000 + $j];
				} else {
					lappend secondaryNodesColumn [expr ($i + 1)*1000 + $j];
				}
			}
		}
		
	}
	
	set Xi [expr $Xi + [lindex $L [expr $i + 1]]];
}

if {$nEleCol == 1} {
unset secondaryNodesColumn;
}

set ctrCol 1;
set ctrSec 1;
source "$model_files_path/createColElem.tcl";
source "$model_files_path/createBaseHingeElem.tcl";
source "$model_files_path/connectColBaseToDeck.tcl";
source "$model_files_path/integrationPointsandWeights_nIntPts_$nIntPts.tcl";
set ColumnElementList {};
for {set i 1} {$i <= $nBent} {incr i 1} {
	# # # # # # # # set Lp [expr (0.08*[lindex $HCol $i]/$mm + 0.022*($fyNom/$MPa)*($DbarSingle/$mm))*$mm];
	# # # # # # # # if {$write_model_files == 1} {
		# # # # # # # # set LPFileID [open "$model_info_directory/lengthPH.txt" "a"];
		# # # # # # # # puts $LPFileID $Lp
		# # # # # # # # close $LPFileID
	# # # # # # # # }
	set deckNode [expr ($i + 1)*100];
	
	if {[expr $nCols%2] == 0} {
	
		for {set k 1} {$k <= [expr $nCols/2]} {incr k 1} {
		
			set iNodeTag [expr ($i + 1)*1000 + (2*$k - 1)*100 + 1]
			set jNodeTag [expr ($i + 1)*1000 + (2*$k - 1)*100 + $nEleCol + 1]
			if {$write_model_files == 1} {
				set ColEleSecFileID [open "$model_info_directory/col_elem_sec_info_col_$ctrCol.txt" "w"];	
			}
			set ColElements {};
			set BHElements {};
			set secTags {};
			for {set secIter 1} {$secIter <= $num_secdef_per_col} {incr secIter 1} {
				lappend secTags [expr $ColSecTag + $ctrSec]
				set ctrSec [expr $ctrSec + 1];
			}
			if {$write_model_files == 1} {
			set elems [connectColBaseToDeck $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $transfTagColDeckConnection $maxIters $tol $ColElements \
				$existsBH $BHEleType $BHintPts [expr $BHSecTag + $ctrCol] $BHintType $BHElements $deckNode $Ubig $write_model_files $ExportID $ColEleSecFileID];
			} else {
			set elems [connectColBaseToDeck $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $transfTagColDeckConnection $maxIters $tol $ColElements \
				$existsBH $BHEleType $BHintPts [expr $BHSecTag + $ctrCol] $BHintType $BHElements $deckNode $Ubig $write_model_files];
			}
			set BHElements [lindex $elems 0];
			set ColElements [lindex $elems 1];
			lappend ColumnElementList $ColElements;
			if {$write_model_files == 1} {
				close $ColEleSecFileID
			}
			set ctrCol [expr $ctrCol + 1];
			
			
			
			
			set iNodeTag [expr ($i + 1)*1000 + (2*$k)*100 + 1]
			set jNodeTag [expr ($i + 1)*1000 + (2*$k)*100 + $nEleCol + 1]
			if {$write_model_files == 1} {
				set ColEleSecFileID [open "$model_info_directory/col_elem_sec_info_col_$ctrCol.txt" "w"];	
			}
			set ColElements {};
			set BHElements {};
			set secTags {};
			for {set secIter 1} {$secIter <= $num_secdef_per_col} {incr secIter 1} {
				lappend secTags [expr $ColSecTag + $ctrSec]
				set ctrSec [expr $ctrSec + 1];
			}
			if {$write_model_files == 1} {
			set elems [connectColBaseToDeck $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $transfTagColDeckConnection $maxIters $tol $ColElements \
				$existsBH $BHEleType $BHintPts [expr $BHSecTag + $ctrCol] $BHintType $BHElements $deckNode $Ubig $write_model_files $ExportID $ColEleSecFileID];
			} else {
			set elems [connectColBaseToDeck $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $transfTagColDeckConnection $maxIters $tol $ColElements \
				$existsBH $BHEleType $BHintPts [expr $BHSecTag + $ctrCol] $BHintType $BHElements $deckNode $Ubig $write_model_files];
			}
			set BHElements [lindex $elems 0];
			set ColElements [lindex $elems 1];
			lappend ColumnElementList $ColElements;
			if {$write_model_files == 1} {
				close $ColEleSecFileID
			}
			set ctrCol [expr $ctrCol + 1];
		}
	
	} else {
	
		set iNodeTag [expr ($i + 1)*1000 + 1]
		set jNodeTag [expr ($i + 1)*1000 + $nEleCol + 1]
		if {$write_model_files == 1} {
			set ColEleSecFileID [open "$model_info_directory/col_elem_sec_info_col_$ctrCol.txt" "w"];	
		}
		set ColElements {};
		set BHElements {};
		set secTags {};
		for {set secIter 1} {$secIter <= $num_secdef_per_col} {incr secIter 1} {
			lappend secTags [expr $ColSecTag + $ctrSec]
			set ctrSec [expr $ctrSec + 1];
		}
		if {$write_model_files == 1} {
		set elems [connectColBaseToDeck $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $transfTagColDeckConnection $maxIters $tol $ColElements \
			$existsBH $BHEleType $BHintPts [expr $BHSecTag + $ctrCol] $BHintType $BHElements $deckNode $Ubig $write_model_files $ExportID $ColEleSecFileID];
		} else {
		set elems [connectColBaseToDeck $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $transfTagColDeckConnection $maxIters $tol $ColElements \
			$existsBH $BHEleType $BHintPts [expr $BHSecTag + $ctrCol] $BHintType $BHElements $deckNode $Ubig $write_model_files];
		}
		set BHElements [lindex $elems 0];
		set ColElements [lindex $elems 1];
		lappend ColumnElementList $ColElements;
		if {$write_model_files == 1} {
			close $ColEleSecFileID
		}
		set ctrCol [expr $ctrCol + 1];
		
		
		
		for {set k 1} {$k <= [expr ($nCols-1)/2]} {incr k 1} {
		
			set iNodeTag [expr ($i + 1)*1000 + (2*$k - 1)*100 + 1]
			set jNodeTag [expr ($i + 1)*1000 + (2*$k - 1)*100 + $nEleCol + 1]
			if {$write_model_files == 1} {
				set ColEleSecFileID [open "$model_info_directory/col_elem_sec_info_col_$ctrCol.txt" "w"];	
			}
			set ColElements {};
			set BHElements {};
			set secTags {};
			for {set secIter 1} {$secIter <= $num_secdef_per_col} {incr secIter 1} {
				lappend secTags [expr $ColSecTag + $ctrSec]
				set ctrSec [expr $ctrSec + 1];
			}
			if {$write_model_files == 1} {
			set elems [connectColBaseToDeck $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $transfTagColDeckConnection $maxIters $tol $ColElements \
				$existsBH $BHEleType $BHintPts [expr $BHSecTag + $ctrCol] $BHintType $BHElements $deckNode $Ubig $write_model_files $ExportID $ColEleSecFileID];
			} else {
			set elems [connectColBaseToDeck $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $transfTagColDeckConnection $maxIters $tol $ColElements \
				$existsBH $BHEleType $BHintPts [expr $BHSecTag + $ctrCol] $BHintType $BHElements $deckNode $Ubig $write_model_files];
			}
			set BHElements [lindex $elems 0];
			set ColElements [lindex $elems 1];
			lappend ColumnElementList $ColElements;
			if {$write_model_files == 1} {
				close $ColEleSecFileID
			}
			set ctrCol [expr $ctrCol + 1];
			
			
			
			
			set iNodeTag [expr ($i + 1)*1000 + (2*$k)*100 + 1]
			set jNodeTag [expr ($i + 1)*1000 + (2*$k)*100 + $nEleCol + 1]
			if {$write_model_files == 1} {
				set ColEleSecFileID [open "$model_info_directory/col_elem_sec_info_col_$ctrCol.txt" "w"];	
			}
			set ColElements {};
			set BHElements {};
			set secTags {};
			for {set secIter 1} {$secIter <= $num_secdef_per_col} {incr secIter 1} {
				lappend secTags [expr $ColSecTag + $ctrSec]
				set ctrSec [expr $ctrSec + 1];
			}
			if {$write_model_files == 1} {
			set elems [connectColBaseToDeck $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $transfTagColDeckConnection $maxIters $tol $ColElements \
				$existsBH $BHEleType $BHintPts [expr $BHSecTag + $ctrCol] $BHintType $BHElements $deckNode $Ubig $write_model_files $ExportID $ColEleSecFileID];
			} else {
			set elems [connectColBaseToDeck $iNodeTag $jNodeTag $nEleCol $colEleType $nIntPts $locations $weights $secTags $transfTagCol $transfTagColDeckConnection $maxIters $tol $ColElements \
				$existsBH $BHEleType $BHintPts [expr $BHSecTag + $ctrCol] $BHintType $BHElements $deckNode $Ubig $write_model_files];
			}
			set BHElements [lindex $elems 0];
			set ColElements [lindex $elems 1];
			lappend ColumnElementList $ColElements;
			if {$write_model_files == 1} {
				close $ColEleSecFileID
			}
			set ctrCol [expr $ctrCol + 1];
		}
	}	
}
