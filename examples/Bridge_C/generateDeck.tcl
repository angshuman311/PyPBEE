# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

set allNodes {};

#Global X - along longitudinal direction of bridge 
#Global Y - along transverse direction of bridge
#Global Z - along height of bridge

#Local X - along longitudinal direction of bridge - along Global X
#Local Y - along transverse direction of bridge - along Global Y
#Local Z - along height of bridge - along Global Z

#Node numbering for a given span starts with (span number)*100 and increases in increment of 1 along the span till the first node of the next span is reached

set primaryNodesDeck {};
set secondaryNodesDeck {};
set Xi 0;		# abscissa of primary nodes  
set xi 0;		# abscissa of secondary nodes 
for {set i 1} {$i <= [expr $nSpans + 1]} {incr i 1} {
	node [expr $i*100] $Xi 0. [expr $HColumnMax + $dcg];
	if {$write_model_files == 1} {
	puts $ExportID "node [expr $i*100] $Xi 0. [expr $HColumnMax + $dcg]";
	}
	lappend primaryNodesDeck [expr $i*100]
	lappend allNodes [expr $i*100]
	set xi $Xi;
	if {$i < [expr $nSpans + 1]} {
		for {set j 1} {$j < $nEleSpan} {incr j 1} {
			set xi [expr $xi + [lindex $L $i]/$nEleSpan];
			node [expr $i*100 + $j] $xi 0. [expr $HColumnMax + $dcg];
			if {$write_model_files == 1} {
			puts $ExportID "node [expr $i*100 + $j] $xi 0. [expr $HColumnMax + $dcg]";
			}
			lappend secondaryNodesDeck [expr $i*100 + $j];
			lappend allNodes [expr $i*100 + $j];
		}
		set Xi [expr $Xi + [lindex $L $i]];
	}
}


set deckE [expr 29400.*$MN/$m**2];
set deckG [expr $deckE/(2.*(1 + $Uc))];

set DeckElements {};
for {set i 1} {$i < [expr $nSpans + 1]} {incr i 1} {
		for {set j 1} {$j < $nEleSpan} {incr j 1} {
			element elasticBeamColumn [expr $i*100 + $j] [expr $i*100 + $j - 1] [expr $i*100 + $j] $Adeck $deckE $deckG $Jdeck $Iydeck $Izdeck $transfTagDeck; 
			if {$write_model_files == 1} {
			puts $ExportID "element elasticBeamColumn [expr $i*100 + $j] [expr $i*100 + $j - 1] [expr $i*100 + $j] $Adeck $deckE $deckG $Jdeck $Iydeck $Izdeck $transfTagDeck;";
			}
			lappend DeckElements [expr $i*100 + $j]
		}
		element elasticBeamColumn [expr $i*100 + $j] [expr $i*100 + $j - 1] [expr ($i+1)*100] $Adeck $deckE $deckG $Jdeck $Iydeck $Izdeck $transfTagDeck;
		if {$write_model_files == 1} {
		puts $ExportID "element elasticBeamColumn [expr $i*100 + $j] [expr $i*100 + $j - 1] [expr ($i+1)*100] $Adeck $deckE $deckG $Jdeck $Iydeck $Izdeck $transfTagDeck;";
		}
		lappend DeckElements [expr $i*100 + $j]
}

if {$write_model_files == 1} {
set DeckNodeFileID [open "$model_info_directory/primary_nodes_deck.txt" "w"];
puts $DeckNodeFileID $primaryNodesDeck
close $DeckNodeFileID
}