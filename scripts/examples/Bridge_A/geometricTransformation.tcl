# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------


set transfTagDeck 1000;
set transfTagCol 2000;
set transfTagColDeckConnection 3000;
set transfTagAbut 4000;

#geomTransf PDelta $transfTag     $vecxzX $vecxzY $vecxzZ <-jntOffset $dXi $dYi $dZi $dXj $dYj $dZj>
 geomTransf $gT  $transfTagDeck                                       0       0       1;
 geomTransf $gT  $transfTagAbut                                       0       0       1;
 geomTransf $gT  $transfTagCol                                        [expr -1.*sin($skewr)]    [expr cos($skewr)]    0.;
 geomTransf $gT  $transfTagColDeckConnection                          [expr cos($skewr)]        [expr sin($skewr)]    0.;
 
if {$write_model_files == 1} {
	puts $ExportID "geomTransf $gT  $transfTagDeck                                       0       0       1;"
	puts $ExportID "geomTransf $gT  $transfTagAbut                                       0       0       1;"
	puts $ExportID "geomTransf $gT  $transfTagCol                                        [expr -1.*sin($skewr)]    [expr cos($skewr)]    0.;"
	puts $ExportID "geomTransf $gT  $transfTagColDeckConnection                          [expr cos($skewr)]        [expr sin($skewr)]    0.;"
}