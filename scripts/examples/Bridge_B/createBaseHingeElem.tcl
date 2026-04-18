# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------


# createBaseHingeElem $iNodeTag $BHEleType $BHintPts $BHSecTag $transfTagCol $BHintType $BHElements $ExportID
proc createBaseHingeElem {iNodeTag BHEleType BHintPts BHSecTag transfTagCol BHintType BHElements bond BondSecTag skewr ExportID} {
	
    if {$bond == 1} {
        set BondSecTagThis [expr $BondSecTag + ($BHSecTag - ($BHSecTag/1000)*1000)]
        element $BHEleType [expr $iNodeTag - 1] [expr $iNodeTag - 1] $iNodeTag $BHintPts $BHSecTag $transfTagCol -integration $BHintType;
        element zeroLengthSection [expr ($iNodeTag - 1)*10] [expr ($iNodeTag - 1)*10] [expr $iNodeTag - 1] $BondSecTagThis -orient 0 0 1 [expr -1.*cos($skewr)] [expr -1.*sin($skewr)] 0. -doRayleigh 0;
        if {$ExportID != -1} {
            puts $ExportID "element $BHEleType [expr $iNodeTag - 1] [expr $iNodeTag - 1] $iNodeTag $BHintPts $BHSecTag $transfTagCol -integration $BHintType;"
            puts $ExportID "element zeroLengthSection [expr ($iNodeTag - 1)*10] [expr ($iNodeTag - 1)*10] [expr $iNodeTag - 1] $BondSecTagThis -orient 0 0 1 [expr -1.*cos($skewr)] [expr -1.*sin($skewr)] 0. -doRayleigh 0;"
        }
    } else {
        element $BHEleType [expr $iNodeTag - 1] [expr $iNodeTag - 1] $iNodeTag $BHintPts $BHSecTag $transfTagCol -integration $BHintType;
        if {$ExportID != -1} {
            puts $ExportID "element $BHEleType [expr $iNodeTag - 1] [expr $iNodeTag - 1] $iNodeTag $BHintPts $BHSecTag $transfTagCol -integration $BHintType;"
        }
    }
	lappend BHElements [expr $iNodeTag - 1];
	return $BHElements;
}
