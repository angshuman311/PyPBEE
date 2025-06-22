# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------


# createCirclePatches $matTag $yCenter $zCenter $radiusInner $radiusOuter $minFiber $maxFiber $startAng $endAng ?$exportID?
proc createCirclePatches {matTag yCenter zCenter radiusInner radiusOuter \
	minFiber maxFiber startAng endAng exportID} {
	
  # Set constant PI
	set PI [expr 2.0 * asin(1.0)];

  # Determine number of total concentric rings and ring width;
	set numCirc [expr int(($radiusOuter - $radiusInner) / $maxFiber + 0.5)];
	set ringWidth [expr ($radiusOuter - $radiusInner) / $numCirc];

  # Determine ratio of min fiber to max fiber;
	set fibRatio [expr $minFiber / $maxFiber];

  # Do for loop over concentric rings to determine where to put patches;
	set jRadiusOuter $radiusOuter;
	set jRadiusCutoff [expr $radiusOuter * $fibRatio];
	set jRadiusInner $radiusInner;
	for {set iRing $numCirc} {$iRing >= 1} {incr iRing -1} {
		set iRadius [expr $radiusInner + $iRing * $ringWidth];

	  # Use fiber ratio to determine when to place a patch;
		if {$iRadius < $jRadiusCutoff} {

		  # Set inner radius and subdivisions of circular patch;
			set jRadiusInner $iRadius;
			set jNumC [expr max(int($jRadiusOuter * ($endAng - $startAng) * \
				($PI / 180.0) / $maxFiber + 0.5), 1)];
			set jNumR [expr max(int(($jRadiusOuter - $jRadiusInner) / \
				$ringWidth + 0.5), 1)];

		  # Create patch; check first to see if default angles are used;
			patch circ $matTag $jNumC $jNumR $yCenter $zCenter $jRadiusInner \
				$jRadiusOuter $startAng $endAng;

		  # If channel is provided, export information to file;
			if {$exportID != -1} {
				puts $exportID "	patch circ $matTag $jNumC $jNumR $yCenter\
					$zCenter $jRadiusInner $jRadiusOuter $startAng $endAng";
			}

			set jRadiusOuter $jRadiusInner;
			set jRadiusCutoff [expr $jRadiusOuter * $fibRatio];
		}

		if {$iRing == 1} {
		
		  # Set inner radius and subdivisions of circular patch;
			set jNumC [expr max(int($jRadiusOuter * ($endAng - $startAng) * \
				($PI / 180.0) / $maxFiber + 0.5), 1)];
			set jNumR [expr max(int(($jRadiusOuter - $radiusInner) / \
				$ringWidth + 0.5), 1)];

		  # Create patch; check first to see if default angles are used;
			patch circ $matTag $jNumC $jNumR $yCenter $zCenter \
				$radiusInner $jRadiusOuter $startAng $endAng;

		  # If channel is provided, export information to file;
			if {$exportID != -1} {
				puts $exportID "	patch circ $matTag $jNumC $jNumR $yCenter\
					$zCenter $radiusInner $jRadiusOuter $startAng $endAng\n";
			}
		}
	}
}