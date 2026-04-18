# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

# UNITS -----------------------------------------------------------------------------
set in 1.; 				# define basic units -- output units
set kip 1.; 			# define basic units -- output units
set sec 1.; 			# define basic units -- output units
set LunitTXT "inch";		# define basic-unit text for output
set FunitTXT "kip";			# define basic-unit text for output
set TunitTXT "sec";			# define basic-unit text for output
set ft 	[expr 12.*$in]; 		# define engineering units
set ksi [expr $kip/pow($in,2)];
set psi [expr $ksi/1000.];
set lbf [expr $psi*$in*$in];		# pounds force
set pcf [expr $lbf/pow($ft,3)];		# pounds per cubic foot
set psf [expr $lbf/pow($ft,2)];		# pounds per square foot
set in2 [expr $in*$in]; 			# inch^2
set in4 [expr $in*$in*$in*$in]; 	# inch^4
set cm 	[expr $in/2.54];		# centimeter, needed for displacement input in MultipleSupport excitation
set cmsec2 [expr $cm/pow($sec,2)];	# cm/sec^2, needed for some ground accelerations
set m [expr $cm*100];			# meter
set mm [expr $cm/10];           # millimeter
set mm2 [expr $mm*$mm];         # millimeter^2
set kN [expr 0.2247*$kip];      # kilo newton
set N [expr 1.e-3*$kN];	        # newton
set MN [expr 1.e6*$N];	        # mega newton
set MPa [expr 0.1450*$ksi];     # mega pascal
set GPa [expr 1000*$MPa];     	# giga pascal
set pi [expr 2*asin(1.0)]; 		# define constants
set g [expr 32.2*$ft/pow($sec,2)]; 	# gravitational acceleration
set Ubig 10.e10; 				# a really large number
set Usmall [expr 1/$Ubig]; 		# a really small number

puts "Basic Units - $LunitTXT, $FunitTXT, $TunitTXT"