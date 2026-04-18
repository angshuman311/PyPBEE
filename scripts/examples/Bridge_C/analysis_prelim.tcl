# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

set omega {};
set F {};
set T {};
set PI [expr 2.0*asin(1.0)];
set lambda [eigen -fullGenLapack $modes];
# set lambda [eigen $modes];
foreach lam $lambda {
         lappend omega [expr sqrt($lam)];
         lappend F [expr sqrt($lam)/(2*$PI)];
         lappend T [expr (2.0*$PI)/sqrt($lam)];
}

######################################################################################
############################## SAVE TIME PERIODS #####################################
######################################################################################

set PeriodFile [open "$model_info_directory/periods.txt" "w"];
puts $PeriodFile "$T";
close $PeriodFile;

######################################################################################
############################### SAVE MODE SHAPES #####################################
######################################################################################

for {set iMode 1} {$iMode <= $modes} {incr iMode} {

	set modeFile [open "$model_info_directory/mode_shape_$iMode.txt" "w"];
    
    for {set i 0} {$i < [llength $allNodes]} {incr i 1} {
        puts $modeFile "[lindex $allNodes $i] [nodeEigenvector [lindex $allNodes $i] $iMode]";
    }

	close $modeFile;
}

puts "Modal Analysis COMPLETE!"
