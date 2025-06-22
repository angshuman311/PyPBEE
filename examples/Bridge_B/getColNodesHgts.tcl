# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

proc getColNodesHgts {colEleDistPattern colEleDistHgts} {
    
    set colNodesHgts {};
    
    lappend colNodesHgts 0.
    
    set n_e_0 [lindex $colEleDistPattern 0]
    set h_e_0 [lindex $colEleDistHgts 0]
    for {set i 1} {$i <= $n_e_0} {incr i 1} {
        lappend colNodesHgts [expr ($i)*($h_e_0)/$n_e_0]
    }
    
    set n_e_1 [lindex $colEleDistPattern 1]
    set h_e_1 [lindex $colEleDistHgts 1]
    for {set i 1} {$i <= $n_e_1} {incr i 1} {
        lappend colNodesHgts [expr $h_e_0 + ($i)*($h_e_1)/$n_e_1]
    }
    
    set n_e_2 [lindex $colEleDistPattern 2]
    set h_e_2 [lindex $colEleDistHgts 2]
    for {set i 1} {$i <= $n_e_2} {incr i 1} {
        lappend colNodesHgts [expr $h_e_0 + $h_e_1 + ($i)*($h_e_2)/$n_e_2]
    }

	return $colNodesHgts;
}
