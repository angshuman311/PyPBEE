proc evaluate_edp_4 deformation_vals {
    set temp [min_list $deformation_vals]
    return [expr -([lindex $temp 0])]
}
