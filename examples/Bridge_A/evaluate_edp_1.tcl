proc evaluate_edp_1 strain_vals {
    set temp [min_list $strain_vals]
    return [expr -([lindex $temp 0])]
}
