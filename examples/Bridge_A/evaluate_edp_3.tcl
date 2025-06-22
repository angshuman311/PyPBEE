proc evaluate_edp_3 strain_vals {
    set temp [max_list $strain_vals]
    set max_strain [lindex $temp 0]
    set strain_vals [lrange $strain_vals [lindex $temp 1] end]
    set temp [min_list $strain_vals]
    set min_strain [lindex $temp 0]
    return [expr $max_strain - $min_strain]
}
