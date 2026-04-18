proc evaluate_edp_2 strain_vals {
    set temp [max_list $strain_vals]
    return [lindex $temp 0]
}
