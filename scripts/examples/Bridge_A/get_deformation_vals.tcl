proc get_deformation_vals file_path {
    set data_fid [open $file_path "r"]
    set data [read $data_fid]
    close $data_fid
    set data_new [split $data "\n"]
    set deformation_vals {}
    for {set k 0} {$k <= [expr [llength $data_new] - 2]} {incr k 1} {
        set data_t [lindex $data_new $k]
        lappend deformation_vals [lindex $data_t [expr 1]]
    }
    return $deformation_vals
}
