proc max_list list {
    set index 0
    set max_index $index
    set max_val [lindex $list 0]
    foreach val $list {
        if {$val > $max_val} {
            set max_index $index
            set max_val $val
        }
        incr index
    }
    return [list $max_val $max_index]
}
