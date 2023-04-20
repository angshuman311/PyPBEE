proc min_list list {
    set index 0
    set min_index $index
    set min_val [lindex $list 0]
    foreach val $list {
        if {$val < $min_val} {
            set min_index $index
            set min_val $val
        }
        incr index
    }
    return [list $min_val $min_index]
}
