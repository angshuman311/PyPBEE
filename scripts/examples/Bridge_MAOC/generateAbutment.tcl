# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------


set AbutmentFixedNodes {};
set AbutmentFreeNodes {};

set AbutmentSpringElementsLong {};
set AbutmentSpringElementsVert {};
set AbutmentSpringElementsTrans {};

if {$write_model_files == 1} {
    set matDataFileID [open "$model_info_directory/material_data.txt" "a"];
}

# ###################################################################################################################################################
# ############################################################### ABUTMENT ##########################################################################
# ###################################################################################################################################################

# #######################################
# ########## NODE ASSIGNMENT ############
# #######################################

for {set i 1} {$i <= 2} {incr i 1} {

	if {$i == 1} {
		set X 0.
		set SpanNumber 0;
	} else {
		set X $L_tot;
		set SpanNumber $nSpans;
	}

	node [expr ($SpanNumber + 1)*1000 + 100 + 0]    [expr $X+(2*$dw/10)*tan($skew*$pi/180)] [expr (-2*$dw/10)] [expr $HColumnMax + $dcg];
	node [expr ($SpanNumber + 1)*1000 + 200 + 0]    [expr $X-(2*$dw/10)*tan($skew*$pi/180)] [expr (2*$dw/10)] [expr $HColumnMax + $dcg];
	node [expr ($SpanNumber + 1)*1000 + 300 + 0]    [expr $X+(4*$dw/10)*tan($skew*$pi/180)] [expr (-4*$dw/10)] [expr $HColumnMax + $dcg];
	node [expr ($SpanNumber + 1)*1000 + 400 + 0]    [expr $X-(4*$dw/10)*tan($skew*$pi/180)] [expr (4*$dw/10)] [expr $HColumnMax + $dcg];
	node [expr ($SpanNumber + 1)*1000 + 500 + 0]    [expr $X+$dx] [expr (-$dw/2)] [expr $HColumnMax + $dcg];
	node [expr ($SpanNumber + 1)*1000 + 600 + 0]    [expr $X-$dx] [expr ($dw/2)] [expr $HColumnMax + $dcg];
    
	node [expr ($SpanNumber + 1)*1000 + 100 + 1]    [expr $X+(2*$dw/10)*tan($skew*$pi/180)] [expr (-2*$dw/10)] [expr $HColumnMax + $dcg];
	node [expr ($SpanNumber + 1)*1000 + 200 + 1]    [expr $X-(2*$dw/10)*tan($skew*$pi/180)] [expr (2*$dw/10)] [expr $HColumnMax + $dcg];
	node [expr ($SpanNumber + 1)*1000 + 300 + 1]    [expr $X+(4*$dw/10)*tan($skew*$pi/180)] [expr (-4*$dw/10)] [expr $HColumnMax + $dcg];
	node [expr ($SpanNumber + 1)*1000 + 400 + 1]    [expr $X-(4*$dw/10)*tan($skew*$pi/180)] [expr (4*$dw/10)] [expr $HColumnMax + $dcg];
	node [expr ($SpanNumber + 1)*1000 + 500 + 1]    [expr $X+$dx] [expr (-$dw/2)] [expr $HColumnMax + $dcg];
	node [expr ($SpanNumber + 1)*1000 + 600 + 1]    [expr $X-$dx] [expr ($dw/2)] [expr $HColumnMax + $dcg];

	node [expr ($SpanNumber + 1)*1000 + 1]          [expr $X] [expr 0.00] [expr $HColumnMax + $dcg];
    
	lappend AbutmentFixedNodes	[expr ($SpanNumber + 1)*1000 + 100 + 1]\
								[expr ($SpanNumber + 1)*1000 + 200 + 1]\
								[expr ($SpanNumber + 1)*1000 + 300 + 1]\
								[expr ($SpanNumber + 1)*1000 + 400 + 1]\
								[expr ($SpanNumber + 1)*1000 + 500 + 1]\
								[expr ($SpanNumber + 1)*1000 + 600 + 1]\
								[expr ($SpanNumber + 1)*1000 + 1];

	lappend AbutmentFreeNodes   [expr ($SpanNumber + 1)*1000 + 100 + 0]\
								[expr ($SpanNumber + 1)*1000 + 200 + 0]\
								[expr ($SpanNumber + 1)*1000 + 300 + 0]\
								[expr ($SpanNumber + 1)*1000 + 400 + 0]\
								[expr ($SpanNumber + 1)*1000 + 500 + 0]\
								[expr ($SpanNumber + 1)*1000 + 600 + 0];
                                
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 100 + 0]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 200 + 0]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 300 + 0]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 400 + 0]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 500 + 0]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 600 + 0]      
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 100 + 1]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 200 + 1]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 300 + 1]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 400 + 1]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 500 + 1]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 600 + 1]
    lappend allNodes [expr ($SpanNumber + 1)*1000 + 1]      
    
	if {$write_model_files == 1} {
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 100 + 0] [expr $X+(2*$dw/10)*tan($skew*$pi/180)] [expr (-2*$dw/10)] [expr $HColumnMax + $dcg];"
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 200 + 0] [expr $X-(2*$dw/10)*tan($skew*$pi/180)] [expr (2*$dw/10)] [expr $HColumnMax + $dcg];"
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 300 + 0] [expr $X+(4*$dw/10)*tan($skew*$pi/180)] [expr (-4*$dw/10)] [expr $HColumnMax + $dcg];"
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 400 + 0] [expr $X-(4*$dw/10)*tan($skew*$pi/180)] [expr (4*$dw/10)] [expr $HColumnMax + $dcg];"
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 500 + 0] [expr $X+$dx] [expr (-$dw/2)] [expr $HColumnMax + $dcg];"
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 600 + 0] [expr $X-$dx] [expr ($dw/2)] [expr $HColumnMax + $dcg];"
        puts $ExportID ""
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 100 + 1] [expr $X+(2*$dw/10)*tan($skew*$pi/180)] [expr (-2*$dw/10)] [expr $HColumnMax + $dcg];"
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 200 + 1] [expr $X-(2*$dw/10)*tan($skew*$pi/180)] [expr (2*$dw/10)] [expr $HColumnMax + $dcg];"
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 300 + 1] [expr $X+(4*$dw/10)*tan($skew*$pi/180)] [expr (-4*$dw/10)] [expr $HColumnMax + $dcg];"
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 400 + 1] [expr $X-(4*$dw/10)*tan($skew*$pi/180)] [expr (4*$dw/10)] [expr $HColumnMax + $dcg];"
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 500 + 1] [expr $X+$dx] [expr (-$dw/2)] [expr $HColumnMax + $dcg];"
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 600 + 1] [expr $X-$dx] [expr ($dw/2)] [expr $HColumnMax + $dcg];"
        puts $ExportID ""
        puts $ExportID "node [expr ($SpanNumber + 1)*1000 + 1] [expr $X] [expr 0.00] [expr $HColumnMax + $dcg];"
	}

}

# ###################################################
# ############### RIGID ELEMENTS ####################
# ###################################################

for {set i 1} {$i <= 2} {incr i 1} {

	if {$i == 1} {
		set X 0.
		set SpanNumber 0;
	} else {
		set X $L_tot;
		set SpanNumber $nSpans;
	}

    # rigidLink beam [expr ($SpanNumber + 1)*100] [expr ($SpanNumber + 1)*1000 + 100 + 0];
    # rigidLink beam [expr ($SpanNumber + 1)*100] [expr ($SpanNumber + 1)*1000 + 200 + 0];
    # rigidLink beam [expr ($SpanNumber + 1)*100] [expr ($SpanNumber + 1)*1000 + 300 + 0];
    # rigidLink beam [expr ($SpanNumber + 1)*100] [expr ($SpanNumber + 1)*1000 + 400 + 0];
    # rigidLink beam [expr ($SpanNumber + 1)*100] [expr ($SpanNumber + 1)*1000 + 500 + 0];
    # rigidLink beam [expr ($SpanNumber + 1)*100] [expr ($SpanNumber + 1)*1000 + 600 + 0];

    element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 100 + 0] [expr ($SpanNumber + 1)*1000 + 100 + 0] [expr ($SpanNumber + 1)*100]            1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;
    element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 200 + 0] [expr ($SpanNumber + 1)*1000 + 200 + 0] [expr ($SpanNumber + 1)*100]            1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;
    element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 300 + 0] [expr ($SpanNumber + 1)*1000 + 300 + 0] [expr ($SpanNumber + 1)*1000 + 100 + 0] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;
    element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 400 + 0] [expr ($SpanNumber + 1)*1000 + 400 + 0] [expr ($SpanNumber + 1)*1000 + 200 + 0] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;
    element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 500 + 0] [expr ($SpanNumber + 1)*1000 + 500 + 0] [expr ($SpanNumber + 1)*1000 + 300 + 0] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;
    element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 600 + 0] [expr ($SpanNumber + 1)*1000 + 600 + 0] [expr ($SpanNumber + 1)*1000 + 400 + 0] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;

    if {$write_model_files == 1} {
        puts $ExportID "element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 100 + 0] [expr ($SpanNumber + 1)*1000 + 100 + 0] [expr ($SpanNumber + 1)*100]            1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;"
        puts $ExportID "element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 200 + 0] [expr ($SpanNumber + 1)*1000 + 200 + 0] [expr ($SpanNumber + 1)*100]            1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;"
        puts $ExportID "element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 300 + 0] [expr ($SpanNumber + 1)*1000 + 300 + 0] [expr ($SpanNumber + 1)*1000 + 100 + 0] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;"
        puts $ExportID "element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 400 + 0] [expr ($SpanNumber + 1)*1000 + 400 + 0] [expr ($SpanNumber + 1)*1000 + 200 + 0] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;"
        puts $ExportID "element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 500 + 0] [expr ($SpanNumber + 1)*1000 + 500 + 0] [expr ($SpanNumber + 1)*1000 + 300 + 0] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;"
        puts $ExportID "element elasticBeamColumn [expr ($SpanNumber + 1)*1000 + 600 + 0] [expr ($SpanNumber + 1)*1000 + 600 + 0] [expr ($SpanNumber + 1)*1000 + 400 + 0] 1. $Ubig $Ubig 1. 1. 1. $transfTagAbut;"
    }
}
# #####################################################################################################################################################################
# ################################################################### MATERIALS FOR ABUTMENT MODELING #################################################################
# #####################################################################################################################################################################

# ##################################################################################
# ################################ SHEAR KEY #######################################
# ##################################################################################

set SK_PF1 [expr 0.45*$MN];		set SK_PD1 [expr 6.350e-03*$m];
set SK_PF2 [expr 1.8074*$MN];	set SK_PD2 [expr 95.25e-03*$m];
set SK_PF3 [expr 0.3615*$MN]; 	set SK_PD3 [expr 120.6e-03*$m];
set SK_PF4 [expr 0.3615*$MN]; 	set SK_PD4 [expr 127.0e-03*$m];

for {set i 1} {$i <= 4} {incr i 1} {

    uniaxialMaterial Concrete02 [expr 8000 + $i] [expr -$SK_PF2] [expr -$SK_PD2] [expr -$SK_PF3] [expr -$SK_PD3] 0.99 [expr 0.0001*$SK_PF2] 1.0;

    if {$write_model_files == 1} {
        puts $matDataFileID "uniaxialMaterial Concrete02 [expr 8000 + $i] [expr -$SK_PF2] [expr -$SK_PD2] [expr -$SK_PF3] [expr -$SK_PD3] 0.99 [expr 0.0001*$SK_PF2] 1.0;"
    }
}

# ##################################################################################
# ########################### EMBANKMENT BACKFILL ##################################
# ##################################################################################

set fac {0.2 0.2 0.2 0.2 0.2}
set R    [expr exp(-$skew/45.)]
set ymax [expr 0.05 * $dd]
set n_bf 5
for {set i 1} {$i <= 2} {incr i 1} {
    
    set varName       "K50_u";
	set K50_u         [expr $[append varName "_$i"]*$kip/$in];
    
    set varName       "Fult_u";
	set Fult_u         [expr $[append varName "_$i"]*$kip];
    
    set K50 [expr $K50_u * $R]
    set Fult [expr $Fult_u * $R]
    
    set C_hyp [expr 2.*$K50 - $Fult/$ymax]
    set D_hyp [expr 2.*($K50/$Fult - 1./$ymax)]
    
    for {set j 1} {$j <= $n_bf} {incr j 1} {
        set C_hyp_j [expr $C_hyp * [lindex $fac [expr $n_bf - $j]]]
        uniaxialMaterial HyperbolicGapMaterial [expr 100 + ($i - 1)*$n_bf + $j] $C_hyp_j $C_hyp_j 1.0 [expr -$C_hyp_j/$D_hyp] -$gapL;
        if {$write_model_files == 1} {
            puts $matDataFileID "uniaxialMaterial HyperbolicGapMaterial [expr 100 + ($i - 1)*$n_bf + $j] $C_hyp_j $C_hyp_j 1.0 [expr -$C_hyp_j/$D_hyp] -$gapL;"
        }
    }
}

# ##################################################################################
# ################################### OTHERS #######################################
# ##################################################################################

uniaxialMaterial Elastic 	  20000 	[expr 1.e-10];
uniaxialMaterial Elastic           203       [expr 1.e6];

#
# ################################################################################################################################################################################################################################################################
# ################################################################################################################################################################################################################################################################
# ################################################################################################################################################################################################################################################################
#

set doRayleighHyp 0;
set doRayleighFict 1;
set doRayleighBearing 0;
set doRayleighSK 0;

# 1st Abutment
element zeroLength 1301 1301 1300 -mat 105 -dir 1 -doRayleigh $doRayleighHyp -orient [expr cos($skew*$pi/180)] [expr sin($skew*$pi/180)] 0. 0. 0. 1.;
element zeroLength 1101 1101 1100 -mat 104 -dir 1 -doRayleigh $doRayleighHyp -orient [expr cos($skew*$pi/180)] [expr sin($skew*$pi/180)] 0. 0. 0. 1.;
element zeroLength 1001 1001 100  -mat 103 -dir 1 -doRayleigh $doRayleighHyp -orient [expr cos($skew*$pi/180)] [expr sin($skew*$pi/180)] 0. 0. 0. 1.;
element zeroLength 1201 1201 1200 -mat 102 -dir 1 -doRayleigh $doRayleighHyp -orient [expr cos($skew*$pi/180)] [expr sin($skew*$pi/180)] 0. 0. 0. 1.;
element zeroLength 1401 1401 1400 -mat 101 -dir 1 -doRayleigh $doRayleighHyp -orient [expr cos($skew*$pi/180)] [expr sin($skew*$pi/180)] 0. 0. 0. 1.;

element elastomericBearingPlasticity 1303 1301 1300 $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;
element elastomericBearingPlasticity 1103 1101 1100 $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;
element elastomericBearingPlasticity 1003 1001 100  $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;
element elastomericBearingPlasticity 1203 1201 1200 $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;
element elastomericBearingPlasticity 1403 1401 1400 $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;

# Last Abutment
element zeroLength [expr ($nSpans + 1)*1000 + 300 + 1] [expr ($nSpans + 1)*1000 + 300 + 1] [expr ($nSpans + 1)*1000 + 300 + 0] -mat 106 -dir 1 -doRayleigh $doRayleighHyp -orient [expr -cos($skew*$pi/180)] [expr -sin($skew*$pi/180)] 0. 0. 0. 1.;
element zeroLength [expr ($nSpans + 1)*1000 + 100 + 1] [expr ($nSpans + 1)*1000 + 100 + 1] [expr ($nSpans + 1)*1000 + 100 + 0] -mat 107 -dir 1 -doRayleigh $doRayleighHyp -orient [expr -cos($skew*$pi/180)] [expr -sin($skew*$pi/180)] 0. 0. 0. 1.;
element zeroLength [expr ($nSpans + 1)*1000 + 1]       [expr ($nSpans + 1)*1000 + 1]       [expr ($nSpans + 1)*100]            -mat 108 -dir 1 -doRayleigh $doRayleighHyp -orient [expr -cos($skew*$pi/180)] [expr -sin($skew*$pi/180)] 0. 0. 0. 1.;
element zeroLength [expr ($nSpans + 1)*1000 + 200 + 1] [expr ($nSpans + 1)*1000 + 200 + 1] [expr ($nSpans + 1)*1000 + 200 + 0] -mat 109 -dir 1 -doRayleigh $doRayleighHyp -orient [expr -cos($skew*$pi/180)] [expr -sin($skew*$pi/180)] 0. 0. 0. 1.;
element zeroLength [expr ($nSpans + 1)*1000 + 400 + 1] [expr ($nSpans + 1)*1000 + 400 + 1] [expr ($nSpans + 1)*1000 + 400 + 0] -mat 110 -dir 1 -doRayleigh $doRayleighHyp -orient [expr -cos($skew*$pi/180)] [expr -sin($skew*$pi/180)] 0. 0. 0. 1.;

element elastomericBearingPlasticity [expr ($nSpans + 1)*1000 + 300 + 3] [expr ($nSpans + 1)*1000 + 300 + 1] [expr ($nSpans + 1)*1000 + 300 + 0] $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;
element elastomericBearingPlasticity [expr ($nSpans + 1)*1000 + 100 + 3] [expr ($nSpans + 1)*1000 + 100 + 1] [expr ($nSpans + 1)*1000 + 100 + 0] $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;
element elastomericBearingPlasticity [expr ($nSpans + 1)*1000 + 3]       [expr ($nSpans + 1)*1000 + 1]       [expr ($nSpans + 1)*100]            $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;
element elastomericBearingPlasticity [expr ($nSpans + 1)*1000 + 200 + 3] [expr ($nSpans + 1)*1000 + 200 + 1] [expr ($nSpans + 1)*1000 + 200 + 0] $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;
element elastomericBearingPlasticity [expr ($nSpans + 1)*1000 + 400 + 3] [expr ($nSpans + 1)*1000 + 400 + 1] [expr ($nSpans + 1)*1000 + 400 + 0] $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;

# Shear Keys
element zeroLength 1501 1501 1500 -mat 8001 -dir 2 -doRayleigh $doRayleighSK -orient  1 0 0 0 1 0;
element zeroLength 1601 1601 1600 -mat 8002 -dir 2 -doRayleigh $doRayleighSK -orient -1 0 0 0 -1 0;
element zeroLength [expr ($nSpans + 1)*1000 + 500 + 1] [expr ($nSpans + 1)*1000 + 500 + 1] [expr ($nSpans + 1)*1000 + 500 + 0] -mat 8003 -dir 2 -doRayleigh $doRayleighSK -orient 1 0 0 0 1 0;
element zeroLength [expr ($nSpans + 1)*1000 + 600 + 1] [expr ($nSpans + 1)*1000 + 600 + 1] [expr ($nSpans + 1)*1000 + 600 + 0] -mat 8004 -dir 2 -doRayleigh $doRayleighSK -orient -1 0 0 0 -1 0;

if {$write_model_files == 1} {
set shearKeySpringElemFileID [open "$model_info_directory/shear_key_spring_elem_tags.txt" "w"];
puts $shearKeySpringElemFileID "1501";
puts $shearKeySpringElemFileID "1601";
puts $shearKeySpringElemFileID "[expr ($nSpans + 1)*1000 + 500 + 1]";
puts $shearKeySpringElemFileID "[expr ($nSpans + 1)*1000 + 600 + 1]";
close $shearKeySpringElemFileID;
}

lappend  AbutmentSpringElementsLong 1301\
                                    1101\
                                    1001\
                                    1201\
                                    1401\
									[expr ($nSpans + 1)*1000 + 300 + 1]\
                                    [expr ($nSpans + 1)*1000 + 100 + 1]\
                                    [expr ($nSpans + 1)*1000 + 1]      \
                                    [expr ($nSpans + 1)*1000 + 200 + 1]\
                                    [expr ($nSpans + 1)*1000 + 400 + 1];


lappend  AbutmentSpringElementsTrans 1501\
                                     1601\
                                     [expr ($nSpans + 1)*1000 + 500 + 1]\
                                     [expr ($nSpans + 1)*1000 + 600 + 1];


lappend  AbutmentSpringElementsVert 1303\
                                    1103\
                                    1003\
                                    1203\
                                    1403\
									[expr ($nSpans + 1)*1000 + 300 + 3]\
                                    [expr ($nSpans + 1)*1000 + 100 + 3]\
                                    [expr ($nSpans + 1)*1000 + 3]      \
                                    [expr ($nSpans + 1)*1000 + 200 + 3]\
                                    [expr ($nSpans + 1)*1000 + 400 + 3];

if {$write_model_files == 1} {
puts $ExportID "element zeroLength 1301 1301 1300 -mat 105 -dir 1 -doRayleigh $doRayleighHyp -orient [expr cos($skew*$pi/180)] [expr sin($skew*$pi/180)] 0. 0. 0. 1.;"
puts $ExportID "element zeroLength 1101 1101 1100 -mat 104 -dir 1 -doRayleigh $doRayleighHyp -orient [expr cos($skew*$pi/180)] [expr sin($skew*$pi/180)] 0. 0. 0. 1.;"
puts $ExportID "element zeroLength 1001 1001 100  -mat 103 -dir 1 -doRayleigh $doRayleighHyp -orient [expr cos($skew*$pi/180)] [expr sin($skew*$pi/180)] 0. 0. 0. 1.;"
puts $ExportID "element zeroLength 1201 1201 1200 -mat 102 -dir 1 -doRayleigh $doRayleighHyp -orient [expr cos($skew*$pi/180)] [expr sin($skew*$pi/180)] 0. 0. 0. 1.;"
puts $ExportID "element zeroLength 1401 1401 1400 -mat 101 -dir 1 -doRayleigh $doRayleighHyp -orient [expr cos($skew*$pi/180)] [expr sin($skew*$pi/180)] 0. 0. 0. 1.;"

puts $ExportID "element elastomericBearingPlasticity 1303 1301 1300 $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;"
puts $ExportID "element elastomericBearingPlasticity 1103 1101 1100 $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;"
puts $ExportID "element elastomericBearingPlasticity 1003 1001 100  $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;"
puts $ExportID "element elastomericBearingPlasticity 1203 1201 1200 $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;"
puts $ExportID "element elastomericBearingPlasticity 1403 1401 1400 $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;"

puts $ExportID "element zeroLength [expr ($nSpans + 1)*1000 + 300 + 1] [expr ($nSpans + 1)*1000 + 300 + 1] [expr ($nSpans + 1)*1000 + 300 + 0] -mat 106 -dir 1 -doRayleigh $doRayleighHyp -orient [expr -cos($skew*$pi/180)] [expr -sin($skew*$pi/180)] 0. 0. 0. 1.;"
puts $ExportID "element zeroLength [expr ($nSpans + 1)*1000 + 100 + 1] [expr ($nSpans + 1)*1000 + 100 + 1] [expr ($nSpans + 1)*1000 + 100 + 0] -mat 107 -dir 1 -doRayleigh $doRayleighHyp -orient [expr -cos($skew*$pi/180)] [expr -sin($skew*$pi/180)] 0. 0. 0. 1.;"
puts $ExportID "element zeroLength [expr ($nSpans + 1)*1000 + 1]       [expr ($nSpans + 1)*1000 + 1]       [expr ($nSpans + 1)*100]            -mat 108 -dir 1 -doRayleigh $doRayleighHyp -orient [expr -cos($skew*$pi/180)] [expr -sin($skew*$pi/180)] 0. 0. 0. 1.;"
puts $ExportID "element zeroLength [expr ($nSpans + 1)*1000 + 200 + 1] [expr ($nSpans + 1)*1000 + 200 + 1] [expr ($nSpans + 1)*1000 + 200 + 0] -mat 109 -dir 1 -doRayleigh $doRayleighHyp -orient [expr -cos($skew*$pi/180)] [expr -sin($skew*$pi/180)] 0. 0. 0. 1.;"
puts $ExportID "element zeroLength [expr ($nSpans + 1)*1000 + 400 + 1] [expr ($nSpans + 1)*1000 + 400 + 1] [expr ($nSpans + 1)*1000 + 400 + 0] -mat 110 -dir 1 -doRayleigh $doRayleighHyp -orient [expr -cos($skew*$pi/180)] [expr -sin($skew*$pi/180)] 0. 0. 0. 1.;"

puts $ExportID "element elastomericBearingPlasticity [expr ($nSpans + 1)*1000 + 300 + 3] [expr ($nSpans + 1)*1000 + 300 + 1] [expr ($nSpans + 1)*1000 + 300 + 0] $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;"
puts $ExportID "element elastomericBearingPlasticity [expr ($nSpans + 1)*1000 + 100 + 3] [expr ($nSpans + 1)*1000 + 100 + 1] [expr ($nSpans + 1)*1000 + 100 + 0] $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;"
puts $ExportID "element elastomericBearingPlasticity [expr ($nSpans + 1)*1000 + 3]       [expr ($nSpans + 1)*1000 + 1]       [expr ($nSpans + 1)*100]            $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;"
puts $ExportID "element elastomericBearingPlasticity [expr ($nSpans + 1)*1000 + 200 + 3] [expr ($nSpans + 1)*1000 + 200 + 1] [expr ($nSpans + 1)*1000 + 200 + 0] $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;"
puts $ExportID "element elastomericBearingPlasticity [expr ($nSpans + 1)*1000 + 400 + 3] [expr ($nSpans + 1)*1000 + 400 + 1] [expr ($nSpans + 1)*1000 + 400 + 0] $kInit1 $qd1 0.1 0.0 3.0 -P 203 -T 20000 -My 20000 -Mz 20000 -orient 0. 0. 1. 1. 0. 0. -doRayleigh $doRayleighBearing;"

puts $ExportID "element zeroLength 1501 1501 1500 -mat 202 -dir 2 -doRayleigh $doRayleighSK -orient  1 0 0 0 1 0;"
puts $ExportID "element zeroLength 1601 1601 1600 -mat 202 -dir 2 -doRayleigh $doRayleighSK -orient -1 0 0 0 -1 0;"
puts $ExportID "element zeroLength [expr ($nSpans + 1)*1000 + 500 + 1] [expr ($nSpans + 1)*1000 + 500 + 1] [expr ($nSpans + 1)*1000 + 500 + 0] -mat 202 -dir 2 -doRayleigh $doRayleighSK -orient 1 0 0 0 1 0;"
puts $ExportID "element zeroLength [expr ($nSpans + 1)*1000 + 600 + 1] [expr ($nSpans + 1)*1000 + 600 + 1] [expr ($nSpans + 1)*1000 + 600 + 0] -mat 202 -dir 2 -doRayleigh $doRayleighSK -orient -1 0 0 0 -1 0;"
}

if {$write_model_files == 1} {
close $matDataFileID
}
