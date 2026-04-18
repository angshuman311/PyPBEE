# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

if {$write_model_files == 1} {
set matDataFileID [open "$model_info_directory/material_data.txt" "w"];
set colRebarMatTagInfoFileID [open "$model_info_directory/col_rebar_mat_info.txt" "w"];
}

set IDconcCore  1000;
set IDconcCover 2000;
set IDSteel     3000;
if {$bond == 1} {
    set IDconcCoreBond 10000;
    set IDconcCoverBond 20000;
    set IDSteelBond 30000;
}
set IDShear 	4000;
set IDTorsion 	5000;
set IDShearBH 	6000;
set IDTorsionBH 7000;
set ID_SSI 9000;
set ID_Rigid 99999;

uniaxialMaterial Elastic $ID_Rigid $Ubig

set wconc         [expr $wconc_all*$kip/($in**3)];                                    
set mconc	      [expr ($wconc)/$g];                           # Normal Mass Weight per Volume

set ctr 1;
for {set i 1} {$i <= [expr $nCols * $nBent]} {incr i 1} {
	
	for {set j 1} {$j <= [expr $num_secdef_per_col]} {incr j 1} {
		
		# if j == 1 or j == num_secdef_per_col, only then it's a possible plastic hinge
		if {$write_model_files == 1} {
			if {$j == 1} {
				set predFileID1 [open "$model_info_directory/predictor_info_col_rebar_strain_damage_col_$i\_edge_1.txt" "w"];
			}
			if {$j == $num_secdef_per_col} {
				set predFileID2 [open "$model_info_directory/predictor_info_col_rebar_strain_damage_col_$i\_edge_2.txt" "w"];
			}
		}
		
		set varName       "fpc";
		set fpc           [expr $[append varName "_col_$i\_secdef_$j"]*$ksi];                                     # CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
		
		set varName       "Ec";
		set Ec            [expr $[append varName "_col_$i\_secdef_$j"]*$ksi];
		set Uc            0.2;                                                         		   # Poisson's ratio
		set Gc            [expr $Ec/(2*(1+$Uc))];                                      		   # Shear Modulus of Elasticity
		
		# Properties of confined and unconfined concrete
		
		set randCtr 	[expr int(rand()*($num_bar_clusters - 1 + 1)) + 1]
		set randCtr     1
		set varName		"fy";
		set fyH         [expr $[append varName "_col_$i\_secdef_$j\_barcluster_$randCtr"]*$ksi];
		set varName		"Es";
		set Es         [expr $[append varName "_col_$i\_secdef_$j\_barcluster_$randCtr"]*$ksi];
		
		set fl			[expr 2.*$AbarH1*$fyH/($ds1*$sH1)];
		set Ke 			[expr (1. - $sHprime1/(2.*$ds1))/(1. - $rho_cc1)];
		set fl_prime 	[expr $fl*$Ke];
		set Kc  		[expr 9./sqrt(1. + 8.*$fl_prime/$fpc) - 2.]
		set n 			[expr 0.8 + $fpc/(2500.*$psi)];
		set epc 		[expr $fpc/$Ec*($n/($n-1.))];
		set r 			[expr 1./(1. - ($fpc/$epc)/$Ec)];
		set ecu 		[expr 0.004];
		set fcu 		[expr $r*(($ecu/$epc)/($r - 1. + ($ecu/$epc)**$r))*$fpc];
		set slope_UC	[expr ($fpc - $fcu)/($epc - $ecu)];
		set ecu_UC		[expr $epc - ($fpc - 0.2*$fpc)/$slope_UC];
		
		set fpcc 		[expr $fpc + $Kc*$fl_prime];
		set Ecc 		[expr (40000.*sqrt($fpcc/$psi) + 1.e6)*$psi];
		set Gcc          [expr $Ecc/(2*(1+$Uc))];                                      		   # Shear Modulus of Elasticity
		
        if {$j == 1} {
            set Gcc1 $Gcc
        }
        
		set epcc 		[expr $epc*(1. + 5.*(($fpcc/$fpc) - 1.))];
		set epcc1 		[expr $fpcc*$epc/$fpc];
		set rc 			[expr 1./(1. - ($fpcc/$epcc)/$Ecc)];
		set rho_s 		[expr ($AbarH1*$pi*$ds1)/($pi*($ds1**2)*$sH1/4.)]; 					   # Volumetric confinement ratio
		set eccu 		[expr 0.004 + 3.*$Ke*$rho_s];
		set fccu 		[expr $rc*(($eccu/$epcc)/($rc - 1. + ($eccu/$epcc)**$rc))*$fpcc];
		set slope_C		[expr ($fpcc - $fccu)/($epcc1 - $eccu)];
		set eccu_C		[expr $epcc1 - ($fpcc - 0.2*$fpcc)/$slope_C];
		
		uniaxialMaterial Concrete01 [expr $IDconcCore + $ctr]  [expr -$fpcc] [expr -$epcc1] [expr -0.2*$fpcc] [expr -$eccu_C];       # column core-concrete
		uniaxialMaterial Concrete01 [expr $IDconcCover + $ctr] [expr -$fpc] [expr -$epc] [expr -0.2*$fpc] [expr -$ecu_UC];       	   # column cover-concrete
        if {$bond == 1} {
            if {$j == 1 || $j == $num_secdef_per_col} {
                uniaxialMaterial Concrete01 [expr $IDconcCoreBond + $ctr]  [expr -$fpcc] [expr -$epcc1] [expr -0.8*$fpcc] [expr -$eccu_C];       # column core-concrete
                uniaxialMaterial Concrete01 [expr $IDconcCoverBond + $ctr] [expr -$fpc] [expr -$epc] [expr -0.8*$fpc] [expr -$ecu_UC];       	   # column cover-concrete
            }
        }
		if {$write_model_files == 1} {
			puts $matDataFileID "uniaxialMaterial Concrete01 [expr $IDconcCore + $ctr]  [expr -$fpcc] [expr -$epcc1] [expr -0.2*$fpcc] [expr -$eccu_C];"
			puts $matDataFileID "uniaxialMaterial Concrete01 [expr $IDconcCover + $ctr] [expr -$fpc] [expr -$epc] [expr -0.2*$fpc] [expr -$ecu_UC];"
            if {$bond == 1} {
                if {$j == 1 || $j == $num_secdef_per_col} {
                    puts $matDataFileID "uniaxialMaterial Concrete01 [expr $IDconcCoreBond + $ctr]  [expr -$fpcc] [expr -$epcc1] [expr -0.8*$fpcc] [expr -$eccu_C];"
                    puts $matDataFileID "uniaxialMaterial Concrete01 [expr $IDconcCoverBond + $ctr] [expr -$fpc] [expr -$epc] [expr -0.8*$fpc] [expr -$ecu_UC];"
                }
            }
		}
		
		# For column section aggregation
		uniaxialMaterial Elastic [expr $IDShear + $ctr]   [expr (9./10.)*$Gcc*[lindex $Acol 1]];		# 0.9A is the effective shear area of circular cross section
		uniaxialMaterial Elastic [expr $IDTorsion + $ctr] [expr 1.0*$Gcc*[lindex $Jcol 1]];        	# Based on SDC-1.6 sec. 5.6.2 and sec. 1.1,reduction is NOT required 
		
        if {$write_model_files == 1} {
            puts $matDataFileID "uniaxialMaterial Elastic [expr $IDShear + $ctr]   [expr (9./10.)*$Gcc*[lindex $Acol 1]];"
            puts $matDataFileID "uniaxialMaterial Elastic [expr $IDTorsion + $ctr] [expr 1.0*$Gcc*[lindex $Jcol 1]];"
        }
        
		set ctr [expr $ctr + 1];
		
		# if j == 1 or j == num_secdef_per_col, only then it's a possible plastic hinge
		if {$write_model_files == 1} {
			if {$j == 1} {
				puts $predFileID1 "$rhoTransAll";
				puts $predFileID1 "$fyH";
				puts $predFileID1 "$Es";
				puts $predFileID1 "$fpc";
				puts $predFileID1 "[expr (($pi*$DcolAll**2)/4)]";
				puts $predFileID1 "$numBarCol1";
				close $predFileID1;
			}
			if {$j == $num_secdef_per_col} {
				puts $predFileID2 "$rhoTransAll";
				puts $predFileID2 "$fyH";
				puts $predFileID2 "$Es";
				puts $predFileID2 "$fpc";
				puts $predFileID2 "[expr (($pi*$DcolAll**2)/4)]";
				puts $predFileID2 "$numBarCol1";
				close $predFileID2;
			}
		}
	}
    
    # For BH section aggregation
	# each col has one base hinge, so use correct ctr
	if {$existsBH == 1} {
		uniaxialMaterial Elastic [expr $IDShearBH + ($i - 1)*$num_secdef_per_col + 1]     [expr (9./10.)*$Gcc1*[lindex $ABH 1]];
		uniaxialMaterial Elastic [expr $IDTorsionBH + ($i - 1)*$num_secdef_per_col + 1]   [expr $Gcc1*[lindex $JBH 1]];
        if {$write_model_files == 1} {
            puts $matDataFileID "uniaxialMaterial Elastic [expr $IDShearBH + ($i - 1)*$num_secdef_per_col + 1]     [expr (9./10.)*$Gcc1*[lindex $ABH 1]];"
            puts $matDataFileID "uniaxialMaterial Elastic [expr $IDTorsionBH + ($i - 1)*$num_secdef_per_col + 1]   [expr $Gcc1*[lindex $JBH 1]];"
        }
    }
	
}
#
#
# ---------------------------------------------REINFORCING STEEL Parameters--------------------------------------------
#
# Nominal Steel Properties Derived from SDC ver. 1.4 
#
set ctr 1;
set fy_e 0;
for {set i 1} {$i <= [expr $nCols * $nBent]} {incr i 1} {
	for {set j 1} {$j <= [expr $num_secdef_per_col]} {incr j 1} {
		set ToverY_max 0.;
		for {set k 1} {$k <= [expr $num_bar_clusters]} {incr k 1} {
		
			set varName       "fy";
			set fy           [expr $[append varName "_col_$i\_secdef_$j\_barcluster_$k"]*$ksi];
			set fy_e         [expr $fy_e + $fy]
            
			set varName       "fu";
			set fu           [expr $[append varName "_col_$i\_secdef_$j\_barcluster_$k"]*$ksi];
			
			set ToverY 		 [expr $fu/$fy];
			
			if {$ToverY > $ToverY_max} {
				set ToverY_max $ToverY;
			}
			
			set varName       "b";
			set b            [expr $[append varName "_col_$i\_secdef_$j\_barcluster_$k"]*1.0];
			
			set varName       "Es";
			set Es           [expr $[append varName "_col_$i\_secdef_$j\_barcluster_$k"]*$ksi];
			
			set Us           0.2;                                       # Poisson's ratio
			set Gs           [expr $Es/(2*(1+$Us))];                    # Shear Modulus of Elasticity
			
			set Fy $fy;							        # STEEL yield stress
			
			set bs           $b;
			
			set R0 20.;									# control the transition from elastic to plastic branches
			set cR1 0.925;								# control the transition from elastic to plastic branches
			set cR2 0.15;								# control the transition from elastic to plastic branches
			
			# uniaxialMaterial Steel02 $IDSteel $Fy $Es $bs $R0 $cR1 $cR2;
			if {$steelMaterial == "SteelMPF"} {
                uniaxialMaterial SteelMPF [expr $IDSteel + $ctr] $Fy $Fy $Es $bs $bs $R0 $cR1 $cR2;
            } elseif {$steelMaterial == "Steel01"} {
                uniaxialMaterial Steel01 [expr $IDSteel + $ctr] $Fy $Es $bs 
            }
            if {$bond == 1} {
                if {$j == 1 || $j == $num_secdef_per_col} {
                    set fy_bond $fy
                    set fu_bond $fu
                    set tempVarName "fpc";
                    set fpc_bond [expr $[append tempVarName "_col_$i\_secdef_$j"]*$ksi]; 
                    set sy_bond [expr (0.1 * (($DbarSingle)/4000. * ($fy_bond/$psi)/(($fpc_bond/$psi)**0.5) * (2. * 0.4 + 1)) ** (1/0.4)) + 0.013]
                    set su_bond [expr 30.*$sy_bond];
                    set b_bond 0.4;
                    set R_bond 0.7;
                    uniaxialMaterial Bond_SP01 [expr $IDSteelBond + $ctr] $fy_bond $sy_bond $fu_bond $su_bond $b_bond $R_bond
                }
            }
            if {$write_model_files == 1} {
				if {$steelMaterial == "SteelMPF"} {
                    puts $matDataFileID "uniaxialMaterial SteelMPF [expr $IDSteel + $ctr] $Fy $Fy $Es $bs $bs $R0 $cR1 $cR2;"
				} elseif {$steelMaterial == "Steel01"} {
                    puts $matDataFileID "uniaxialMaterial Steel01 [expr $IDSteel + $ctr] $Fy $Es $bs;"
                }
                puts $colRebarMatTagInfoFileID "[expr $IDSteel + $ctr]";
                if {$bond == 1} {
                    if {$j == 1 || $j == $num_secdef_per_col} {
                        puts $matDataFileID "uniaxialMaterial Bond_SP01 [expr $IDSteelBond + $ctr] $fy_bond $sy_bond $fu_bond $su_bond $b_bond $R_bond"
                    }
                }
			}
            
			set ctr [expr $ctr + 1];
		}
		
		# if j == 1 or j == num_secdef_per_col, only then it's a possible plastic hinge
		if {$write_model_files == 1} {
			if {$j == 1} {
				set predFileID1 [open "$model_info_directory/predictor_info_col_rebar_strain_damage_col_$i\_edge_1.txt" "a"];
				puts $predFileID1 "$ToverY_max";
				close $predFileID1;
			}
			if {$j == $num_secdef_per_col} {
				set predFileID2 [open "$model_info_directory/predictor_info_col_rebar_strain_damage_col_$i\_edge_2.txt" "a"];
				puts $predFileID2 "$ToverY_max";
				close $predFileID2;
			}
		}
	}
}


# SSI
if {$ssi == 1} {
    uniaxialMaterial Elastic [expr $ID_SSI + 1] $kx_ssi
    uniaxialMaterial Elastic [expr $ID_SSI + 2] $ky_ssi
    uniaxialMaterial Elastic [expr $ID_SSI + 3] $kz_ssi
    uniaxialMaterial Elastic [expr $ID_SSI + 4] $mx_ssi
    uniaxialMaterial Elastic [expr $ID_SSI + 5] $my_ssi
    uniaxialMaterial Elastic [expr $ID_SSI + 6] $mz_ssi
    
    uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 1] $cx_ssi  1.
    uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 2] $cy_ssi  1.
    uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 3] $cz_ssi  1.
    uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 4] $cmx_ssi 1.
    uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 5] $cmy_ssi 1.
    uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 6] $cmz_ssi 1.
    
    if {$write_model_files == 1} {
        puts $matDataFileID "uniaxialMaterial Elastic [expr $ID_SSI + 1] $kx_ssi"
        puts $matDataFileID "uniaxialMaterial Elastic [expr $ID_SSI + 2] $ky_ssi"
        puts $matDataFileID "uniaxialMaterial Elastic [expr $ID_SSI + 3] $kz_ssi"
        puts $matDataFileID "uniaxialMaterial Elastic [expr $ID_SSI + 4] $mx_ssi"
        puts $matDataFileID "uniaxialMaterial Elastic [expr $ID_SSI + 5] $my_ssi"
        puts $matDataFileID "uniaxialMaterial Elastic [expr $ID_SSI + 6] $mz_ssi"
        
        puts $matDataFileID "uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 1] $cx_ssi  1."
        puts $matDataFileID "uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 2] $cy_ssi  1."
        puts $matDataFileID "uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 3] $cz_ssi  1."
        puts $matDataFileID "uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 4] $cmx_ssi 1."
        puts $matDataFileID "uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 5] $cmy_ssi 1."
        puts $matDataFileID "uniaxialMaterial Viscous [expr ($ID_SSI * 10) + 6] $cmz_ssi 1."  
    }
}

set fy_e [expr $fy_e/($ctr - 1)]

if {$write_model_files == 1} {
close $matDataFileID
close $colRebarMatTagInfoFileID
}

#####################################################################################################################################################################
########################################################################### COLUMN SECTION ##########################################################################
#####################################################################################################################################################################

# set maxFiber [expr 0.04*$m];
# set minFiber [expr 0.03*$m];

# set maxFiber [expr 0.1*$m];
# set minFiber [expr 0.08*$m];

set maxFiber [expr 0.8*$m];
set minFiber [expr 0.05*$m];

source "$model_files_path/createCirclePatches.tcl";

# COLUMN SECTION ---------------------------------------------------------------------------------------
set ColSecTag 2000;
set BondSecTag 4000;

set yCenter 			0.;																										# Y - center of section
set zCenter 			0.;																										# Z - center of section
						
set intRadCore			[expr 0.00];																    						# Internal radius of core concrete
set extRadCore			[expr [lindex $Dcol 1]/2. - [lindex $cover 1] - [lindex $DbarH 1]];										# External radius of core concrete
						
set intRadCover			[expr $extRadCore];																						# Internal radius of the cover concrete
set extRadCover			[expr [lindex $Dcol 1]/2.];																				# External radius of the cover concrete
						
set numBar				[expr [lindex $numBarCol 1]];																			# Number of reinforcing bars along layer
set areaBar				[expr [lindex $Abar 1]];																				# Area of individual reinforcing bar (bundles of 2)
set radius				[expr [lindex $Dcol 1]/2. - [lindex $cover 1] - [lindex $DbarH 1] - [lindex $Dbar 1]/2.];				# Radius of reinforcing layer
set theta				[expr 360.0/[lindex $numBarCol 1]];																		# Angle increment between bars

set ctr 1;
set ctrSec 1;
for {set i 1} {$i <= [expr $nCols * $nBent]} {incr i 1} {
	for {set j 1} {$j <= [expr $num_secdef_per_col]} {incr j 1} {
		
        if {$bond == 1} {
            set ctrBond $ctr
        }
        
		if {$write_model_files == 1} {
			set sectionID [expr $ColSecTag + $ctrSec];
			set ColSecDefFileID [open "$model_info_directory/fib_secdef_sec_$sectionID\.txt" "w+"];
			puts $ColSecDefFileID "section Fiber [expr $ColSecTag*10 + $ctrSec] -torsion [expr $IDTorsion + $ctrSec] \{";
		} else {
            set ColSecDefFileID -1;
        }
		
		section Fiber [expr $ColSecTag*10 + $ctrSec] -torsion [expr $IDTorsion + $ctrSec] {
			# createCirclePatches $matTag $yCenter $zCenter $radiusInner $radiusOuter \
				$minFiber $maxFiber $startAng $endAng ?$exportID?
			createCirclePatches [expr $IDconcCover + $ctrSec] $yCenter $zCenter $intRadCover $extRadCover \
				$minFiber $maxFiber 0.0 360.0 $ColSecDefFileID;
			createCirclePatches [expr $IDconcCore + $ctrSec] $yCenter $zCenter $intRadCore $extRadCore \
				$minFiber $maxFiber 0.0 360.0 $ColSecDefFileID;
			
			for {set k 1} {$k <= [expr $num_bar_clusters]} {incr k 1} {
				set startAng [expr ($k - 1)*((int($numBar/$num_bar_clusters) - 1)*$theta + $theta)];
				set endAng [expr $startAng + (int($numBar/$num_bar_clusters) - 1)*$theta];
				layer circ [expr $IDSteel + $ctr] [expr int($numBar/$num_bar_clusters)] $areaBar $yCenter $zCenter $radius $startAng $endAng; 
				if {$write_model_files == 1} {
					puts $ColSecDefFileID "layer circ [expr $IDSteel + $ctr] [expr int($numBar/$num_bar_clusters)] $areaBar $yCenter $zCenter $radius $startAng $endAng;"
				}
				set ctr [expr $ctr + 1];
			}
			
		}
		
		section Aggregator [expr $ColSecTag + $ctrSec] [expr $IDShear + $ctrSec] Vy [expr $IDShear + $ctrSec] Vz -section [expr $ColSecTag*10 + $ctrSec];
		
		if {$write_model_files == 1} {
			puts $ColSecDefFileID "\}";
			puts $ColSecDefFileID "section Aggregator [expr $ColSecTag + $ctrSec] [expr $IDShear + $ctrSec] Vy [expr $IDShear + $ctrSec] Vz -section [expr $ColSecTag*10 + $ctrSec];";
			close $ColSecDefFileID;
		}
		
        if {$bond == 1} {
            if {$j == 1 && $existsBH == 1} {
                # do nothing
            } elseif {($j == 1 && $existsBH == 0) || $j == $num_secdef_per_col} {
                if {$write_model_files == 1} {
                    set sectionID [expr $BondSecTag + $ctrSec];
                    set BondSecDefFileID [open "$model_info_directory/fib_secdef_sec_$sectionID\.txt" "w+"];
                    puts $BondSecDefFileID "section Fiber [expr $BondSecTag*10 + $ctrSec] -torsion [expr $ID_Rigid] \{";
                } else {
                    set BondSecDefFileID -1
                }
                
                section Fiber [expr $BondSecTag*10 + $ctrSec] -torsion [expr $ID_Rigid] {
                    # createCirclePatches $matTag $yCenter $zCenter $radiusInner $radiusOuter \
                        $minFiber $maxFiber $startAng $endAng ?$exportID?
                    createCirclePatches [expr $IDconcCoverBond + $ctrSec] $yCenter $zCenter $intRadCover $extRadCover \
                        $minFiber $maxFiber 0.0 360.0 $BondSecDefFileID;
                    createCirclePatches [expr $IDconcCoreBond + $ctrSec] $yCenter $zCenter $intRadCore $extRadCore \
                        $minFiber $maxFiber 0.0 360.0 $BondSecDefFileID;
                    
                    for {set k 1} {$k <= [expr $num_bar_clusters]} {incr k 1} {
                        set startAng [expr ($k - 1)*((int($numBar/$num_bar_clusters) - 1)*$theta + $theta)];
                        set endAng [expr $startAng + (int($numBar/$num_bar_clusters) - 1)*$theta];
                        layer circ [expr $IDSteelBond + $ctrBond] [expr int($numBar/$num_bar_clusters)] $areaBar $yCenter $zCenter $radius $startAng $endAng; 
                        if {$write_model_files == 1} {
                            puts $BondSecDefFileID "layer circ [expr $IDSteelBond + $ctrBond] [expr int($numBar/$num_bar_clusters)] $areaBar $yCenter $zCenter $radius $startAng $endAng;"
                        }
                        set ctrBond [expr $ctrBond + 1];
                    }
                    
                }
                
                section Aggregator [expr $BondSecTag + $ctrSec] [expr $ID_Rigid] Vy [expr $ID_Rigid] Vz -section [expr $BondSecTag*10 + $ctrSec];
                
                if {$write_model_files == 1} {
                    puts $BondSecDefFileID "\}";
                    puts $BondSecDefFileID "section Aggregator [expr $BondSecTag + $ctrSec] [expr $ID_Rigid] Vy [expr $ID_Rigid] Vz -section [expr $BondSecTag*10 + $ctrSec];"
                    close $BondSecDefFileID;
                }
            }
        }
		set ctrSec [expr $ctrSec + 1];
	
    }
}	

# BASE HINGE --------------------------------------------------------------------------------------
set BHSecTag 3000;
if {$existsBH == 1} {
	
	set yCenter 			0.;																				# Y - center of section
	set zCenter 			0.;																				# Z - center of section
	set intRadCore			[expr 0.00*$in];																# Internal radius of core concrete
	set extRadCore			[expr [lindex $DBH 1]/2.-[lindex $coverBH 1]];									# External radius of core concrete
	
	set intRadCover			[expr $extRadCore];																# Internal radius of the cover concrete
	set extRadCover			[expr [lindex $DBH 1]/2.];														# External radius of the cover concrete
	
	set numBar				[expr [lindex $numBarBH 1]];													# Number of reinforcing bars along layer
	set areaBar				[expr ($pi*[lindex $DbarBH 1]**2/4.)*([lindex $bundlesBH 1])];					# Area of individual reinforcing bar (bundles of 2)
	set radius				[expr [lindex $DBH 1]/2. - [lindex $coverBH 1] - [lindex $DbarBH 1]/2.];		# Radius of reinforcing layer
	set theta				[expr 360.0/[lindex $numBarBH 1]];												# Angle increment between bars
	
	set ctrSec 1;
	for {set i 1} {$i <= [expr $nCols * $nBent]} {incr i 1} {
		if {$write_model_files == 1} {
			# naming section id per col
			set sectionID [expr $BHSecTag + $ctrSec];
			set BHSecDefFileID [open "$model_info_directory/fib_secdef_sec_$sectionID\.txt" "w+"];
			puts $BHSecDefFileID "section Fiber [expr $BHSecTag*10 + $ctrSec] -torsion [expr $IDTorsionBH + $ctrSec] \{"
		} else {
            set BHSecDefFileID -1;
        }
		
		section Fiber [expr $BHSecTag*10 + $ctrSec] -torsion [expr $IDTorsionBH + $ctrSec] {
			# createCirclePatches $matTag $yCenter $zCenter $radiusInner $radiusOuter \
				$minFiber $maxFiber $startAng $endAng ?$exportID?
			createCirclePatches [expr $IDconcCore + $ctrSec] $yCenter $zCenter $intRadCover $extRadCover \
				$minFiber $maxFiber 0.0 360.0 $BHSecDefFileID;
			createCirclePatches [expr $IDconcCore + $ctrSec] $yCenter $zCenter $intRadCore $extRadCore \
				$minFiber $maxFiber 0.0 360.0 $BHSecDefFileID;
			
			set randCtr [expr int(rand()*($num_bar_clusters*$num_secdef_per_col - 1 + 1)) + 1]
			set randCtr 1
			layer circ [expr $IDSteel + ($i - 1)*$num_bar_clusters*$num_secdef_per_col + $randCtr] $numBar $areaBar $yCenter $zCenter $radius 0. [expr 360. - $theta]
			if {$write_model_files == 1} {
				puts $BHSecDefFileID "layer circ [expr $IDSteel + ($i - 1)*$num_bar_clusters*$num_secdef_per_col + $randCtr] $numBar $areaBar $yCenter $zCenter $radius 0. [expr 360. - $theta]"
			}
		}
		
		section Aggregator [expr $BHSecTag + $ctrSec] [expr $IDShearBH + $ctrSec] Vy [expr $IDShearBH + $ctrSec] Vz -section [expr $BHSecTag*10 + $ctrSec];
		
		if {$write_model_files == 1} {
			puts $BHSecDefFileID "\}";
			puts $BHSecDefFileID "section Aggregator [expr $BHSecTag + $ctrSec] [expr $IDShearBH + $ctrSec] Vy [expr $IDShearBH + $ctrSec] Vz -section [expr $BHSecTag*10 + $ctrSec];";
			close $BHSecDefFileID
		}
		
        if {$bond == 1} {
            if {$write_model_files == 1} {
                set sectionID [expr $BondSecTag + $ctrSec];
                set BondSecDefFileID [open "$model_info_directory/fib_secdef_sec_$sectionID\.txt" "w+"];
                puts $BondSecDefFileID "section Fiber [expr $BondSecTag*10 + $ctrSec] -torsion [expr $ID_Rigid] \{";
            } else {
                set BondSecDefFileID -1;
            }
            
            section Fiber [expr $BondSecTag*10 + $ctrSec] -torsion [expr $ID_Rigid] {
                # createCirclePatches $matTag $yCenter $zCenter $radiusInner $radiusOuter \
                    $minFiber $maxFiber $startAng $endAng ?$exportID?
                createCirclePatches [expr $IDconcCoreBond + $ctrSec] $yCenter $zCenter $intRadCover $extRadCover \
                    $minFiber $maxFiber 0.0 360.0 $BondSecDefFileID;
                createCirclePatches [expr $IDconcCoreBond + $ctrSec] $yCenter $zCenter $intRadCore $extRadCore \
                    $minFiber $maxFiber 0.0 360.0 $BondSecDefFileID;
                
                layer circ [expr $IDSteelBond + ($i - 1)*$num_bar_clusters*$num_secdef_per_col + $randCtr] $numBar $areaBar $yCenter $zCenter $radius 0. [expr 360. - $theta]
                if {$write_model_files == 1} {
                    puts $BondSecDefFileID "layer circ [expr $IDSteelBond + ($i - 1)*$num_bar_clusters*$num_secdef_per_col + $randCtr] $numBar $areaBar $yCenter $zCenter $radius 0. [expr 360. - $theta]"
                }
            }
            
            section Aggregator [expr $BondSecTag + $ctrSec] [expr $ID_Rigid] Vy [expr $ID_Rigid] Vz -section [expr $BondSecTag*10 + $ctrSec];
            
            if {$write_model_files == 1} {
                puts $BondSecDefFileID "\}";
                puts $BondSecDefFileID "section Aggregator [expr $BondSecTag + $ctrSec] [expr $ID_Rigid] Vy [expr $ID_Rigid] Vz -section [expr $BondSecTag*10 + $ctrSec];";
                close $BondSecDefFileID;
            }
        }
        
		# jump num_secdef_per_col
		set ctrSec [expr $ctrSec + $num_secdef_per_col];
	}
}