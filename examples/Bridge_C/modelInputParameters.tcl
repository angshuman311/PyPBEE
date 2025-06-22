# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

set skew 36.
set skewr [expr $skew*$pi/180.0]

set vert             $vertical_gm;
set ssi              $ssi_springs;
set existsBH         $base_hinge;
set steelMaterial    $steel_material
set colEleType		 $col_elem_type;								# (dispBeamColumn, forceBeamColumn)
set bf_skew          $backfill_resp_skew

if {$col_elem_type == "1"} {
    set colEleType "forceBeamColumn";
    set colEleDistPattern  {1 0 0};
    set colEleDistHgts  {"\[expr \$colHgtCurrentBent\]" "\[expr 0.\]" "\[expr 0.\]"};
    set nEleCol [expr [lindex $colEleDistPattern 0] + [lindex $colEleDistPattern 1] + [lindex $colEleDistPattern 2]];
    set nIntPts		7;
    set intType	"Lobatto";
    set bond    0;
}

if {$col_elem_type == "2"} {
    set colEleType "forceBeamColumn";
    set colEleDistPattern  {1 0 0};
    set colEleDistHgts  {"\[expr \$colHgtCurrentBent\]" "\[expr 0.\]" "\[expr 0.\]"};
    set nEleCol [expr [lindex $colEleDistPattern 0] + [lindex $colEleDistPattern 1] + [lindex $colEleDistPattern 2]];
    set nIntPts		7;
    set intType	"Lobatto";
    set bond    1;
}

if {$col_elem_type == "3"} {
    set colEleType "dispBeamColumn";
    set colEleDistPattern  {0 1 4};
    set colEleDistHgts  {"\[expr 0.\]" "\[expr \$colHgtCurrentBent - (\$Lp1 + \$Lp2)\]" "\[expr \$Lp1 + \$Lp2\]"};
    set nEleCol [expr [lindex $colEleDistPattern 0] + [lindex $colEleDistPattern 1] + [lindex $colEleDistPattern 2]];
    set nIntPts		3;
    set intType	"Lobatto";
    set bond    1;
}

if {$col_elem_type == "4"} {
    set colEleType "forceBeamColumn";
    set colEleDistPattern  {0 1 1};
    set colEleDistHgts  {"\[expr 0.\]" "\[expr \$colHgtCurrentBent - 6.*(\$Lp1 + \$Lp2)\]" "\[expr 6.*(\$Lp1 + \$Lp2)\]"};
    set nEleCol [expr [lindex $colEleDistPattern 0] + [lindex $colEleDistPattern 1] + [lindex $colEleDistPattern 2]];
    set nIntPts		3;
    set intType	"Lobatto";
    set bond    0;
}

if {$col_elem_type == "5"} {
    set colEleType "forceBeamColumn";
    set colEleDistPattern  {0 2 1};
    set colEleDistHgts  {"\[expr 0.\]" "\[expr \$colHgtCurrentBent - 2.*(\$Lp1)\]" "\[expr 2.*(\$Lp1)\]"};
    set nEleCol [expr [lindex $colEleDistPattern 0] + [lindex $colEleDistPattern 1] + [lindex $colEleDistPattern 2]];
    set nIntPts		2;
    set intType	"Lobatto";
    set bond    1;
}

# Deck
# ----------------------
set L1				[expr 156.12*$ft];								# Length of span 1 
set L2				[expr 144.*$ft];								# Length of span 2
set L3				[expr 118.08*$ft];								# Length of span 3
set L				{spanLengths $L1 $L2 $L3};							# list of lengths of spans
set L_tot			[expr $L1 + $L2 + $L3];								# Total length
set nSpans			[expr [llength $L] - 1];						# number of bridge spans
set nBent			[expr $nSpans - 1];								# number of column bents
set nEleSpan		10;												# number of elements per span of bridge

set dw				[expr 77.*$ft];									# Deck Width
set dd				[expr 6.3*$ft];									# Deck Depth
set Adeck			[expr 12230811*pow($mm,2)];							# Area of deck
set Jdeck			[expr 2.212e13*pow($mm,4)];								# Torsional constant of deck
set Iydeck			[expr 6.809e12*pow($mm,4)];							# Iy - local
set Izdeck			[expr 5.158e14*pow($mm,4)];							# Iz - local
set gapL			[expr 1.0*$in];												# 1" Gap between the backwall and deck

# Cap Beam
# ----------------------
set nCapBeamEle 1;
set Lcap1L			[expr 32.*$ft];									# Cap Beam Length
set Lcap1R			[expr 32.*$ft];									# Cap Beam Length
set Lcap2L			[expr 34.*$ft];									# Cap Beam Length
set Lcap2R			[expr 34.*$ft];									# Cap Beam Length

set bcap1  			[expr 7.48*$ft];										# Cap Beam Breadth
set hcap1  			[expr 6.3*$ft];										# Cap Beam Height
set Acap1  			[expr $bcap1*$hcap1];									# Area of cap beam
set Jcap1  			[expr ($bcap1*$hcap1/12)*($bcap1**2+$hcap1**2)];		# Torsional constant of cap beam
set Iycap1 			[expr ($bcap1*($hcap1**3))/12];						# Iy - local
set Izcap1 			[expr ($hcap1*($bcap1**3))/12];						# Iz - local

set bcap2  			[expr 7.48*$ft];										# Cap Beam Breadth
set hcap2  			[expr 6.3*$ft];										# Cap Beam Height
set Acap2  			[expr $bcap2*$hcap2];									# Area of cap beam
set Jcap2  			[expr ($bcap2*$hcap2/12)*($bcap2**2+$hcap2**2)];		# Torsional constant of cap beam
set Iycap2 			[expr ($bcap2*($hcap2**3))/12];						# Iy - local
set Izcap2 			[expr ($hcap2*($bcap2**3))/12];						# Iz - local


# Lists of above parameters bent wise
set Lcap			{CapBeamLengths $Lcap1R $Lcap1L $Lcap2R $Lcap2L};					# Cap Beam Lengths List {Bent(Right, Left), ..}
set Acap			{CapBeamAreas $Acap1 $Acap2};
set Jcap			{CapBeamJ $Jcap1 $Jcap2};		
set Iycap			{CapBeamIy $Iycap1 $Iycap2};		
set Izcap			{CapBeamIz $Izcap1 $Izcap2};		
set bcap 			{CapBeamBreadths $bcap1 $bcap2};
set hcap 			{CapBeamHeights $hcap1 $hcap2};

# Columns
# ----------------------
set nCols			3;												# number of columns per bent
set dcg				[expr 3.23*$ft];								# Deck Centeriod-Column Top Dist.

set DcolAll                 [expr $all_col_dia_in_ft*$ft];
set all_rho_long              [expr $all_rho_long*1.0];
set num_bar_clusters          [expr int($num_bar_clusters*1.0)];
set num_secdef_per_col			[expr int($num_secdef_per_col*1.0)];												# number of section definitions per entire column (NOT the number of integration points in a column element)

set rhoTransAll		    [expr 0.5*$all_rho_long];													# transverse reinforcement ratio
# set rhoTransAll		    [expr 0.75*$all_rho_long];													# transverse reinforcement ratio
# set rhoTransAll		    [expr 0.01];													# transverse reinforcement ratio
set Dcol1				[expr $DcolAll];												# Column Diameter
set Acol1				[expr (($pi*$Dcol1**2)/4)];										# Area of column
set Jcol1				[expr ($pi*($Dcol1/2)**4)/2];									# Torsional constant of column
set Izcol1				[expr ($pi*($Dcol1/2)**4)/4];									# Iz - local
set Iycol1				[expr ($pi*($Dcol1/2)**4)/4];									# Iy - local
set rhoLong1			[expr $all_rho_long];                 							# long. reinf. ratio
set rhoTrans1			[expr $rhoTransAll];                 							# trans. reinf. ratio
set numBarCol1			24;																# Number of rebar per column
set Dbar1				[expr sqrt(4.*$rhoLong1*$Acol1/($numBarCol1*$pi))];				# Long. bar dia
set Abar1				[expr ($pi*$Dbar1**2)/4.];										# Long. bar area
set cover1				[expr 2.*$in];													# concrete cover
set sH1					[expr 85*$mm];												# spacing of hoops
set h2hOuter1			[expr $Dcol1 - 2.*$cover1];
set AbarH1				[expr $rhoTrans1*$h2hOuter1*$sH1/4.];
set DbarH1				[expr sqrt(4.*$AbarH1/$pi)];
set sHprime1			[expr $sH1 - $DbarH1];
set ds1					[expr $h2hOuter1 - $DbarH1];
set Ac1					[expr ($pi*$ds1**2)/4.];
set rho_cc1				[expr $numBarCol1*$Abar1/$Ac1];
set A_cc1				[expr $Ac1*(1. - $rho_cc1)];

set Dcol			{ColumnDiameters $Dcol1 $Dcol1};
set Acol			{ColumnAreas $Acol1 $Acol1};
set Jcol			{ColumnJ $Jcol1 $Jcol1};
set Izcol		    {ColumnIz $Izcol1 $Izcol1};
set Iycol		    {ColumnIy $Iycol1 $Iycol1};
set rhoLong			{ColumnReinfRatioLong $rhoLong1 $rhoLong1};
set rhoTrans		{ColumnReinfRatioTrans $rhoTrans1 $rhoTrans1};
set numBarCol		{ColumnReinfNumberofBars $numBarCol1 $numBarCol1};
set Dbar			{ColumnReinfDia $Dbar1 $Dbar1};
set Abar			{ColumnReinfArea $Abar1 $Abar1};
set cover			{ColumnCover $cover1 $cover1};
set sH			    {ColumnHoopSpacingC2C $sH1 $sH1};
set h2hOuter	    {ColumnHoopToHoopSpacingOuter $h2hOuter1 $h2hOuter1}
set AbarH		    {ColumnHoopArea $AbarH1 $AbarH1};
set DbarH		    {ColumnHoopDia $DbarH1 $DbarH1};
set sHprime		    {ColumnHoopSpacingClear $sHprime1 $sHprime1};
set ds			    {ColumnSpiralDiaC2C $ds1 $ds1};
set Ac			    {ColumnAreaCore $Ac1 $Ac1};
set rho_cc		    {ColumnReinfRatioLongCore $rho_cc1 $rho_cc1};
set A_cc		    {ColumnAreaCoreConc $A_cc1 $A_cc1};

if {$all_col_dia_in_ft < 6.} {
set DbarSingle		[expr 1.410*$in];
} else {
set DbarSingle		[expr 1.693*$in];
}

set HColumnMax		[expr 24.6*$ft];								# Column Height
set HCol			{ColHeights $HColumnMax $HColumnMax};			# list of column heights in each bent
                                
set maxIters		50;
set tol				1e-9;

set DBH1			[expr 0.6*$Dcol1];								# Column BH Diameter
set DbarBH1			[expr 1.270*$in];                 				# using #35 reinforcing bars
set numBarBH1		8;												# Number of rebar per column BH
set coverBH1		[expr $DBH1/4.];
set ABH1			[expr (($pi*$DBH1**2)/4)];						# Area of column BH
set JBH1			[expr ($pi*($DBH1/2)**4)/2];					# Torsional constant of column BH
set IzBH1			[expr ($pi*($DBH1/2)**4)/4];					# Iz - local BH
set IyBH1			[expr ($pi*($DBH1/2)**4)/4];					# Iy - local BH
set bundleBH1		1;												# Bundles of reinforcing bar to be placed at each fiber location

# Lists of above parameters bent wise
set DBH				{BaseHingeDiameters $DBH1 $DBH1};						# Column BH Diameter
set DbarBH			{BaseHingeReinfDia $DbarBH1 $DbarBH1};					# using #35 reinforcing bars
set numBarBH		{BaseHingeReinfNumberofBars $numBarBH1 $numBarBH1};		# Number of rebar per column BH
set coverBH			{BaseHingeCover $coverBH1 $coverBH1};
set ABH				{BaseHingeAreas $ABH1 $ABH1};							# Area of column BH
set JBH				{BaseHingeJ $JBH1 $JBH1};								# Torsional constant of column BH
set IzBH			{BaseHingeIz $IzBH1 $IzBH1};							# Iz - local BH
set IyBH			{BaseHingeIy $IyBH1 $IyBH1};							# Iy - local BH
set bundlesBH		{BaseHingeReinfBundle $bundleBH1 $bundleBH1};               # Bundles of reinforcing bar to be placed at each fiber location

set BaseHingeHeight1 [expr 25.*$mm];
set BaseHingeHeight {BaseHingeHeights $BaseHingeHeight1 $BaseHingeHeight1};

set BHEleType "dispBeamColumn";

# SSI
# ----------------------
set kx_ssi [expr 660. * $MN / $m]
set ky_ssi [expr 698. * $MN / $m]
set kz_ssi [expr 545. * $MN / $m]
set mx_ssi [expr 3065.* $MN * $m / 1.]
set my_ssi [expr 5322.* $MN * $m / 1.]
set mz_ssi [expr 6920.* $MN * $m / 1.]

set cx_ssi [expr 13. * $MN * $sec / $m]
set cy_ssi [expr 13. * $MN * $sec / $m]
set cz_ssi [expr 14. * $MN * $sec / $m]
set cmx_ssi [expr 21. * $MN * $m * $sec / 1.]
set cmy_ssi [expr 39. * $MN * $m * $sec / 1.]
set cmz_ssi [expr 42. * $MN * $m * $sec / 1.]

# Bearings
# ----------------------
set nBearings		9;												# number of bearing pads
set hb [expr 115.*$mm];
set hb [expr $hb/1.5];
set Gb [expr 0.69*$MPa];
set Ab [expr 460.*$mm*460.*$mm];
set kInit1 [expr ($nBearings/5.)*$Gb*$Ab/$hb];
set qd1 [expr ($nBearings/5.)*$Gb*$Ab/2.];

# Shear Keys
# --------------------------------------------------------------
set b_SK [expr 46.12*$in];
set d_SK [expr 64.9*$in];
set d1_SK [expr 65.94*$in];
set d_out_cracked_SK [expr 32.84*$in];
set h_SK [expr 92.5*$in];
set sv_SK [expr 11.82*$in];
set sh_SK [expr 17.73*$in];
set sh_SK_back [expr 150.*$mm];
set s_SK [expr max($sh_SK,$sv_SK)];
set a_SK [expr 8.865*$in];
set beta_SK 27.;


set db_22 [expr 22.225*$mm];
set Ab_22 [expr $pi*pow($db_22,2)/4.];

set db_16 [expr 15.875*$mm];
set Ab_16 [expr $pi*pow($db_16,2)/4.];

set db_36 [expr 35.81*$mm];
set Ab_36 [expr $pi*pow($db_36,2)/4.];

set db_19 [expr 19.05*$mm];
set Ab_19 [expr $pi*pow($db_19,2)/4.];

# Abutment modeling
# ----------------------
set abutvar			0.3;											# maximum percentage of strength/stiffness increase for skew angle of 60 degrees

# Derived quantities for skewed bridges
# ---------------------------------------

set dxc {};
for {set i 1} {$i <= $nBent} {incr i 1} {
lappend dxc	[expr ([lindex $Lcap [expr 2*$i - 1]])*tan($skew*$pi/180.0)];			# Cap Beam Longitudinal Offset Regarding Skew
lappend dxc	[expr ([lindex $Lcap [expr 2*$i]])*tan($skew*$pi/180.0)];				# Cap Beam Longitudinal Offset Regarding Skew
}
set dxc [linsert $dxc 0 "CapBeamLongitudinalOffsets"];

set dx				[expr ($dw/2)*tan($skew*$pi/180.0)];			# Abutment Longitudinal Offset Regarding Skew


# Geometric Transformation
# ------------------------------
set gT				"Linear";
# set gT				"PDelta";

# Intergration
# ----------------------
set BHintPts 2;
set BHintType "Lobatto";

# Analysis
# ----------------------
set GravityAnalysisDone "No";
