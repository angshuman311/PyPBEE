# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------

set addZerosFor 0.0;
set scaleFac_L 1.0;
set scaleFac_T -1.0;
set scaleFac_V 1.0;
set dtAnalysis 0.009;
set gmSkew 0.0;
set algorithmBasic [split "NewtonLineSearch"];
set testBasic "NormDispIncr";
set showTestBasic 0;
set showTest 0;
set tolDynBasic 1e-04;
set tolDynDisp 1e-04;
set tolDynUnb 1e-01;
set maxNumIterDynBasic 500;
set maxNumIterDyn 2500;
set maxDimKrylov 50;



