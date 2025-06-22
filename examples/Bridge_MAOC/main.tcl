# ------------------------------------------------------------------------------------------------------------------------
# Uncertainty Quantification in Performance Based Earthquake Engineering
# Department of Structural Engineering
# University of California, San Diego
# ------------------------------------------------------------------------------------------------------------------------
# Author: Angshuman Deb, UC San Diego
# ------------------------------------------------------------------------------------------------------------------------

model BasicBuilder -ndm 3 -ndf 6;
puts "Model Setup OK"

source "$model_files_path/units.tcl";

source "$model_files_path/modelInputParameters.tcl";
puts "Model Parameters OK"

if {$write_model_files == 1} {
set ExportID [open "$model_info_directory/model_data.txt" "w"];
} else {
set ExportID -1;
}

source "$model_files_path/materials.tcl";
puts "Materials OK"

source "$model_files_path/geometricTransformation.tcl";
puts "Geometric Transformations OK"

source "$model_files_path/generateDeck.tcl";
puts "Deck Generation OK"

source "$model_files_path/generateColumn.tcl";
puts "Column Generation OK"
puts "Cap Beam Generation OK"

source "$model_files_path/generateAbutment.tcl";
puts "Abutment Generation OK"

source "$model_files_path/generateBC.tcl";
puts "BC Generation OK.tcl"

source "$model_files_path/generateMass.tcl";
puts "Mass Generation OK"

source "$model_files_path/generateGravityLoads.tcl";
puts "Gravity Loads Generation OK"

if {$write_model_files == 1} {
close $ExportID;
}