# PAckage installation
# If you have the package already installed, you can skip this step. If you have made changes to the package and want to use the updated version, you can run this code to reinstall the package from the local path.
using Pkg

Pkg.add(url="/Users/saersamanipour/Desktop/dev/pkgs/jAligner4GCxGC/jAligner4GCxGC.jl")


#################################################################################
################################# EXECUTION #####################################
#################################################################################
using jAligner4GCxGC

# Specify the path to input files
println("Specifying a folder with .txt files for processing. ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
path2files = "test"   # Path to the working folder with files to be aligned
println("Done. There are ", length(filter(f -> endswith(f, ".txt"), readdir(path2files))), " .txt files. ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "\n")



# Convert TXT files to CSV (only needed the first time)
convert_txt2csv(path2files)



# Filter out rows where "Classifications" from ChromaTOF contains 'Bleed', 'Toluene', or 'DCM'
filter_classifications(path2files)



# Filter out rows where "Name" contains unwanted values, e.g., siloxanes & other bleed

filter_names(path2files, keywords)



# Filter out rows where "Name" contains, e.g., compounds with spectra similar to toluene (m/z 91 or 92),
# if their spectra do not contain other m/z values, like molecular ions

filter_names_if_missing_mz(path2files, filters)


# Filter out rows where "Name" contains "Peak" if their spectra contain bleed m/z values

filter_peak_names_with_bleed_mz(path2files, bleed_filters)



# SKIP: no need for the new files/new alignment
# Reset masses classified/filtered by filter_signal_external (previous classification)
# reset_filters(path2files)
# SKIP: no need for the new files/new alignment



# Remove m/z:intensity pairs with intensity values containing '0.' from spectra
remove_zero_intensity_pairs(path2files)



# Remove m/z:intensity pairs with intensity values <= 3% of the highest intensity value
remove_low_intensity_pairs(path2files)



# Classify features by base masses or filter out features by base masses
filtType = "Bleed"   # String with keyword of what is being classified/filtered; will be assigned and stored in 'Origin' column
base = [73, 147, 207, 267, 281, 341, 355, 429, 479, 59, 78, 135, 156, 193, 195, 197, 209, 215, 221, 251, 253, 255, 269, 
    327, 329, 331, 333, 339, 343, 377, 401, 403, 405, 415, 417, 439, 451, 475, 477, 489, 503, 549, 553, 563, 623
]   # Base mass(es) used for scanning 'BaseMass' column (can be a vector [xx, xx, xx])
mz = []   # Other signal(s)/mass(es) that need to be present
int = []    # Relative intensity of those signals with 1 = 100% (if 1 is not provided the base mass intensity is used)
tr1_range = [0, 0]   # 1st dimension retention range mass(es) can occur in
tr2_range = [0, 0]   # 2nd dimension retention range mass(es) can occur in
mz_tol = 0.1    # m/z tolerance
filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)

# Classify features by base masses or filter out features by base masses
filtType = "Toluene"   # String with keyword of what is being classified/filtered; will be assigned and stored in 'Origin' column
base = [91]   # Base mass(es) used for scanning 'BaseMass' column (can be a vector [xx, xx, xx])
mz = [92]   # Other signal(s)/mass(es) that need to be present
int = [0.5]    # Relative intensity of those signals with 1 = 100% (if 1 is not provided the base mass intensity is used)
tr1_range = [0, 0]   # 1st dimension retention range mass(es) can occur in
tr2_range = [0, 0]   # 2nd dimension retention range mass(es) can occur in
mz_tol = 0.1    # m/z tolerance
filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)

# Classify features by base masses or filter out features by base masses
filtType = "DCM"   # String with keyword of what is being classified/filtered; will be assigned and stored in 'Origin' column
base = [49]   # Base mass(es) used for scanning 'BaseMass' column (can be a vector [xx, xx, xx])
mz = [84]   # Other signal(s)/mass(es) that need to be present
int = [0.5]    # Relative intensity of those signals with 1 = 100% (if 1 is not provided the base mass intensity is used)
tr1_range = [0, 0]   # 1st dimension retention range mass(es) can occur in
tr2_range = [0, 0]   # 2nd dimension retention range mass(es) can occur in
mz_tol = 0.1    # m/z tolerance
filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)

# Classify features by base masses or filter out features by base masses
filtType = "Alkane"   # String with keyword of what is being classified/filtered; will be assigned and stored in 'Origin' column
base = [57]   # Base mass(es) used for scanning 'BaseMass' column (can be a vector [xx, xx, xx])
mz = [71, 85]   # Other signal(s)/mass(es) that need to be present
int = [0.3, 0.2]    # Relative intensity of those signals with 1 = 100% (if 1 is not provided the base mass intensity is used)
tr1_range = [0, 0]   # 1st dimension retention range mass(es) can occur in
tr2_range = [0, 0]   # 2nd dimension retention range mass(es) can occur in
mz_tol = 0.1    # m/z tolerance
filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)

# Classify features by base masses or filter out features by base masses
filtType = "Acid"   # String with keyword of what is being classified/filtered; will be assigned and stored in 'Origin' column
base = [55, 60, 73]   # Base mass(es) used for scanning 'BaseMass' column (can be a vector [xx, xx, xx])
mz = [60, 73]    # Other signal(s)/mass(es) that need to be present
int = [0.5, 1]    # Relative intensity of those signals with 1 = 100% (if 1 is not provided the base mass intensity is used)
tr1_range = [0, 0]   # 1st dimension retention range mass(es) can occur in
tr2_range = [0, 0]   # 2nd dimension retention range mass(es) can occur in
mz_tol = 0.1    # m/z tolerance
filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)



# SKIP: merges terpenes and other closely eluting compounds with similar spectra that should not be merged
# # Merge internal drifting feature
# tr1_dev = [3, 4.9]      # Retention deviation right and left of the peak in 1st dimension
# tr2_dev = [0.05, 0.1]   # Retention deviation up and down of the peak in 2nd dimension
# sim_thresh = 0.9      # Spectral similarity threshold for merging
# mz_tol = 0.025        # m/z tolerance for setting up the spectrum for dot product calculations 
# internal_feature_merging(path2files, tr1_dev, tr2_dev, sim_thresh, mz_tol)
# SKIP: merges terpenes etc. that should not be merged



# Align samples
h_thresh = 10000000    # Area above which the highA retention time (RT) threshold is used
rt_tol_lowA = 11 #11     # 1st dimension RT tolerance for low intensity peaks
rt_tol_highA = 11 #11      # 1st dimension RT tolerance for high intensity peaks
rt2_tol_lowA = 0.1      # 2nd dimension RT tolerance for low intensity peaks
rt2_tol_highA = 0.1      # 2nd dimension RT tolerance for high intensity peaks
sim_thresh = 0.8     # Minimum similarity threshold
mz_tol = 0.025      # m/z tolerance for setting up the spectrum for dot product calculations
overviewMerging = "height" # Method how the consensus spectrum is setup; "height" weights the spectrum based on the base peak height, "equal" assumes equal contribution
weights = [5, 1, 1]      # Weights used for score calculation [similarity, library probability, difference RI (library vs. experimental)]
tr1_parts = 5           # Number of parts the 1st RT range is split into
numfrags = 150  # Number of fragments in the spectra to consider for comparison
similarity_method = "DISCO"    # Spectral similarity method: "DISCO" (DIstance & Spectrum Correlation Optimization) or "NDP" (Normalized Dot Product)
feature_align_external(path2files, rt2_tol_lowA, rt2_tol_highA, rt_tol_lowA, h_thresh, rt_tol_highA, sim_thresh, mz_tol, overviewMerging, weights, tr1_parts; numfrags, similarity_method)

#"""

# Filter out rows after alignment, where "Name" contains, e.g., compounds with spectra similar to toluene (m/z 91 or 92),
# if their spectra do not contain other m/z values, like molecular ions
filter_names_mz_postalign(path2files, filters)



# Merge features with the same name and base m/z if their dRI < N
merge_postalign(path2files, 7.0)



# Retrieve InChIKey, SMILES & CAS# values from PubChem
PubChemRetriever(path2files)



# Retrieve chemical class data from ClassyFire Batch
# Specify the path to chromedriver.exe on your PC 
#path2chromedriver = "C:\\Program Files (x86)\\Google\\chromedriver-win64\\chromedriver.exe"
# Call the function
#ClassyFire_Batch(path2files, path2chromedriver)



# Retrieve functional use data from EPA FUse, NORMAN dust trial & in-house harmonized use databases
# Specify paths to file and 3 files with classification data
path_file = "D:\\Projects\\Test\\Aligned_File.csv"
path_harm = "D:\\Projects\\Test\\Databases\\Classes_Harmonized.txt"
path_fuse = "D:\\Projects\\Test\\Databases\\FUse_Database.txt"
path_norm = "D:\\Projects\\Test\\Databases\\NORMAN_Dust_Compounds.txt"
# Call the function
ClassFinder(path_file, path_harm, path_fuse, path_norm)



# Library search
# Load library and aligned overview files
println("Specifying library/database file... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
pathDB = "D:\\Projects\\Test\\Databases\\Dust_Database.xlsx"
println("Done specifying library/database. ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "\n")
println("Specifying aligned overview file... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
pathFile = "D:\\Projects\\Test\\Aligned_File.csv"
println("Done specifying overview file. ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "\n")

# Run library search
addRIwin = [30, 50, 50, 100]       # Additional RI window on top of the aligned minimum and maximum
specID = "CS"       # Which of the spectra should be matched: "BM" = best match, "CS" = consensus
mz_tol = 0.1    # m/z tolerance
numfrags = 50   # Minimum number of highest intensity fragments used from each spectrum (default = 15); can be left out for running the function
index = "all"      # Which features should be used for library matching: "all" uses every entry, a vector of numbers [1, 2, 6] uses only those indices, and collect(1:100) uses only indices/entries from 1 to 100
similarity_method = "DISCO"    # Spectral similarity method: "DISCO" (DIstance & Spectrum Correlation Optimization) or "NDP" (Normalized Dot Product)
librarySearch(pathDB, pathFile, addRIwin, specID, mz_tol, numfrags, index, similarity_method)



# Library search visualization
println("Specifying file for visualization... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
# pathFile = "D:\\Projects\\Test\\Aligned_File_LibSearch_CS.csv"
pathFile = "D:\\Projects\\Test\\251210_Aligned_Fltrd_Mrgd_PubChem_Class_LibSearch_CS.csv"
println("Done specifying file for visualization... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "\n")
index = "all"       # Which features should be plotted: "all" plots every entry, a vector of numbers [1, 2, 6] plots only those indices, and collect(1:100) plots only indices/entries from 1 to 100
libraryVisualization(pathFile, index)



# A function to resort columns in a specific order and delete specified columns from the aligned table(s)
# Specify file to be processed
fileColDrop = "D:\\Projects\\Test\\Aligned_File.csv"
# Specify metadata column order
meta_order = [
    "Nr","MinRt1D","MaxRt1D","AveRt1D","MinRt2D","MaxRt2D","AveRt2D",
    "MinRI","MaxRI","MedRI","MinRIlib","MaxRIlib","MedRIlib",
    "Name","Formula","BaseMass","InChIKey","InChIKey_PubChem","InChIKey_Consensus",
    "CAS","CAS_PubChem","CAS_Consensus","SMILES_PubChem",
    "Class_Harmonized","Class_FUse","Class_NORMAN",
    "Kingdom","Superclass","Class_CFB","Subclass",
    "Parent_Level_1","Parent_Level_2","Parent_Level_3","Parent_Level_4","Parent_Level_5",
    "Origin","Class","DetectionNum","DetectionFreq",
    "Area","TotalArea","Similarity","Probability","dRI","Score",
    "All_Names","All_Formulae","All_InChIKey","All_TotalArea","All_BaseMass",
    "All_Origin","All_Class","All_CAS","All_Similarity",
    "Spectrum_Consensus","Spectrum_BestMatch"
]
# Specify low-importance or QC columns to drop
col_names = [
    "MinRt1D", "MaxRt1D", "MinRt2D", "MaxRt2D", "MinRI", "MaxRI", "MinRIlib", "MaxRIlib",
    "All_Names", "All_Formulae", "All_InChIKey", "All_TotalArea", "All_BaseMass",
    "All_Origin", "All_Class", "All_CAS", "All_Similarity"
]
# Call the function
resort_delete_cols(path2files, fileColDrop, col_names, meta_order)

