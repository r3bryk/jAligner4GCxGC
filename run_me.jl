#################################################################################
################################ Installation ###################################
#################################################################################

# These lines (10 and 14) should be run in Julia REPL the first time to install the package and load it; 
# after that, you can just run the execution part (starting with "using jAligner4GCxGC") 
# without running the installation part again, unless you want to update the package or have made changes to the code.

# Enable using packages
using Pkg


# Get package
Pkg.add(url="https://github.com/r3bryk/jAligner4GCxGC")



#################################################################################
################################# Execution #####################################
#################################################################################

# From this line onwards, the code can be run as a script (e.g., in VS Code) 
# without running the installation part again, unless you want to update the package or have made changes to the code.

# Before running any part of the script, specify all the file paths and/or parameters where needed.

# Import package
using jAligner4GCxGC



# Specify the path to input files
println("Specifying a folder with .txt files for processing. ", 
    jAligner4GCxGC.Dates.format(jAligner4GCxGC.now(), "yyyy-mm-dd HH:MM:SS"))

path2files = "test" # Path to the working folder with files to be aligned

println("Done. There are ", length(filter(f -> endswith(f, ".txt"), readdir(path2files))), " .txt files. ", 
    jAligner4GCxGC.Dates.format(jAligner4GCxGC.now(), "yyyy-mm-dd HH:MM:SS"), "\n")



# Convert TXT files to CSV (only needed the first time)
jAligner4GCxGC.convert_txt2csv(path2files)



# Filter out rows where "Classifications" from ChromaTOF contains 'Bleed', 'Toluene', or 'DCM'
jAligner4GCxGC.filter_classifications(path2files)



# Specify the keywords (unwanted compound names, parts of names, etc.) to filter out
keywords = [
    "1H-Tetrazol-5-amine", "1H-Tetrazole, 1-methyl-", "Acetic acid, mercapto-", 
    "1H-Pyrrole-2-carbonitrile", "5-Diazo-1,3-cyclopentadiene", "1H-Pyrrole-3-carbonitrile", 
    "1H-1,2,3-Triazole-4-", "1H-1,2,4-Triazole, 3-", "1-Propanamine, 3-dibenzo", 
    "2,4,6-Cycloheptatrien", "2-Picoline, 6-nitro-", "4-Benzyloxy-3-hydroxy", 
    "9-[(S)-2-(Hydroxymethyl)", "Arsan", "arsan", "Arsen", "arsen", "arsin",
    "Benzyl (1,2,3-thiadiazol", "Benzyl-4,9,15-trioxa", "Bora", "bora", 
    "Borin", "borin", "boro", "Boro", "boryl", "Chromium", "Cobalt", 
    "Iron,", "Mercury", "sila", "Sila", "Silic", "silic", "Silo", "silo", 
    "silver", "silyl", "TBDMS", "TMS", "triazolo[", "Tricyclo[3.2.2.0(2,4)]", 
    "Zinc", "Zirconium", "Beryllium"
]

# Filter out rows where "Name" contains unwanted values, e.g., siloxanes & other bleed
jAligner4GCxGC.filter_names(path2files, keywords)



# Specify the filters as a list of tuples (keyword, spectrum_number (specific integer m/z in spectrum))
filters = [
    ("(4aS,4bS,10aS)-1,1,4a-Trimethyl-7-(propan-2-ylidene)-1,2,3,4,4a,4b,5,6,7,9,10,10a-dodecahydrophenanthrene", 272),
    ("(Benzyloxy)(methyl)amine", 105),
    ("(Isopropoxymethy)lbenzene", 107),
    ("(R)-9-[(S)-2-(Hydroxymethyl)pyrrolidin-1-yl]-3-methyl-3,4-dihydro-2H-benzo[b][1,4,5]oxathiazepine 1,1-dioxide", 312),
    ("(S)-9-[(S)-2-(Hydroxymethyl)pyrrolidin-1-yl]-3-methyl-3,4-dihydro-2H-benzo[b][1,4,5]oxathiazepine 1,1-dioxide", 311),
    ("1-Benzyl-2-(trifluoromethyl)aziridine", 201),
    ("1-Benzyl-3-hydroxypyridinium hydroxide", 185),
    ("1-Benzylcyclopentanol-1", 129),
    ("1H-[1,2,3]Triazole-4-carboxylic acid, 5-acetylamino-1-benzyl-, phenylamide", 335),
    ("1-Phenyl-2-propanol", 136),
    ("2(3H)-Benzofuranone, 6-ethenylhexahydro-6-methyl-3-methylene-7-(1-methylethenyl)-, [3aS-(3aa,6a,7ß,7aß)]-", 232),
    ("2-(Benzyloxy)-1-chloro-3-fluorobenzene", 117),
    ("2-(Benzyloxy)-4-methoxybenzaldehyde", 242),
    ("2-Benzyl-3-isopropyl-cyclopentanone", 216),
    ("2-Benzyloxyethylamine", 105),
    ("2-Butanol, 3-benzyloxy-", 135),
    ("3-(Iodomethyl)pyridine", 219),
    ("3-Benzyl-4-chloro-1,2,3-triazole 1-oxide", 130),
    ("3-Benzylsulfanyl-3-fluoro-2-trifluoromethyl-acrylonitrile", 258),
    ("3H-Pyrazole, 5-ethynyl-3,3-dimethyl-", 120),
    ("3-Picoline, 2-nitro-", 108),
    ("4-(Benzyloxy)-3-fluorophenol, trifluoroacetate", 314),
    ("4-(Benzyloxy)pyridine 1-oxide", 201),
    ("4-Azido-2-phenylmethanesulfinyl-benzonitrile", 115),
    ("4-Azido-2-phenylmethanesulfonyl-benzonitrile", 234),
    ("4-Benzyloxy-2-methyl-2-buten-1-ol", 108),
    ("4H-1,2,4-triazol-3-ol, 5-[(phenylmethyl)thio]-", 207),
    ("5-Benzyloxy-2-nitrotoluene", 243),
    ("5-Diazo-1,3-cyclopentadiene", 63),
    ("7,8-Diazabicyclo[4.2.2]deca-2,4,7,9-tetraen-7-oxide", 118),
    ("7-Chloro-2,3-dihydro-3-(4-N,N-dimethylaminobenzylidene)-5-phenyl-1H-1,4-benzodiazepin-2-one", 159),
    ("Benzaldehyde, 4-methoxy-3-(phenylmethoxy)-", 242),
    ("Benzene, (1,1-dimethylnonyl)-", 232),
    ("Benzene, (1-azidoethyl)-", 147),
    ("Benzene, (1-butylheptyl)-", 147),
    ("Benzene, (1-butylnonyl)-", 147),
    ("Benzene, (1-ethyldecyl)-", 246),
    ("Benzene, (1-ethylundecyl)-", 260),
    ("Benzene, (1-methyldecyl)-", 232),
    ("Benzene, (1-methyldodecyl)-", 260),
    ("Benzene, (1-methylundecyl)-", 246),
    ("Benzene, (1-pentylheptyl)-", 246),
    ("Benzene, (1-pentyloctyl)-", 260),
    ("Benzene, (1-propyldecyl)-", 133),
    ("Benzene, (1-propylnonyl)-", 133),
    ("Benzene, (2,2-dichloroethyl)-", 174),
    ("Benzene, (2-chloroethyl)-", 140),
    ("Benzene, (2-chloropropyl)-", 154),
    ("Benzene, (2-cyclohexylethyl)-", 188),
    ("Benzene, (2-methylpropyl)-", 134),
    ("Benzene, (3-methylpentyl)-", 162),
    ("Benzene, (bromomethyl)-", 170),
    ("Benzene, (butoxymethyl)-", 107),
    ("Benzene, (ethoxymethyl)-", 135),
    ("Benzene, (iodomethyl)-", 127),
    ("Benzene, (phenoxymethyl)-", 184),
    ("Benzene, (propoxymethyl)-", 107),
    ("Benzene, [(2-propenyloxy)methyl]-", 107),
    ("Benzene, [(methylsulfinyl)methyl]-", 154),
    ("Benzene, [(methylsulfonyl)methyl]-", 170),
    ("Benzene, 1,1'-(1,1,2,2-tetramethyl-1,2-ethanediyl)bis-", 119),
    ("Benzene, 1,1'-[oxybis(methylene)]bis-", 107),
    ("Benzene, 1-methyl-2-nitroso-", 121),
    ("Benzene, 1-methyl-3-(1-methylethenyl)-", 132),
    ("Benzene, n-butyl-", 134),
    ("Benzeneacetaldehyde", 120),
    ("Benzeneacetamide", 135),
    ("Benzeneacetamide, N", 149),
    ("Benzeneacetic acid", 136),
    ("Benzeneacetic acid 1-methylethyl ester", 176),
    ("Benzeneacetic acid, 2-propenyl ester", 103),
    ("Benzeneethanol", 178),
    ("Benzenemethanesulfonamide", 107),
    ("Benzenemethanesulfonyl chloride", 126),
    ("Benzenemethanethiol", 124),
    ("Benzenepropanenitrile, a-phenyl-", 107),
    ("Benzenesulfonamide, 4-methyl-", 171),
    ("Benzonitrile, m-phenethyl-", 201),
    ("Benzyl (1,2,3-thiadiazol-4-y)carbamate", 108),
    ("Benzyl 2-chloroethyl sulfone", 218),
    ("Benzyl 4-nitrophenyl carbonate", 139),
    ("Benzyl butyl phthalate", 149),
    ("Benzyl chloride", 126),
    ("Benzyl chloroformate", 170),
    ("Benzyl isopentyl ether", 107),
    ("Benzyl lactate", 180),
    ("Benzyl methyl disulfide", 170),
    ("Benzyl methyl ketone", 134),
    ("Benzyl N-[4-(4-cyano-3-fluorophenyl)phenyl]carbamate, TFA", 238),
    ("Benzylcyclopentane", 160),
    ("Bibenzyl", 182),
    ("Bicyclo[2.2.2]oct-7-en-2-one, 5-methylene-", 134),
    ("Bicyclo[3.1.1]hept-2-ene, 3,6,6-trimethyl-", 121),
    ("Butane, 1-(benzyloxy)-2-[(benzyloxy)methyl]-", 193),
    ("Cycloheptatrienylium, iodide", 78),
    ("Decane, 1-chloro-", 105),
    ("Dispiro[cyclopropane-1,3'-tricyclo[5.2.1.0(2,6)]decane-10',1''-cyclopropane]-4',8'-diene", 115),
    ("Dodecane", 170),
    ("Dodecane, 1-chloro-", 105),
    ("Dodecane, 2,6,11-trimethyl-", 113),
    ("Ethane, hexachloro-", 201),
    ("Eucalyptol", 154),
    ("Hydrazinecarbothioamide", 60),
    ("Hydroxylamine, O-(phenylmethyl)-", 105),
    ("Isophorone", 138),
    ("MGK-264", 164),
    ("N-(Phenylacetyl)glycine", 193),
    ("N-Benzyloxy-2-carbomethoxyaziridine", 105),
    ("N-Hydroxymethyl-2-phenylacetamide", 165),
    ("Pentadecane", 212),
    ("Pentalene, 1,2,4,5,6,6a-hexahydro-2-methylene-", 120),
    ("Phenylacetamide", 118),
    ("Phenylacetamide, N-propyl-", 177),
    ("Phenylethyl Alcohol", 122),
    ("Phosphine oxide, bis(pentamethylphenyl)-", 342),
    ("Phthalan", 120),
    ("Pyridine, 2-methyl-, 1-oxide", 109),
    ("Pyridine, 4-methyl-2-nitro-", 138),
    ("Sabinyl, 2-methylbutanoate", 119),
    ("Spiro[4.4]non-3-en-2-one, 4-methyl-3-(1H-tetrazol-5-yl)-1-oxa-", 123),
    ("ß-Myrcene", 136),
    ("Sydnone, 3-(phenylmethyl)-", 176),
    ("tert-Nonylphenol, Ac derivative", 177),
    ("Thiocyanic acid", 149),
    ("Thiocyanic acid, phenylmethyl ester", 149),
    ("Tricyclo[3.2.2.0(2,4)]non-8-ene-6,6,7,7-tetracarbonitrile", 128)    
]

# Filter out rows where "Name" contains, e.g., compounds with spectra similar to toluene (m/z 91 or 92),
# if their spectra do not contain other m/z values, like molecular ions
jAligner4GCxGC.filter_names_if_missing_mz(path2files, filters)



# Specify the filters for "Peak" and bleed m/z as a list of tuples (keyword, spectrum_number (specific integer m/z in spectrum))
bleed_filters = [
    ("Peak", 73), ("Peak", 147), ("Peak", 207), ("Peak", 267), ("Peak", 281), ("Peak", 341), ("Peak", 355), ("Peak", 429),
    ("Peak", 59), ("Peak", 78), ("Peak", 135), ("Peak", 156), ("Peak", 193), ("Peak", 195), ("Peak", 209), ("Peak", 215),
    ("Peak", 221), ("Peak", 251), ("Peak", 253), ("Peak", 255), ("Peak", 269), ("Peak", 327), ("Peak", 329), ("Peak", 331),
    ("Peak", 333), ("Peak", 339), ("Peak", 343), ("Peak", 377), ("Peak", 401), ("Peak", 403), ("Peak", 405), ("Peak", 415),
    ("Peak", 417), ("Peak", 439), ("Peak", 451), ("Peak", 475), ("Peak", 477), ("Peak", 489), ("Peak", 503), ("Peak", 549),
    ("Peak", 553), ("Peak", 563), ("Peak", 623), ("Peak", 91), ("Peak", 92), ("Peak", 137), ("Peak", 479)
]

# Filter out rows where "Name" contains "Peak" if their spectra contain bleed m/z values
jAligner4GCxGC.filter_peak_names_with_bleed_mz(path2files, bleed_filters)



# SKIP: no need for the new files/new alignment
# Reset masses classified/filtered by filter_signal_external (previous classification)
# jAligner4GCxGC.reset_filters(path2files)
# SKIP: no need for the new files/new alignment



# Remove m/z:intensity pairs with intensity values containing '0.' from spectra
jAligner4GCxGC.remove_zero_intensity_pairs(path2files)



# Remove m/z:intensity pairs with intensity values <= 3% of the highest intensity value
jAligner4GCxGC.remove_low_intensity_pairs(path2files)



# Classify features by base masses or filter out features by base masses
# Specify parameters
filtType = "Bleed"   # String with keyword of what is being classified/filtered; will be assigned and stored in 'Origin' column
base = [73, 147, 207, 267, 281, 341, 355, 429, 479, 59, 78, 135, 156, 193, 195, 197, 209, 215, 221, 251, 253, 255, 269, 
    327, 329, 331, 333, 339, 343, 377, 401, 403, 405, 415, 417, 439, 451, 475, 477, 489, 503, 549, 553, 563, 623
]   # Base mass(es) used for scanning 'BaseMass' column (can be a vector [xx, xx, xx])
mz = []   # Other signal(s)/mass(es) that need to be present
int = []    # Relative intensity of those signals with 1 = 100% (if 1 is not provided the base mass intensity is used)
tr1_range = [0, 0]   # 1st dimension retention range mass(es) can occur in
tr2_range = [0, 0]   # 2nd dimension retention range mass(es) can occur in
mz_tol = 0.1    # m/z tolerance

# Call the function
jAligner4GCxGC.filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)


# Classify features by base masses or filter out features by base masses
# Specify parameters
filtType = "Toluene"   # String with keyword of what is being classified/filtered; will be assigned and stored in 'Origin' column
base = [91]   # Base mass(es) used for scanning 'BaseMass' column (can be a vector [xx, xx, xx])
mz = [92]   # Other signal(s)/mass(es) that need to be present
int = [0.5]    # Relative intensity of those signals with 1 = 100% (if 1 is not provided the base mass intensity is used)
tr1_range = [0, 0]   # 1st dimension retention range mass(es) can occur in
tr2_range = [0, 0]   # 2nd dimension retention range mass(es) can occur in
mz_tol = 0.1    # m/z tolerance

# Call the function
jAligner4GCxGC.filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)


# Classify features by base masses or filter out features by base masses
# Specify parameters
filtType = "DCM"   # String with keyword of what is being classified/filtered; will be assigned and stored in 'Origin' column
base = [49]   # Base mass(es) used for scanning 'BaseMass' column (can be a vector [xx, xx, xx])
mz = [84]   # Other signal(s)/mass(es) that need to be present
int = [0.5]    # Relative intensity of those signals with 1 = 100% (if 1 is not provided the base mass intensity is used)
tr1_range = [0, 0]   # 1st dimension retention range mass(es) can occur in
tr2_range = [0, 0]   # 2nd dimension retention range mass(es) can occur in
mz_tol = 0.1    # m/z tolerance

# Call the function
jAligner4GCxGC.filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)


# Classify features by base masses or filter out features by base masses
# Specify parameters
filtType = "Alkane"   # String with keyword of what is being classified/filtered; will be assigned and stored in 'Origin' column
base = [57]   # Base mass(es) used for scanning 'BaseMass' column (can be a vector [xx, xx, xx])
mz = [71, 85]   # Other signal(s)/mass(es) that need to be present
int = [0.3, 0.2]    # Relative intensity of those signals with 1 = 100% (if 1 is not provided the base mass intensity is used)
tr1_range = [0, 0]   # 1st dimension retention range mass(es) can occur in
tr2_range = [0, 0]   # 2nd dimension retention range mass(es) can occur in
mz_tol = 0.1    # m/z tolerance

# Call the function
jAligner4GCxGC.filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)


# Classify features by base masses or filter out features by base masses
# Specify parameters
filtType = "Acid"   # String with keyword of what is being classified/filtered; will be assigned and stored in 'Origin' column
base = [55, 60, 73]   # Base mass(es) used for scanning 'BaseMass' column (can be a vector [xx, xx, xx])
mz = [60, 73]    # Other signal(s)/mass(es) that need to be present
int = [0.5, 1]    # Relative intensity of those signals with 1 = 100% (if 1 is not provided the base mass intensity is used)
tr1_range = [0, 0]   # 1st dimension retention range mass(es) can occur in
tr2_range = [0, 0]   # 2nd dimension retention range mass(es) can occur in
mz_tol = 0.1    # m/z tolerance

# Call the function
jAligner4GCxGC.filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)



# SKIP: merges terpenes and other closely eluting compounds with similar spectra that should not be merged
# # Merge internal drifting feature
# # Specify parameters
# tr1_dev = [3, 4.9]    # Retention deviation right and left of the peak in 1st dimension
# tr2_dev = [0.05, 0.1] # Retention deviation up and down of the peak in 2nd dimension
# sim_thresh = 0.9  # Spectral similarity threshold for merging
# mz_tol = 0.025    # m/z tolerance for setting up the spectrum for dot product calculations

# Call the function
# jAligner4GCxGC.internal_feature_merging(path2files, tr1_dev, tr2_dev, sim_thresh, mz_tol)
# SKIP: merges terpenes etc. that should not be merged



# Align samples
# Specify parameters
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

# Call the function
jAligner4GCxGC.feature_align_external(path2files, rt2_tol_lowA, rt2_tol_highA, rt_tol_lowA, h_thresh, rt_tol_highA, 
    sim_thresh, mz_tol, overviewMerging, weights, tr1_parts; numfrags, similarity_method)



# Filter out rows after alignment, where "Name" contains, e.g., compounds with spectra similar to toluene (m/z 91 or 92),
# if their spectra do not contain other m/z values, like molecular ions
jAligner4GCxGC.filter_names_mz_postalign(path2files, filters)



# Merge features with the same name and base m/z if their dRI < N
jAligner4GCxGC.merge_postalign(path2files, 7.0)



# Retrieve InChIKey, SMILES & CAS# values from PubChem
jAligner4GCxGC.PubChemRetriever(path2files)



# Retrieve functional use data from EPA FUse, NORMAN dust trial & in-house harmonized use databases
# Specify paths to working file and 3 files with classification data
path_file = "D:\\Projects\\Test\\Aligned_File.csv"
isfile(path_file) ? println("Found working file:\n", path_file, "\n") : @warn("Working file not found:\n", path_file)

path_harm = "D:\\Projects\\Test\\Databases\\Classes_Harmonized.txt"
isfile(path_harm) ? println("Found harmonized class file:\n", path_harm, "\n") : @warn("Harmonized class file not found:\n", path_harm)

path_fuse = "D:\\Projects\\Test\\Databases\\FUse_Database.txt"
isfile(path_fuse) ? println("Found FUse file:\n", path_fuse, "\n") : @warn("FUse file not found:\n", path_fuse)

path_norm = "D:\\Projects\\Test\\Databases\\NORMAN_Dust_Compounds.txt"
isfile(path_norm) ? println("Found NORMAN file:\n", path_norm, "\n") : @warn("NORMAN file not found:\n", path_norm)

# Call the function
jAligner4GCxGC.ClassFinder(path_file, path_harm, path_fuse, path_norm)



# Library search
# Specify paths to library/database and aligned overview file
pathDB = "D:\\Projects\\Test\\Databases\\Dust_Database.xlsx"
isfile(pathDB) ? println("Found library file:\n", pathDB, "\n") : @warn("Library file not found:\n", pathDB)

pathFile = "D:\\Projects\\Test\\Aligned_File.csv"
isfile(pathFile) ? println("Found aligned file:\n", pathFile, "\n") : @warn("Aligned file not found:\n", pathFile)

# Specify parameters
addRIwin = [30, 50, 50, 100]    # Additional RI window on top of the aligned minimum and maximum
specID = "CS"   # Which of the spectra should be matched: "BM" = best match, "CS" = consensus
mz_tol = 0.1    # m/z tolerance
numfrags = 50   # Minimum number of highest intensity fragments used from each spectrum (default = 15); can be left out for running the function
index = "all"   # Which features should be used for library matching: "all" uses every entry, a vector of numbers [1, 2, 6] uses only those indices, and collect(1:100) uses only indices/entries from 1 to 100
similarity_method = "DISCO" # Spectral similarity method: "DISCO" (DIstance & Spectrum Correlation Optimization) or "NDP" (Normalized Dot Product)

# Run library search
jAligner4GCxGC.librarySearch(pathDB, pathFile, addRIwin, specID, mz_tol, numfrags, index, similarity_method)



# Library search visualization
# Specify file for visualization
pathFile = "D:\\Projects\\Test\\Aligned_File_LibSearch_CS.csv"
isfile(pathFile) ? println("Found file for visualization:\n", pathFile, "\n") : @warn("File for visualization not found:\n", pathFile)

# Run library search visualization
index = "all"   # Which features should be plotted: "all" plots every entry, a vector of numbers [1, 2, 6] plots only those indices, and collect(1:100) plots only indices/entries from 1 to 100
jAligner4GCxGC.libraryVisualization(pathFile, index)



# A function to resort columns in a specific order and delete specified columns from the aligned table(s)
# Specify file to be processed
fileColDrop = "D:\\Projects\\Test\\Aligned_File.csv"
isfile(fileColDrop) ? println("Found file:\n", fileColDrop, "\n") : @warn("File not found:\n", fileColDrop)

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
jAligner4GCxGC.resort_delete_cols(path2files, fileColDrop, col_names, meta_order)
