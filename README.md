# jAligner4GCxGC

## Co-authors

Denice van Herwerden, https://github.com/DvHerwerden.

Andriy Rebryk, https://github.com/r3bryk.

## Description

`jAligner4GCxGC` is an application for inter-sample GC×GC-MS peak alignment, data filtration, and data retrieval.

## Prerequisites

Before using the script, several applications/tools have to be installed:
1. Visual Studio Code; https://code.visualstudio.com/download.
2. Julia Programming Language; https://julialang.org/downloads/.
3. Julia Extension in Visual Studio Code > Extensions > Search “julia” > Press `Install`.
4. Python 3 Programming Language; https://www.python.org/downloads/.
5. Conda Package and Environment Management System; https://www.anaconda.com/download.
6. Pip (package installer for Python); open Command Prompt > Type or copy python.exe -m pip install --upgrade pip > Press Enter.
7. Selenium for Python; open Command Prompt > Type or copy pip install selenium > Press Enter.
8. Google Chrome; https://www.google.com/chrome/.
9. ChromeDriver; use the same version as your Google Chrome version; https://googlechromelabs.github.io/chrome-for-testing/.

Then, the packages and functions must be loaded as follows:
1. Open the script file, e.g., `jAligner4GCxGC.jl`, with Visual Studio Code and wait until Julia environment and extension are loaded.
2. Enable using packages by highlighting `using Pkg` in **line 2** and pressing `“Ctrl + Enter”`.
3. Install packages by highlighting **lines 5-19** and pressing `“Ctrl + Enter”`. This procedure will take some time. When done, you will see a message `julia>` in the terminal field; this message will appear after any operation is completed and Julia is ready to proceed.
4. Install PyCall and Conda related packages by highlighting **lines 22-24** and pressing `“Ctrl + Enter”`.
5. Set up PyCall by specifying the path to python.exe on your computer in **line 33**, e.g., `ENV["PYTHON"] = "C:\\Program Files (x86)\\Python311\\python.exe"`, where `"C:\\Python\\Python3\\python.exe”` is the path to your Python3 folder. Then, load it by highlighting **lines 33-35** and pressing `“Ctrl + Enter”`.
6. Mute installation of packages by highlighting **lines 5-35** and pressing `“Ctrl + /”`. This step is needed to avoid repeated installation of the packages every time the script is used.
7. Load/import packages for use by highlighting **lines 38-56** and pressing `“Ctrl + Enter”`.
8. Load all functions by highlighting **lines 65-2873** (everything between the headers `FUNCTIONS` and `EXECUTION`) and pressing `“Ctrl + Enter”`.
**NB!** All the steps, except for **steps 3-6**, must be repeated every time you open the script file.

## How to use the script

### To use the script, the following steps should be executed:
1. All the input files, that should be aligned, should be copied to an empty working folder, the path to which should be specified in the script at **line 2883** as `path2files = "D:\\Projects\\Test"`, where `"D:\\Projects\\Test”` is the path to your working folder. Load input files by highlighting **line 2882-2884** and pressing `“Ctrl + Enter”`.

2. Convert `TXT` files to `CSV` format by highlighting **line 2889** (`convert_txt2csv(path2files)`) and pressing `“Ctrl + Enter”`.

3. Filter out rows for which `"Classifications"` column contains `'Bleed'`, `'Toluene'`, or `'DCM'` by highlighting **line 2894** (`filter_classifications(path2files)`) and pressing `“Ctrl + Enter”`.

4. Filter out rows where `"Name"` contains unwanted values, e.g., `siloxanes` and other bleed:
- Specify keywords (unwanted compound names or parts of names) to filter out in the keywords array. The names should be in square brackets, double quoted, and separated by commas, e.g., `keywords = ["1H-Tetrazol-5-amine", "Silic", "TMS"]`. Then, load the keywords array by highlighting **lines 2900-2911** and pressing `“Ctrl + Enter”`.
- Filter out rows where `"Name"` contains unwanted values by highlighting **line 2913** (`filter_names(path2files, keywords)`) and pressing `“Ctrl + Enter”`.

5. Filter out rows where `"Name"` contains, e.g., compounds with spectra similar to toluene (`m/z 91` or `m/z 92`), if their spectra do not contain other m/z values (`m/z > 91/92`) such as molecular ions:
- Specify filters as a list of tuples (keyword, spectrum_number (specific integer m/z in spectrum)), e.g., `filters = [("(Isopropoxymethy)lbenzene", 107), ("2-Butanol, 3-benzyloxy-", 135), ("Bibenzyl", 182)]`. Then, load the filters by highlighting **lines 2920-3047** and pressing `“Ctrl + Enter”`.
- Filter out rows by highlighting **line 3049** (`filter_names_if_missing_mz(path2files, filters)`) and pressing `“Ctrl + Enter”`.

6. Filter out rows where `"Name"` contains `"Peak"` (a marker for unknown feature) if their spectra contain bleed m/z values such as `73` or `207`:
- Specify bleed filters as a list of tuples (keyword, spectrum_number (specific integer m/z in spectrum)), e.g., `bleed_filters = [("Peak", 73), ("Peak", 147), "Peak", 207)]`. Then, read the bleed filters by highlighting **lines 3054-3061** and pressing `“Ctrl + Enter”`.
- Filter out rows by highlighting **line 3063** (`filter_peak_names_with_bleed_mz(path2files, bleed_filters)`) and pressing `“Ctrl + Enter”`.

7. If needed, reset masses classified/filtered by `filter_signal_external` function (previous classification; see below) by unmuting **line 3069** (`reset_filters(path2files)`) by highlighting it and pressing `“Ctrl + /”` and run it by pressing `“Ctrl + Enter”`.

8. Remove `m/z:intensity` pairs with intensity values containing `'0.'` from spectra by highlighting **line 3075** (`remove_zero_intensity_pairs(path2files)`) and pressing `“Ctrl + Enter”`.

9. Remove `m/z:intensity` pairs with intensity values ≤ 3% of the highest intensity value by highlighting **line 3080** (`remove_low_intensity_pairs(path2files)`) and pressing `“Ctrl + Enter”`.

10. Classify features as `"Bleed"` by base masses or filter out features by base masses or mass/intensity combinations: 
- Specify in **lines 3085-3093** the input parameters, such as: 
(i) Base mass(es) (`base = [73, 147, 207]`).
(ii) Other signals that need to be present (mz).
(iii) Relative intensity of those signals with 1 = 100% (int).
(iv) 1st and 2nd dimension retention ranges mass(es) can occur in (`tr1_range`, `tr2_range`).
(v) m/z tolerance (`mz_tol = 0.1`).
(vi) Load the parameters by highlighting **lines 3085-3093** and pressing `“Ctrl + Enter”`.
- Classify features by highlighting **line 3094** (`filter_signal_external(path2files, filtType, base, mz,int,  mz_tol, tr1_range, tr2_range)`) and pressing `“Ctrl + Enter”`.

11. Classify features as `"Toluene"` by base masses or filter out features by base masses or mass/intensity combinations: 
- Specify in **lines 3097-3103** the input parameters, such as: 
(i) Base mass(es) (`base = [91]`).
(ii) Other signals that need to be present (`mz = [92]`).
(iii) Relative intensity of those signals with 1 = 100% (`int = [0.5]`).
(iv) 1st and 2nd dimension retention ranges mass(es) can occur in (`tr1_range`, `tr2_range`).
(v) m/z tolerance (`mz_tol = 0.1`).
(vi) Load the parameters by highlighting **lines 3097-3103** and pressing `“Ctrl + Enter”`.
- Classify features by highlighting **line 3104** (`filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)`) and pressing `“Ctrl + Enter”`.

12. Classify features as `"DCM"` by base masses or filter out features by base masses or mass/intensity combinations: 
- Specify in **lines 3107-3113** the input parameters, such as: 
(i) Base mass(es) (`base = [49]`).
(ii) Other signals that need to be present (`mz = [84]`).
(iii) Relative intensity of those signals with 1 = 100% (`int = [0.5]`).
(iv) 1st and 2nd dimension retention ranges mass(es) can occur in (`tr1_range`, `tr2_range`).
(v) m/z tolerance (`mz_tol = 0.1`).
(vi) Load the parameters by highlighting **lines 3107-3113** and pressing `“Ctrl + Enter”`.
- Classify features by highlighting **line 3114** (`filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)`) and pressing `“Ctrl + Enter”`.

13. Classify features as `"Alkane"` by base masses or filter out features by base masses or mass/intensity combinations: 
- Specify in **lines 3117-3123** the input parameters, such as:
(i) Base mass(es) (`base = [57]`).
(ii) Other signals that need to be present (`mz = [71, 85]`).
(iii) Relative intensity of those signals with 1 = 100% (`int = [0.3, 0.2]`).
(iv) 1st and 2nd dimension retention ranges mass(es) can occur in (`tr1_range`, `tr2_range`).
(v) m/z tolerance (`mz_tol = 0.1`).
(vi) Load the parameters by highlighting **lines 3117-3123** and pressing `“Ctrl + Enter”`.
- Classify features by highlighting **line 3124** (`filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)`) and pressing `“Ctrl + Enter”`.

14. Classify features as `"Acid"` by base masses or filter out features by base masses or mass/intensity combinations: 
- Specify in **lines 3127-3133** the input parameters, such as: 
(i) Base mass(es) (`base = [55, 60, 73]`).
(ii) Other signals that need to be present (`mz = [60, 73]`).
(iii) Relative intensity of those signals with 1 = 100% (`int = [0.5, 1]`).
(iv) 1st and 2nd dimension retention ranges mass(es) can occur in (`tr1_range`, `tr2_range`).
(v) m/z tolerance (`mz_tol = 0.1`).
(vi) Load the parameters by highlighting **lines 3127-3133** and pressing `“Ctrl + Enter”`.
- Classify features by highlighting **line 3134** (`filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)`) and pressing `“`Ctrl + Enter”`.

15. If needed (not recommended), merge internal drifting features (features of the same chemical origin, e.g., potential main peak and tailing peak):
- Specify in **lines 3140-3143** the input parameters, such as: 
(i) Retention deviation right and left of the peak in the 1st dimension (`tr1_dev = [3, 5]`).
(ii) Retention deviation up and down of the peak in 2nd dimension (`tr2_dev = [0.05, 0.1]`).
(iii) Spectral similarity threshold for merging (`sim_thresh = 0.8`).
(iv) m/z tolerance (`mz_tol = 0.1`).
(v) Load the parameters by highlighting **lines 3140-3143** and pressing `“Ctrl + Enter”`.
- Merge internal drifting features by highlighting **line 3144** (`internal_feature_merging(path2files, tr1_dev, tr2_dev, sim_thresh, mz_tol)`) and pressing `“Ctrl + Enter”`.

16. Align samples:
- Specify in **lines 3150-3161** the input parameters, such as: 
(i) Area value above which the highA RT threshold is used (`h_thresh = 10000000`).
(ii) 1st dimension RT tolerance for low intensity peaks (`rt_tol_lowA = 7`).
(iii) 1st dimension RT tolerance for high intensity peaks (`rt_tol_highA = 11`).
(iv) 2nd dimension RT tolerance for low intensity peaks (`rt2_tol_lowA = 0.1`).
(v) 2nd dimension RT tolerance for high intensity peaks (`rt2_tol_highA = 0.15`).
(vi) Spectral similarity threshold for alignment (`sim_thresh = 0.8`).
(vii) m/z tolerance (`mz_tol = 0.1`).
(viii) Method how the consensus spectrum is setup; `"height"` weights the spectrum based on the base peak height, `"equal"` assumes equal contribution (`overviewMerging = "height"`).
(ix) Weights used for score calculation [similarity, library probability, difference RI (library vs. experimental)] (`weights = [1, 1, 1]`).
(x) Number of parts the 1st RT range is split into for alignment (`tr1_parts = 5`).
(xi) Load the parameters by highlighting **lines 3150-3161** and pressing `“Ctrl + Enter”`.
- Align samples by highlighting **line 3162** (`feature_align_external(path2files, rt2_tol_lowA, rt2_tol_highA, rt_tol_lowA, h_thresh, rt_tol_highA, sim_thresh, mz_tol, overviewMerging, weights, tr1_parts; numfrags, similarity_method)`) and pressing `“Ctrl + Enter”`.

17. Filter out after alignment rows for which `"Name"` contains, e.g., compounds with spectra similar to toluene’s (`m/z 91` or `m/z 92`), if their spectra do not contain other m/z values (`m/z > 91/92`) such as molecular ions by highlighting **line 3168** (`filter_names_mz_postalign(path2files, filters)`) and pressing `“Ctrl + Enter”`.

18. Merge features with the same name and base m/z if their `ΔRI < 7.0` by highlighting **line 3173** (`merge_postalign(path2files, 7.0)`) and pressing `“Ctrl + Enter”`.

19. Retrieve `InChIKey`, `SMILES`, and `CAS#` values from PubChem by highlighting **line 3178** (`PubChemRetriever(path2files::String)`) and pressing `“Ctrl + Enter”`.

20. Retrieve chemical class data from ClassyFire Batch webpage: 
- Specify the path to chromedriver.exe in **line 3184** (`path2chromedriver = "C:\\Program Files (x86)\\Google\\chromedriver-win64\\chromedriver.exe"`, where `"C:\\Program Files (x86)\\Google\\chromedriver-win64\\chromedriver.exe"` is the path to the folder with chromedriver.exe). Then, load it by highlighting **line 3184** and pressing `“Ctrl + Enter”`.
- Retrieve chemical class data by highlighting **line 3186** (`ClassyFire_Batch(path2files, path2chromedriver)`) and pressing `“Ctrl + Enter”`.

21. Retrieve functional use data from, e.g., EPA FUse, NORMAN dust trial, and in-house harmonized use databases, or any other database(s) at hand:
- Specify in **lines 3192-3195** paths to the working file and the 3 files with classification data as:
(i) `path_file = "D:\\Projects\\Test\\Aligned_File.csv"`. Required columns: `Name`, `InChIKey_Consensus`, and `CAS_PubChem`.
(ii) `path_harm = "D:\\Projects\\Test\\Databases\\Classes_Harmonized`.txt". Required columns: `Name`, `InChIKey`, `Class`, and `CAS`.
(iii) `path_fuse = "D:\\Projects\\Test\\Databases\\FUse_Database.txt"`. Required columns: `Name`, `InChIKey`, `harmonized_function`, and `CAS`.
(iv) `path_norm = "D:\\Projects\\Test\\Databases\\NORMAN_Dust_Compounds.txt"`. Required columns: `Name`, `InChIKey`, `Class`, and `CAS`.
(v) Load the files by highlighting **lines 3192-3195** and pressing `“Ctrl + Enter”`.
- Retrieve functional use data by highlighting **line 3197** (`ClassFinder(path_file, path_harm, path_fuse, path_norm)`) and pressing `“Ctrl + Enter”`.

22. Search library/database:
- Specify library/database and aligned files as: 
(i) `pathDB = "D:\\Projects\\Test\\Dust_Database.xlsx"` in **line 3204**. Required columns: `"Spectrum"`, `"#"`, `"Name"`, `"RI"`, `"CAS Number"`, `"InChIKey"`, `"Formula"`.
(ii) `path_file = "D:\\Projects\\Test\\Aligned_File.csv"` in **line 3207**. Required columns: `"Name"` and/or `"All_Names"`, `"Spectrum_BestMatch"` and/or `"Spectrum_Consensus"`, `"BaseMass"` and/or `"All_BaseMass"`, `"Origin"` and/or `"All_Origin"`, `"Class"` and/or `"All_Class"`, `"MinRIlib"`, `"MinRI"`, `"MaxRIlib"`, `"MaxRI"`, `"Nr"`.
(iii) Load the files by highlighting **lines 3203-3208** and pressing `“Ctrl + Enter”`.
- Specify in **lines 3211-3216** the input parameters, such as: 
(i) RI windows for semi-standard non-polar RI, standard non-polar RI, standard polar RI, and estimated RI (`addRIwin = [30, 50, 50, 100]`, respectively). If you do not have any of the values for a particular compound, leave it blank; all the RI values can be left blank, but then the match will be based only on spectral similarity.
(ii) Which of the spectra should be matched: `"BM"` = best match, `"CS"` = consensus (`specID = "CS"`), m/z tolerance (`mz_tol = 0.1`).
(iii) Minimum number of highest intensity fragments used from each spectrum; default = 15 (`numfrags = 50`).
(iv) Which features should be used for library matching: `"all"` uses every entry, a vector of numbers `[1, 2, 6]` uses only those indices, and `collect(1:100)` uses only indices/entries from 1 to 100 (`index = "all"`).
(v) Load the parameters by highlighting **lines 3211-3216** and pressing `“Ctrl + Enter”`.
- Run library search by highlighting **line 3217** (`librarySearch(pathDB, pathFile, addRIwin, specID, mz_tol, numfrags, index, similarity_method)`) and pressing `“Ctrl + Enter”`.

23. Create head-to-tail graphs for library vs. experimental (user) spectra:
- Specify the input parameters, such as:
(i) Path to file with library search results, e.g., `pathFile = "D:\\Projects\\Test\\Aligned_File_LibSearch_CS.csv”` in **line 3224**.
(ii) Which features should be used for library matching in **line 3226**: `"all"` uses every entry, a vector of numbers `[1, 2, 6]` uses only those indices, and `collect(1:100)` uses only indices/entries from 1 to 100 (`index = "all"`).
(iii) Load the parameters by highlighting **lines 3222-3226** and pressing `“Ctrl + Enter”`.
- Run visualization by highlighting **line 3227** (`libraryVisualization(pathFile, index)`) and pressing `“Ctrl + Enter”`.

24. If needed, resort columns in a specific order and delete specified columns of low importance or QC columns from the aligned table:
- Specify in **lines 3233-3253** the input parameters, such as:
(i) Path to file to be processed, e.g., `fileColDrop = "D:\\Projects\\Test\\Aligned_File.csv"` in **line 3233**.
(ii) Metadata column order, e.g., `meta_order = ["Nr","MinRt1D","MaxRt1D","AveRt1D","MinRt2D","MaxRt2D","AveRt2D","Spectrum_BestMatch"]`.
(iii) Names of the columns to drop, e.g., `col_names = ["MinRt1D", "MaxRt1D", "MinRt2D", "MaxRt2D"]`.
(iv) Load the parameters by highlighting **lines 3233-3253** and pressing `“Ctrl + Enter”`.
- Run the function by highlighting **line 3255** (`resort_delete_cols(path2files, fileColDrop, col_names, meta_order)`) and pressing `“Ctrl + Enter”`.

## Notes and recommendations


It is recommended to use all the steps, except for `steps 7` and `15`, as some of the functions use the files or columns created by functions in previous steps by their specific filenames or column names, and/or extensions. Hence, if one of those functions was not used and the respective file was not created, the script will crash and return error. Alternatively, in case a function should be skipped, the file(s) can be renamed manually and then used by the respective function. For input filenames and suffices, see the code of the function in question, the part with input file specification, e.g., `if occursin("_Mrgd.csv", f)`.

The following functions can be skipped with no harm to script execution:
1. `filter_classifications(path2files)`
2. `filter_names(path2files, keywords)`
3. `filter_names_if_missing_mz(path2files, filters)`
4. `filter_peak_names_with_bleed_mz(path2files, bleed_filters)`
5. `reset_filters(path2files)`
6. `remove_zero_intensity_pairs(path2files)`
7. `remove_low_intensity_pairs(path2files)`
8. `internal_feature_merging(path2files, tr1_dev, tr2_dev, sim_thresh, mz_tol)`

The script takes `TXT (tab-delimited)` files, such as LECO ChromaTOF result files, as input.

The input files must contain at least the following columns to be processed:
`"Name"`, `"R.T. (s)"`, `"Area"`, `"Spectrum"`, `"Formula"`, `"CAS"`, `"InChIKey"`, `"Retention Index"`, `"Lib. RI"`, `"Height"`, `"Similarity"`, `"Base Mass"`, `"Actual Masses"`, `"FWHH (s)"`, `"Quant Masses"`, `"Probability"`, `"Expected Ion m/z"`, `"Observed Ion m/z"`, `"Classifications"`

**NB!** `"Spectrum"` values must be in LECO ChromaTOF result file format, i.e. `37:1070 38:1930 39:4500 40:3720 41:4570 49:30 50:130 51:90 52:220 62:50 63:120 64:80 67:9999 68:460`.

The library/database file must contain at least the following columns to be processed:
`"#"`, `"Name"`, `"InChIKey"`, `"CAS Number"`, `"Formula"`, `"Spectrum"`, `"RI"`, `"RI StdNP"`, `"RI StdPolar"`, `"RI Est"`, where `"RI"` is semi-standard non-polar RI, `"RI StdNP"` is standard non-polar RI, `"RI StdPolar"` is standard polar RI, and `"RI Est"` is estimated RI. If you do not have any of the RI values for a particular compound, leave the respective column cell blank; all the RI values can be left blank, but then the match will be based only on spectral similarity.

**NB!** `"Spectrum"` values must be in LECO ChromaTOF result file format, i.e. `37:1070 38:1930 39:4500 40:3720 41:4570 49:30 50:130 51:90 52:220 62:50 63:120 64:80 67:9999 68:460`.

## License
[![MIT License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/license/mit)

Intended for academic and research use.