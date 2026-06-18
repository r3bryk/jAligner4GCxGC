module jAligner4GCxGC

# # Enable using packages
# using Pkg

# Import packages
using PlotlyJS
using CSV
using DataFrames
using Glob
using Statistics
using LinearAlgebra
using XLSX
using Dates
using Unicode
using JSON3
using HTTP
#using PubChem
using PubChemCrawler
using JSON
#using WebDriver
#using PyCall

# Import selenium with PyCall
# pyimport_conda("selenium.webdriver", "selenium")



#################################################################################
################################# Functions #####################################
#################################################################################

# A function to import, clean-up, and align data from multiple .csv files, 
# ensuring consistency and handling missing values appropriately
function report_import_external(path2files)
    println("Running external alignment. Make sure that the directory contains only relevant .csv files that require alignment.")
    nn = readdir(path2files)
    Rep = DataFrame()
    namesT = []
    c = 1

    # Load files into a single dataframe
    for i =1:size(nn, 1)
        # println(i)
        m = split(nn[i, :][1], ".")
        if length(m[1]) == 0 || m[end] != "csv" || any(contains.(m, "Aligned"))  || any(contains.(m, "Original")) || any(contains.(m, "REMOVED")) || any(contains.(m, "merged")) || any(contains.(m, "IntMerged")) || any(contains.(m, "_IDresults")) || any(contains.(m, "_IDbest")) || any(contains.(m, "_IDconsens"))
            continue
        end 
        f_n = joinpath(path2files, nn[i, :][1])
        df = CSV.read(f_n, DataFrame, stringtype = String, select = Symbol.(["F_ID", "Name", "R.T. (s)", "Area", "Spectrum", "Formula", "CAS", "InChIKey", "Retention Index", "Lib. RI", "Height", "Similarity", "Base Mass", "Actual Masses", "FWHH (s)", "Quant Masses", "Probability", "Expected Ion m/z", "Observed Ion m/z", "Classifications", "Origin"]))
        if !any(names(df) .== "Origin")
            df[!,"Origin"] = fill("NA", size(df, 1))
        end
        df_ = df
        df_ = df_[vec(.! all(ismissing.(Matrix(df_)), dims = 2)), :]
        df1 = hcat(DataFrame(F_ID = c .* Int.(ones(size(df_, 1)))), df_)
        Rep = vcat(Rep, df1)
        namesT = vcat(namesT, nn[i])
        c += 1 
    end 

    # Remove missing 'Area' entries
    Rep = Rep[ismissing.(Rep[!, "Area"]) .== 0, :]

    # Convert retention times
    rt1, rt2 = extractRTs(Rep)
    Rep[ismissing.(Rep[:, "CAS"]), "CAS"] .= "NA"
    Rep[ismissing.(Rep[:, "InChIKey"]), "InChIKey"] .= "NA"
    Rep[ismissing.(Rep[:, "Formula"]), "Formula"] .= "NA"

    # Fix 'Similarity' score
    simScore = Rep[!, "Similarity"]
    Rep[ismissing.(Rep[!, "Similarity"]), "Similarity"] .= -1

    # Convert 'Peak #' to 'Peak_rt1_rt2'
    ind = findall(contains.(Rep[!, "Name"], "Peak "))
    Rep[ind, "Name"] .= "Peak_" .* string.(Int.(round.(rt1[ind]))) .* "_" .* string.(rt2[ind])

    # Clean-up InChIKeys
    inchi = inchiExtraction(Rep[:, "InChIKey"])
    inchiTrunc = fill("NoID", size(inchi, 1))
    inchiTrunc[inchi .!= "NoID"] = stack(split.(inchi[inchi .!= "NoID"], "-"), dims = 1)[:, 1]
   
    # Ensure library RI values are numbers only
    indRI = findall(.! ismissing.(Rep[!, "Lib. RI"]))
    for r in indRI
       Rep[r, "Lib. RI"] = split(replace(Rep[r, "Lib. RI"], ":" => ""), " ")[1]
    end
    Rep[ismissing.(Rep[!, "Lib. RI"]), "Lib. RI"] .= "NaN"
    libRI = parse.(Float64, Rep[:, "Lib. RI"])
    dRI = Rep[!, "Retention Index"] .- libRI

    # Create a new 'Probability' column with the correct type and copy the values
    new_probability = Vector{Union{Float64, Missing}}(undef, size(Rep, 1))
    for i in 1:size(Rep, 1)
        new_probability[i] = Rep[i, :Probability]
    end
    # Replace the original column with the new column
    Rep[!, :Probability] = new_probability

    Rep[ismissing.(Rep[:, "Probability"]), "Probability"] .= NaN
    Rep[ismissing.(Rep[:, "Expected Ion m/z"]), "Expected Ion m/z"] .= NaN
    Rep[ismissing.(Rep[:, "Observed Ion m/z"]), "Observed Ion m/z"] .= NaN
    if all(ismissing.(Rep[:, "Classifications"]))
        class = fill("NA", size(Rep, 1))
    else
        Rep[ismissing.(Rep[:, "Classifications"]), "Classifications"] .= "NA"
        class = Rep[:, "Classifications"]
    end

    # Total component area calculation
    totArea = Rep[!, "Area"]

    # Convert to general format (Rt = 1st dimension RT and Rt_2D = 2nd dimension RT)
    Rep = DataFrame(F_ID = Rep[:, "F_ID"], Name = Rep[:, "Name"], Int = totArea, Rt = rt1, Rt_2D = rt2, 
                    Spec = Rep[:, "Spectrum"], Form = Rep[:, "Formula"], CAS = Rep[:, "CAS"], Inchi = inchi, 
                    InchiTrunc = inchiTrunc, RI = Rep[!, "Retention Index"], RIlib = libRI, Height = Rep[!, "Height"], 
                    Sim = simScore, dRI = dRI, FWHH = Rep[!, "FWHH (s)"], libProb = Rep[!, "Probability"], 
                    ExpIon = Rep[!, "Expected Ion m/z"], ObsIon = Rep[!, "Observed Ion m/z"], BaseMass = Rep[!, "Base Mass"], ActualMass = Rep[!, "Actual Masses"], QuantMass = Rep[!, "Quant Masses"], Class = class, Origin = Rep[!, "Origin"])
   
    return(Rep, namesT)  
end



# A function to calculate total component area
function calcTotalArea(area, qmass, spec)
    totArea = zeros(length(area))

    for i = 1:length(totArea)
        sp = getSpec(spec[i])
        qm = parse(Float64, split(qmass[i], "\xb10.5")[1][5:end])
        qmAb = sp[argmin(abs.(sp[:, 1] .- qm)), 2]
        totAb = sum(sp[:, 2])
        totArea[i] = area[i] / (qmAb/totAb)
    end

    return totArea
end



# A function to extract and clean-up InChIKey values
function inchiExtraction(inchi)
    for i = 1:length(inchi)
        if inchi[i] == "NA"
            inchi[i] = "NoID"
        elseif contains(inchi[i], "href=")
            inchi[i] = split(inchi[i], "\">")[2][1:end-4]
        else
            println("Unknown InChIKey provided")
            inchi[i] = inchi[i]
        end
    end

    return inchi
end



# A function to extract RT1 and RT2 values
function extractRTs(Rep)
    rts = zeros(size(Rep, 1), 2)

    for i = 1:size(Rep, 1)
        rts[i, :] = parse.(Float64, split(Rep[i, "R.T. (s)"], ", "))
    end

    return rts[:, 1], rts[:, 2]
end



# A function to select candidate features for RT1
function select_candidates_rt1_(rep, ind, rt_tol)
    pf = rep[ind, :]

    lb = max(0, pf.Rt - rt_tol)
    ub = min(maximum(rep.Rt), pf.Rt + rt_tol)

    sel_inds = findall(x -> ub >= x >= lb, rep.Rt)

    return sel_inds
end



# A function to select candidate features for RT2
function select_candidates_rt2_(rep, ind, sel_inds, rt_tol)
    pf = rep[ind, :]

    lb = max(0, pf.Rt_2D - rt_tol)
    ub = min(maximum(rep.Rt_2D), pf.Rt_2D + rt_tol)

    sel_inds_f = sel_inds[findall(x -> ub >= x >= lb, rep.Rt_2D[sel_inds])]

    return sel_inds_f
end



# A function to select candidate features based on spectra similarity
function select_candidates_dot_(rep, ind, sel_inds_f, sim_thresh, mz_tol; numfrags::Int = 150, similarity_method::String = "DISCO")
    pf = rep[ind, :]
    vals = zeros(length(sel_inds_f))
    pfs = getSpec(rep[ind, "Spec"])

    for v = 1:length(sel_inds_f)
        spec = getSpec(rep[sel_inds_f[v], "Spec"])
        master_vect_f = master_mass_gen(pfs[:, 1], pfs[:, 2], spec[:, 1], spec[:, 2], mz_tol)
        reverse_vect_f = master_mass_gen(spec[:, 1], spec[:, 2], pfs[:, 1], pfs[:, 2], mz_tol)

        if size(master_vect_f, 1) > numfrags
            master_vect_fu = master_vect_f[sortperm(vec(sum(master_vect_f[:, 3], dims = 2)), rev = true), :][1:numfrags, :]
            master_vect_fr = master_vect_f[sortperm(vec(sum(master_vect_f[:, 2], dims = 2)), rev = true), :][1:numfrags, :]
            master_vect_f = [master_vect_fu; master_vect_fr]
            master_vect_f = unique(master_vect_f, dims = 1)
        end

        if size(reverse_vect_f, 1) > numfrags
            reverse_vect_fu = reverse_vect_f[sortperm(vec(sum(reverse_vect_f[:, 3], dims = 2)), rev = true), :][1:numfrags, :]
            reverse_vect_fr = reverse_vect_f[sortperm(vec(sum(reverse_vect_f[:, 2], dims = 2)), rev = true), :][1:numfrags, :]
            reverse_vect_f = [reverse_vect_fu; reverse_vect_fr]
            reverse_vect_f = unique(reverse_vect_f, dims = 1)
        end

        # Apply chosen spectral similarity method
        if similarity_method == "DISCO"
            # Mean center each vector
            master_vect_f[:, 2:3] = master_vect_f[:, 2:3] .- mean(master_vect_f[:, 2:3], dims = 1)
            reverse_vect_f[:, 2:3] = reverse_vect_f[:, 2:3] .- mean(reverse_vect_f[:, 2:3], dims = 1)

            # Calculate similarity for both directions and take the minimum
            similarity = min(
                dot(master_vect_f[:, 2], master_vect_f[:, 3]) / (norm(master_vect_f[:, 2]) * norm(master_vect_f[:, 3])),
                dot(reverse_vect_f[:, 2], reverse_vect_f[:, 3]) / (norm(reverse_vect_f[:, 2]) * norm(reverse_vect_f[:, 3]))
            )
        elseif similarity_method == "NDP"
            # Calculate similarity for both directions and take the minimum
            similarity = min(
            sum(master_vect_f[:, 2] .* master_vect_f[:, 3]) / (sqrt(sum(master_vect_f[:, 2] .^ 2)) * sqrt(sum(master_vect_f[:, 3] .^ 2))),
            sum(reverse_vect_f[:, 2] .* reverse_vect_f[:, 3]) / (sqrt(sum(reverse_vect_f[:, 2] .^ 2)) * sqrt(sum(reverse_vect_f[:, 3] .^ 2)))
        )
        else
            error("Unknown similarity method provided: $similarity_method. Use \"DISCO\" or \"NDP\".")
        end
        
        vals[v] = max(0, similarity)    # Set negative similarity scores to zero
    end

    sel_inds_s = sel_inds_f[vals .>= sim_thresh]

    return sel_inds_s
end



# A function to split spectrum and turn it into a 2D array
function getSpec(sp)
    # Remove trailing space if present
    if sp[end] == ' '
        sp = sp[1:end-1]
    end
    # Split the string by spaces and then by colons
    sp = split(sp, " ")
    sp = split.(sp, ":")
    # Parse the split strings into Float64 and reshape into a 2D array
    sp = parse.(Float64, stack(sp, dims = 1))

    return sp
end



# A function to generate a master list of m/z values and their corresponding intensities from two sets of input data
function master_mass_gen(mz_ref, int_ref, mz_user, int_user, mass_tol)
    mz_values = vcat(mz_ref[:], mz_user[:])
    master_vect = zeros(length(mz_values), 3)

    for i = 1:size(master_vect, 1)
        tv1 = abs.(mz_values[i] .- mz_values)
        tv2 = findall(x -> x <= mass_tol / 2, tv1)
        master_vect[i, 1] = mean(mz_values[tv2])
        mz_values[tv2] .= 0
    end

    for i = 1:size(master_vect, 1)
        if master_vect[i, 1] > 0
            tv1 = abs.(master_vect[i, 1] .- mz_ref)
            tv2 = findall(x -> x <= mass_tol / 2, tv1)
            if length(tv2) > 0
                master_vect[i, 2] = maximum(int_ref[tv2])
                int_ref[tv2] .= 0
            end

            ttv1 = abs.(master_vect[i, 1] .- mz_user)
            ttv2 = findall(x -> x <= mass_tol / 2, ttv1)
            if length(ttv2) > 0
                if maximum(ttv2) <= length(int_user)
                    master_vect[i, 3] = maximum(int_user[ttv2])
                    int_user[ttv2] .= 0
                elseif maximum(ttv2) > length(int_user)
                    for j = 1:length(ttv2)
                        if ttv2[j] > length(int_user)
                            ttv2[j] = 0
                        end
                    end
                    master_vect[i, 3] = maximum(int_user[ttv2[ttv2 .> 0]])
                    int_user[ttv2[ttv2 .> 0]] .= 0
                end
            end
        end
    end

    master_vect_f = master_vect[master_vect[:, 1] .> 0, :]

    return master_vect_f
end



# A function to group candidates
function group_candidates_external(reps, ind_files)
    if isempty(reps)
        return Inf, -Inf, NaN, Inf, -Inf, NaN, Inf, -Inf, NaN, Inf, -Inf, NaN, zeros(1, length(ind_files)), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files))), fill("NA", (1, length(ind_files)))
    end

    # Obtain RT and RI values for group
    MinRt = minimum(reps.Rt; init = Inf)
    MaxRt = maximum(reps.Rt; init = -Inf)
    AveRt = round(mean([MinRt,MaxRt]), digits = 2)
    MinMass = minimum(reps.Rt_2D; init = Inf)
    MaxMass = maximum(reps.Rt_2D; init = -Inf)
    AveMass = round(mean(reps.Rt_2D), digits = 4)
    MinRI = minimum(reps.RI; init = Inf)
    MaxRI = maximum(reps.RI; init = -Inf)
    MedRI = isempty(reps.RI) ? NaN : round(median(reps.RI), digits = 1)
    sel_inds_slib = findall(.! isnan.(reps.RIlib))
    if all(isnan.(reps.RIlib[sel_inds_slib]))
        MinRIlib = NaN
        MaxRIlib = NaN
        MedRIlib = NaN
    else
        MinRIlib = minimum(reps.RIlib; init = Inf)
        MaxRIlib = maximum(reps.RIlib; init = -Inf)
        MedRIlib = isempty(reps.RIlib) ? NaN : round(median(reps.RIlib), digits = 1)
    end

    Int_ = zeros(1, length(ind_files))
    ID_ = fill("NA", (1, length(ind_files)))
    SPEC_ = fill("NA", (1, length(ind_files)))
    FR_ = fill("NA", (1, length(ind_files)))
    CAS_ = fill("NA", (1, length(ind_files)))
    IN_ = fill("NA", (1, length(ind_files)))
    INtr_ = fill("NA", (1, length(ind_files)))
    SIM_ = fill("NA", (1, length(ind_files)))
    HIGH_ = fill("NA", (1, length(ind_files)))
    BASM_ = fill("NA", (1, length(ind_files)))
    ACTM_ = fill("NA", (1, length(ind_files)))
    FWHH_ = fill("NA", (1, length(ind_files)))
    PROB_ = fill("NA", (1, length(ind_files)))
    EXPI_ = fill("NA", (1, length(ind_files)))
    OBSI_ = fill("NA", (1, length(ind_files)))
    DRI_ = fill("NA", (1, length(ind_files)))
    CLS_ = fill("NA", (1, length(ind_files)))
    ORI_ = fill("NA", (1, length(ind_files)))
    QUAM_ = fill("NA", (1, length(ind_files)))
    
    for i =1:length(ind_files)
        # println(i)
        selected_f = reps[reps.F_ID .== i, :]

        if size(selected_f, 1) >1
            # Multiple candidates
            Int_[i] = sum(selected_f.Int)
            ID_[i] = prod(selected_f.Name.* " &&& ")[1:end-5]
            all_sp = ""
            for sp in selected_f.Spec
                sp = round.(getSpec(sp))
                sp = mergeSpec(sp)
                all_sp = all_sp * " &&& " * prod(string.(sp[:, 1]) .* ":" .* string.(sp[:, 2]) .* " ")[1:end-1]
            end
            SPEC_[i] = all_sp[6:end]
            FR_[i] = prod(selected_f.Form.* " &&& ")[1:end-5]
            CAS_[i] = prod(selected_f.CAS.* " &&& ")[1:end-5]
            IN_[i] = prod(selected_f.Inchi.* " &&& ")[1:end-5]
            INtr_[i] = prod(selected_f.InchiTrunc.* " &&& ")[1:end-5]
            SIM_[i] = prod(string.(selected_f.Sim).* " &&& ")[1:end-5]
            HIGH_[i] = prod(string.(selected_f.Height).* " &&& ")[1:end-5]
            BASM_[i] = prod(string.(selected_f.BaseMass).* " &&& ")[1:end-5]
            ACTM_[i] = prod(selected_f.ActualMass.* " &&& ")[1:end-5]
            FWHH_[i] = prod(string.(selected_f.FWHH).* " &&& ")[1:end-5]
            PROB_[i] = prod(string.(selected_f.libProb).* " &&& ")[1:end-5]
            EXPI_[i] = prod(string.(selected_f.ExpIon).* " &&& ")[1:end-5]
            OBSI_[i] = prod(string.(selected_f.ObsIon).* " &&& ")[1:end-5]
            DRI_[i] = prod(string.(selected_f.dRI).* " &&& ")[1:end-5]
            CLS_[i] = prod(string.(selected_f.Class).* " &&& ")[1:end-5]
            ORI_[i] = prod(string.(selected_f.Origin).* " &&& ")[1:end-5]
            QUAM_[i] = prod(string.(selected_f.QuantMass).* " &&& ")[1:end-5]
        elseif size(selected_f,1) == 1
            # Found only one entry
            Int_[i] = selected_f.Int[1]
            ID_[i] = selected_f.Name[1]
            sp = round.(getSpec(selected_f.Spec[1]))
            sp = mergeSpec(sp)
            SPEC_[i] = prod(string.(sp[:, 1]) .* ":" .* string.(sp[:, 2]) .* " ")[1:end-1]
            FR_[i] = selected_f.Form[1]
            CAS_[i] = selected_f.CAS[1]
            IN_[i] = selected_f.Inchi[1]
            INtr_[i] = selected_f.InchiTrunc[1]
            SIM_[i] = string(selected_f.Sim[1])
            HIGH_[i] = string(selected_f.Height[1])
            BASM_[i] = string(selected_f.BaseMass[1])
            ACTM_[i] = selected_f.ActualMass[1]
            FWHH_[i] = string(selected_f.FWHH[1])
            PROB_[i] = string(selected_f.libProb[1])
            EXPI_[i] = string(selected_f.ExpIon[1])
            OBSI_[i] = string(selected_f.ObsIon[1])
            DRI_[i] = string(selected_f.dRI[1])
            CLS_[i] = string(selected_f.Class[1])
            ORI_[i] = string(selected_f.Origin[1])
            QUAM_[i] = string(selected_f.QuantMass[1])
        end 
    end

    return(MinRt, MaxRt, AveRt, MinMass, MaxMass, AveMass, MinRI, MaxRI, MedRI, MinRIlib, MaxRIlib, MedRIlib, Int_, ID_, SPEC_, FR_, CAS_, IN_, INtr_, SIM_, HIGH_, BASM_, ACTM_, FWHH_, PROB_, EXPI_, OBSI_, DRI_, CLS_, ORI_, QUAM_)
end



# A function to merge spectra
function mergeSpec(sp)
    sp_m = zeros(length(unique(sp[:, 1])), 2)
    sp_m[:, 1] = unique(sp[:, 1])
    for m = 1:size(sp_m, 1)
        sp_m[m,2] = sum(sp[sp[:, 1] .== sp_m[m], 2])
    end

    return Int.(sp_m)
end



# A function to create aligned overview table
function finalOverviewTable(table_int, table_id, table_sp, table_fr, table_cas, table_in, table_intr, table_sim, table_high, table_prob, table_dri, table_basm, table_actm, table_cls, table_ori, table_qua, overviewMerging = "Height", weights = [1,1,1])
    # Merging can be either Height based or equal contribution
    overview = table_int[:, 1:13]
    # Highest similarity
    overview[!,"DetectionNum"] = vec(sum(Matrix(table_int[:, 14:end]) .> 0, dims = 2))
    overview[!,"DetectionFreq"] = round.(vec((sum(Matrix(table_int[:, 14:end]) .> 0, dims = 2)./ size(table_int[:, 14:end], 2)).*100))
    overview[!,"Name"] = fill("NA", size(overview, 1))
    overview[!,"Formula"] = fill("NA", size(overview, 1))
    overview[!,"InChIKey"] = fill("NA", size(overview, 1))
    # overview[!,"InChIKey_Trunc"] = fill("NA",size(overview,1))
    overview[!,"CAS"] = fill("NA", size(overview, 1))
    overview[!,"Area"] = fill(NaN, size(overview, 1))
    overview[!,"TotalArea"] = fill(NaN, size(overview, 1))
    overview[!,"QuantMasses"] = fill("NA", size(overview, 1))
    overview[!,"BaseMass"] = fill("NA", size(overview, 1))
    overview[!,"ActualMasses"] = fill("NA", size(overview, 1))
    overview[!,"Origin"] = fill("NA", size(overview, 1))
    overview[!,"Class"] = fill("NA", size(overview, 1))
    overview[!,"CAS"] = fill("NA", size(overview, 1))
    overview[!,"Similarity"] = fill(NaN, size(overview, 1))
    overview[!,"Probability"] = fill(NaN, size(overview, 1))
    overview[!,"dRI"] = fill(NaN, size(overview, 1))
    overview[!,"Score"] = fill(NaN, size(overview, 1))
    overview[!,"Spectrum_BestMatch"] = fill("NA", size(overview, 1))
    # All info
    overview[!,"All_Names"] = fill("NA", size(overview, 1))
    overview[!,"All_Formulae"] = fill("NA", size(overview, 1))
    overview[!,"All_InChIKey"] = fill("NA", size(overview, 1))
    # overview[!,"All_InChIKey_Trunc"] = fill("NA",size(overview,1))
    overview[!,"All_TotalArea"] = fill("NA", size(overview, 1))
    overview[!,"All_QuantMasses"] = fill("NA", size(overview, 1))
    overview[!,"All_BaseMass"] = fill("NA", size(overview, 1))
    overview[!,"All_ActualMasses"] = fill("NA", size(overview, 1))
    overview[!,"All_Origin"] = fill("NA", size(overview, 1))
    overview[!,"All_Class"] = fill("NA", size(overview, 1))
    overview[!,"All_CAS"] = fill("NA", size(overview, 1))
    overview[!,"All_Similarity"] = fill("NA", size(overview, 1))
    overview[!,"Spectrum_Consensus"] = fill("NA", size(overview, 1))

    for i = 1:size(overview, 1)
        println(i)
        # Fill in all info
        str = prod(unique(table_id[i, 14:end])[(unique(table_id[i, 14:end]) .!= "NA")] .* " &&& ")[1:end-5]
        overview[i,"All_Names"] = prod(split(str, " &&& ")[split(str, " &&& ") .!= "NA"] .* " &&& ")[1:end-5]
        str = prod(unique(table_fr[i, 14:end])[unique(table_fr[i, 14:end]) .!= "NA"] .* " &&& ")[1:end-5]
        overview[i,"All_Formulae"] = prod(split(str, " &&& ")[split(str, " &&& ") .!= "NA"] .* " &&& ")[1:end-5]
        str = prod(unique(table_in[i, 14:end])[(unique(table_in[i, 14:end]) .!= "NA") .& (unique(table_in[i, 14:end]) .!= "NoID")] .* " &&& ")[1:end-5]
        overview[i,"All_InChIKey"] = prod(split(str, " &&& ")[(split(str, " &&& ") .!= "NA") .& (split(str, " &&& ") .!= "NoID")] .* " &&& ")[1:end-5]
        # str = prod(unique(table_intr[i, 14:end])[(unique(table_intr[i, 14:end]) .!= "NA") .& (unique(table_intr[i, 14:end]) .!= "NoID")] .* " &&& ")[1:end-5]
        # overview[i, "All_InChIKey_Trunc"] = prod(split(str, " &&& ")[(split(str, " &&& ") .!= "NA") .& (split(str, " &&& ") .!= "NoID")] .* " &&& ")[1:end-5]
        str = prod(unique(table_cas[i, 14:end])[unique(table_cas[i, 14:end]) .!= "NA"] .* " &&& ")[1:end-5]
        overview[i,"All_CAS"] = prod(split(str, " &&& ")[split(str, " &&& ") .!= "NA"] .* " &&& ")[1:end-5]
        str = prod(unique(table_sim[i, 14:end])[unique(table_sim[i, 14:end]) .!= "NA"] .* " &&& ")[1:end-5]
        overview[i,"All_Similarity"] = prod(split(str, " &&& ")[split(str, " &&& ") .!= "NA"] .* " &&& ")[1:end-5]
        strq = prod(unique(table_qua[i, 14:end])[unique(table_qua[i, 14:end]) .!= "NA"] .* " &&& ")[1:end-5]
        overview[i,"All_QuantMasses"] = prod(split(strq, " &&& ")[split(strq, " &&& ") .!= "NA"] .* " &&& ")[1:end-5]
        str = prod(unique(table_basm[i, 14:end])[unique(table_basm[i, 14:end]) .!= "NA"] .* " &&& ")[1:end-5]
        overview[i,"All_BaseMass"] = prod(split(str, " &&& ")[split(str, " &&& ") .!= "NA"] .* " &&& ")[1:end-5]
        str = prod(unique(table_actm[i, 14:end])[unique(table_actm[i, 14:end]) .!= "NA"] .* " &&& ")[1:end-5]
        overview[i,"All_ActualMasses"] = prod(split(str, " &&& ")[split(str, " &&& ") .!= "NA"] .* " &&& ")[1:end-5]
        str = prod(unique(table_ori[i, 14:end])[unique(table_ori[i, 14:end]) .!= "NA"] .* " &&& ")[1:end-5]
        overview[i,"All_Origin"] = prod(split(str, " &&& ")[split(str, " &&& ") .!= "NA"] .* " &&& ")[1:end-5]
        str = prod(unique(table_cls[i, 14:end])[unique(table_cls[i, 14:end]) .!= "NA"] .* " &&& ")[1:end-5]
        overview[i,"All_Class"] = prod(split(str, " &&& ")[split(str, " &&& ") .!= "NA"] .* " &&& ")[1:end-5]
        
        # Consensus spectrum
        sps = table_sp[i, 14:end]
        allspec = [0 0]
        for s = 1:length(sps)
            if sps[s] .== "NA"
                continue
            end
            sun = split(sps[s], " &&& ")
            for n = 1:length(s)
                spec = getSpec(sun[n])
                if lowercase(overviewMerging) == "equal"
                    allspec = [allspec; spec]
                elseif lowercase(overviewMerging) == "height"
                    height = parse.(Float64, split(table_high[i, 13+s], " &&& ")[n])
                    spec[:, 2] = spec[:, 2] .* height
                    allspec = [allspec; spec]
                else
                    println("Choose a valid spectra merging weight option: \"height\" or \"equal\"")
                end
            end
        end
        # Merge signals
        allspec = mergeSpec(allspec[2:end, :])
        # Normalize
        allspec[:, 2] = Int.(round.((allspec[:, 2]./maximum(allspec[:, 2])) .* 1000))
        # Remove with intensity of 0
        allspec = allspec[allspec[:, 2] .!= 0, :]
        overview[i,"Spectrum_Consensus"] = prod(string.(allspec[:, 1]) .* ":" .* string.(allspec[:, 2]) .* " ")[1:end-1]

        # Higest similarity match
        indc, inde, val, prob_, dri_, score = findMaxSim(table_sim[i, 14:end], table_prob[i, 14:end], table_dri[i, 14:end], weights)
        if score == -100
            # No identification present
            continue
        end
        # Highest similarity information
        if !contains(table_id[i, 13+indc], " &&& ")
            overview[i, "Name"] = table_id[i, 13+indc]
            overview[i, "Formula"] = table_fr[i, 13+indc]
            overview[i, "InChIKey"] = table_in[i, 13+indc]
            # overview[i, "InChIKey_Trunc"] = table_intr[i, 13+indc]
            overview[i, "CAS"] = table_cas[i, 13+indc]
            overview[i, "Area"] = table_int[i, 13+indc]
            overview[i, "Similarity"] = val
            overview[i, "Probability"] = prob_
            overview[i, "dRI"] = dri_
            overview[i, "Score"] = score
            overview[i, "Spectrum_BestMatch"] = table_sp[i, 13+indc]
            overview[i, "QuantMasses"] = table_qua[i, 13+indc]
            overview[i, "BaseMass"] = table_basm[i, 13+indc]
            overview[i, "ActualMasses"] = table_actm[i, 13+indc]
            overview[i, "Origin"] = table_ori[i, 13+indc]
            overview[i, "Class"] = table_cls[i, 13+indc]
        else
            overview[i, "Name"] = split(table_id[i, 13+indc], " &&& ")[inde]
            overview[i, "Formula"] = split(table_fr[i, 13+indc], " &&& ")[inde]
            overview[i, "InChIKey"] = split(table_in[i, 13+indc], " &&& ")[inde]
            # overview[i, "InChIKey_Trunc"] = split(table_intr[i, 13+indc], " &&& ")[inde]
            overview[i,"CAS"] = split(table_cas[i, 13+indc], " &&& ")[inde]
            overview[i, "Area"] = table_int[i, 13+indc]
            overview[i, "Similarity"] = val
            overview[i, "Probability"] = prob_
            overview[i, "dRI"] = dri_
            overview[i, "Score"] = score
            overview[i, "Spectrum_BestMatch"] = split(table_sp[i, 13+indc], " &&& ")[inde]
            overview[i, "QuantMasses"] = split(table_qua[i, 13+indc], " &&& ")[inde]
            overview[i, "BaseMass"] = split(table_basm[i, 13+indc], " &&& ")[inde]
            overview[i, "ActualMasses"] = split(table_actm[i, 13+indc], " &&& ")[inde]
            overview[i, "Origin"] = split(table_ori[i, 13+indc], " &&& ")[inde]
            overview[i, "Class"] = split(table_cls[i, 13+indc], " &&& ")[inde]
        end

        # Calculate all total areas
        for j = 14:size(table_int, 2)
            if table_int[i, j] .!= 0
                table_int[i, j] = calcTotalArea([table_int[i, j]], [table_qua[i, j]], [overview[i, "Spectrum_Consensus"]])[1]
            end
        end
    end
    # Caculate total component areas
    overview[:, "TotalArea"] = calcTotalArea(overview[!, "Area"], overview[!, "QuantMasses"], overview[!, "Spectrum_BestMatch"])
    # Add table_int
    overview = [overview table_int[:, 14:end]]

    return overview
end



# A function to retrieve the feature with max. similarity
function findMaxSim(sim, prob, dri, w)
    indc = 1
    inde = 1
    val = 0
    prob_ = 0
    dri_ = 0
    score = -100
    for n = 1:length(sim)
        # println("n $n")
        if sim[n] .!= "NA"
            if contains(sim[n], " &&& ")
                valc = parse.(Float64, split(sim[n], " &&& ")) ./ 1000
                valc[valc .< 0] .= 0
                probc = parse.(Float64, split(prob[n], " &&& ")) ./ 100
                probc[isnan.(probc)] .= 0
                dric = abs.(parse.(Float64, split(dri[n], " &&& ")))
                # Low RI diff ==1 high RI (>=250) diff == 0
                dric = 1 .- (dric ./ 250)
                dric[dric .< 0] .= 0
                dric[isnan.(dric)] .= 0
                # Calculate score
                sc = (valc .* w[1] .+ probc .* w[2] .+ dric .* w[3]) / sum(w)
                if any(sc .> score)
                    indc = n
                    inde = argmax(sc)
                    val = valc[inde]
                    prob_ = probc[inde]
                    dri_ = dric[inde]
                    score = sc[inde]
                end
            else
                valc = parse(Float64, sim[n]) / 1000
                if valc < 0
                    valc = 0
                end
                probc = parse(Float64, prob[n]) / 100
                if isnan(probc)
                    probc = 0
                end
                dric = abs(parse(Float64, dri[n]))
                # Low RI diff ==1 high RI (>=250) diff == 0
                dric = 1 - (dric / 250)
                if (dric < 0) || isnan(dric)
                    dric = 0
                end
                # Calculate score
                sc = (valc .* w[1] .+ probc .* w[2] .+ dric .* w[3]) / sum(w)
                if sc > score
                    val = valc
                    prob_ = probc
                    dri_ = dric
                    score = sc
                    indc = n
                    inde = 1
                end                
            end
        end
    end
    # Calculate the adjusted area

    return indc, inde, val, prob_, dri_, score
end



# A function to align RT1 fractions (parts)
function feature_align_fraction(rep, name, ind_files, rt_tol, rt2_tol, rt2_tol_lowA, rt2_tol_highA, rt_tol_lowA,
    h_thresh, rt_tol_highA, sim_thresh, mz_tol, highTol, feature, Inten, Idents, Specs, Forms, Cas, Inchi, InchiT, 
    Sim, High, BasM, ActM, Fwhh, Prob, ExpI, ObsI, dRI, Cls, Ori, QuaM, maxRTmatch; numfrags::Int = 150,
    similarity_method::String = "DISCO")
    outside = sum(rep.Rt .> maxRTmatch)
    i = 1
    while i <= (size(rep,1) - outside)
        println("Processing feature $i of $(size(rep, 1) - outside)")
        
        if i > size(rep, 1)  # safety if rep shrinks
            break
        end

        if rep.Rt[i] == 0 || rep.Rt[i] > maxRTmatch
            i += 1
            continue
        end

        if highTol && (rep.Height[i] < h_thresh)
            rt_tol = rt2_tol_lowA
            rt2_tol = rt2_tol_lowA
            println("Low peak height")
            highTol = false
        end

        ind = i
        sel_inds = select_candidates_rt1_(rep, ind, rt_tol)
        sel_inds_f = select_candidates_rt2_(rep, ind, sel_inds, rt2_tol)
        sel_inds_s = select_candidates_dot_(rep, ind, sel_inds_f, sim_thresh, mz_tol; numfrags = numfrags, similarity_method = similarity_method)

        if !isempty(sel_inds_s)
            MinRt, MaxRt, AveRt, MinMass, MaxMass, AveMass, MinRI, MaxRI, MedRI, MinRIlib, MaxRIlib, MedRIlib, Int_, ID_, SPEC_, FR_, CAS_, IN_, INtr_, SIM_, HIGH_, BASM_, ACTM_, FWHH_, PROB_, EXPI_, OBSI_, DRI_, CLS_, ORI_, QUAM_= group_candidates_external(rep[sel_inds_s, :], ind_files)
            
            rep.Rt[sel_inds_s] .= 0
            
            feature = vcat(feature, hcat(MinRt, MaxRt, AveRt, MinMass, MaxMass, AveMass, MinRI, MaxRI, MedRI, MinRIlib, MaxRIlib, MedRIlib))
            Inten = vcat(Inten, Int_)
            Idents = vcat(Idents, ID_)
            Specs = vcat(Specs, SPEC_)
            Forms = vcat(Forms, FR_)
            Cas = vcat(Cas, CAS_)
            Inchi = vcat(Inchi, IN_)
            InchiT = vcat(InchiT, INtr_)
            Sim = vcat(Sim, SIM_)
            High = vcat(High, HIGH_)
            BasM = vcat(BasM, BASM_)
            ActM = vcat(ActM, ACTM_)
            Fwhh = vcat(Fwhh, FWHH_)
            Prob = vcat(Prob, PROB_)
            ExpI = vcat(ExpI, EXPI_)
            ObsI = vcat(ObsI, OBSI_)
            dRI = vcat(dRI, DRI_)
            Cls = vcat(Cls, CLS_)
            Ori = vcat(Ori, ORI_)
            QuaM = vcat(QuaM, QUAM_)
        end

        if mod(i, 50) == 0
            # Reduce rep
            prev_n = size(rep, 1)
            rep = rep[rep.Rt .!= 0, :]
            outside = sum(rep.Rt .> maxRTmatch)  # keep loop bound in sync
            if size(rep, 1) < prev_n
                println("Part reduced to $(size(rep, 1))")
                i = 0  # restart scan only when something actually changed
            end
        end
        i += 1
    end

    # Get ungrouped features outside of time domain
    rep_rem = rep

    return feature, Inten, Idents, Specs, Forms, Cas, Inchi, InchiT, Sim, High, BasM, ActM, Fwhh, Prob, ExpI, ObsI, dRI, Cls, Ori, QuaM, rep_rem
end



# A general alignment function
function feature_align_external(path2files, rt2_tol_lowA, rt2_tol_highA, rt_tol_lowA, h_thresh, rt_tol_highA, sim_thresh, mz_tol, overviewMerging = "Height", weights = [1,1,1], tr1_parts = 5; numfrags::Int = 150, similarity_method::String = "DISCO")
    rep,name = report_import_external(path2files)
    time_b4_align = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Aligning ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " samples... It will take a while... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "\n")

    if length(name) == 0
        @warn "No file has been imported."
        return res = []
    end 

    feature = zeros(1, 12) # Initial matrix for saving/appending RT and RI values
    Inten = zeros(1, length(name))
    Idents = fill("NA", (1, length(name)))
    Specs = fill("NA", (1, length(name)))
    Forms = fill("NA", (1, length(name)))
    Cas = fill("NA", (1, length(name)))
    Inchi = fill("NA", (1, length(name)))
    InchiT = fill("NA", (1, length(name)))
    Sim = fill("NA", (1, length(name)))
    High = fill("NA", (1, length(name)))
    BasM = fill("NA", (1, length(name)))
    ActM = fill("NA", (1, length(name)))
    Fwhh = fill("NA", (1, length(name)))
    Prob = fill("NA", (1, length(name)))
    ExpI = fill("NA", (1, length(name)))
    ObsI = fill("NA", (1, length(name)))
    dRI = fill("NA", (1, length(name)))
    Cls = fill("NA", (1, length(name)))
    Ori = fill("NA", (1, length(name)))
    QuaM = fill("NA", (1, length(name)))

    # Sort list from xx to xx intensity
    rep = rep[sortperm(rep[!, "Height"], rev = true), :]
    rt_tol = rt_tol_highA
    rt2_tol = rt_tol_highA
    highTol = true

    ind_files = unique(rep.F_ID)

    # Selection of parts based on equal data distribution
    tr1_part_lim = zeros(tr1_parts+1)
    tr1_part_lim[1] = minimum(rep[!, "Rt"])
    tr1_part_lim[end] = maximum(rep[!, "Rt"])
    rtps = sort(rep[!, "Rt"])
    for p = 2:tr1_parts
        tr1_part_lim[p] = rtps[Int(round((length(rtps) / tr1_parts) * (p-1)))] # Select retention time at specific part of data
    end
    # tr1_part_lim = collect(minimum(rep[!, "Rt"]):(maximum(rep[!, "Rt"]) - minimum(rep[!, "Rt"])) / tr1_parts:maximum(rep[!, "Rt"]))

    for t = 1:length(tr1_part_lim)-1
        ind_part = findall((tr1_part_lim[1] - 3*rt_tol) .<= rep[!, "Rt"] .<= (tr1_part_lim[t+1] + 1*rt_tol))
        println("Part $(t)/$tr1_parts: passing $(length(ind_part)) features")
        maxRTmatch = tr1_part_lim[t+1]
        
        feature, Inten, Idents, Specs, Forms, Cas, Inchi, InchiT, Sim, High, BasM, ActM, Fwhh, Prob, ExpI, ObsI, 
        dRI, Cls, Ori, QuaM, rep_rem = feature_align_fraction(rep[ind_part, :], name, ind_files, rt_tol, 
        rt2_tol, rt2_tol_lowA, rt2_tol_highA, rt_tol_lowA, h_thresh, rt_tol_highA, sim_thresh, mz_tol, highTol, 
        feature, Inten, Idents, Specs, Forms, Cas, Inchi, InchiT, Sim, High, BasM, ActM, Fwhh, Prob, ExpI, ObsI, 
        dRI, Cls, Ori, QuaM, maxRTmatch; numfrags = numfrags, similarity_method = similarity_method)
        
        rep = [rep[.!((tr1_part_lim[1] - 3*rt_tol) .<= rep[!, "Rt"] .<= (tr1_part_lim[t+1] + 1*rt_tol)), :]; rep_rem]
        println("Features remaining in rep $(size(rep, 1))")
    end

    # Aligning remaining features
    i = 1
    println("Final part: passing $(size(rep, 1)) features")

    while i <= size(rep, 1)
    # while i = 1:size(rep,1)
        println("Processing feature $i of $(size(rep, 1))")

        if i > size(rep, 1)
            break
        end        
        
        if rep.Rt[i] == 0
            i+=1
            continue
        end

        if highTol && (rep.Height[i] < h_thresh)
            rt_tol = rt_tol_lowA
            rt2_tol = rt2_tol_lowA
            println("Low peak height")
            highTol = false
        end

        ind = i
        # println(i)
        sel_inds = select_candidates_rt1_(rep, ind, rt_tol)
        sel_inds_f = select_candidates_rt2_(rep, ind, sel_inds, rt2_tol)
        sel_inds_s = select_candidates_dot_(rep, ind, sel_inds_f, sim_thresh, mz_tol; numfrags = numfrags, similarity_method = similarity_method)
        # println(length(sel_inds_s))

        if !isempty(sel_inds_s)
            MinRt, MaxRt, AveRt, MinMass, MaxMass, AveMass, MinRI, MaxRI, MedRI, MinRIlib, MaxRIlib, MedRIlib, Int_, ID_, SPEC_, FR_, CAS_, IN_, INtr_, SIM_, HIGH_, BASM_, ACTM_, FWHH_, PROB_, EXPI_, OBSI_, DRI_, CLS_, ORI_, QUAM_= group_candidates_external(rep[sel_inds_s, :], ind_files)
            
            rep.Rt[sel_inds_s] .= 0

            feature = vcat(feature, hcat(MinRt, MaxRt, AveRt, MinMass, MaxMass, AveMass, MinRI, MaxRI, MedRI, MinRIlib, MaxRIlib, MedRIlib))
            Inten = vcat(Inten, Int_)
            Idents = vcat(Idents, ID_)
            Specs = vcat(Specs, SPEC_)
            Forms = vcat(Forms, FR_)
            Cas = vcat(Cas, CAS_)
            Inchi = vcat(Inchi, IN_)
            InchiT = vcat(InchiT, INtr_)
            Sim = vcat(Sim, SIM_)
            High = vcat(High, HIGH_)
            BasM = vcat(BasM, BASM_)
            ActM = vcat(ActM, ACTM_)
            Fwhh = vcat(Fwhh, FWHH_)
            Prob = vcat(Prob, PROB_)
            ExpI = vcat(ExpI, EXPI_)
            ObsI = vcat(ObsI, OBSI_)
            dRI = vcat(dRI, DRI_)
            Cls = vcat(Cls, CLS_)
            Ori = vcat(Ori, ORI_)
            QuaM = vcat(QuaM, QUAM_)
        end

        if mod(i, 10) == 0
            # Reduce rep
            prev_n = size(rep, 1)
            rep = rep[rep.Rt .!= 0, :]
            if size(rep, 1) < prev_n
                println("Reduced to $(size(rep, 1))")
                i = 0  # restart only if something changed
            end
        end
        i += 1
    end

    table1=DataFrame(hcat(1:size(feature, 1)-1, feature[2:end, :]), [:Nr, :MinRt1D, :MaxRt1D, :AveRt1D,
    :MinRt2D, :MaxRt2D, :AveRt2D, :MinRI, :MaxRI, :MedRI, :MinRIlib, :MaxRIlib, :MedRIlib])

    table_int=DataFrame(Inten[2:end, :], :auto)
    table_id=DataFrame(Idents[2:end, :], :auto)
    table_sp=DataFrame(Specs[2:end, :], :auto)
    table_fr=DataFrame(Forms[2:end, :], :auto)
    table_cas=DataFrame(Cas[2:end, :], :auto)
    table_in=DataFrame(Inchi[2:end, :], :auto)
    table_intr=DataFrame(InchiT[2:end, :], :auto)
    table_sim=DataFrame(Sim[2:end, :], :auto)
    table_high=DataFrame(High[2:end, :], :auto)
    table_basm=DataFrame(BasM[2:end, :], :auto)
    table_actm=DataFrame(ActM[2:end, :], :auto)
    table_fwhh=DataFrame(Fwhh[2:end, :], :auto)
    table_prob=DataFrame(Prob[2:end, :], :auto)
    table_expi=DataFrame(ExpI[2:end, :], :auto)
    table_obsi=DataFrame(ObsI[2:end, :], :auto)
    table_dri=DataFrame(dRI[2:end, :], :auto)
    table_cls=DataFrame(Cls[2:end, :], :auto)
    table_ori=DataFrame(Ori[2:end, :], :auto)
    table_qua=DataFrame(QuaM[2:end, :], :auto)
   
    rename!(table_int, Symbol.(name))
    rename!(table_id, Symbol.(name))
    rename!(table_sp, Symbol.(name))
    rename!(table_fr, Symbol.(name))
    rename!(table_cas, Symbol.(name))
    rename!(table_in, Symbol.(name))
    rename!(table_intr, Symbol.(name))
    rename!(table_sim, Symbol.(name))
    rename!(table_high, Symbol.(name))
    rename!(table_basm, Symbol.(name))
    rename!(table_actm, Symbol.(name))
    rename!(table_fwhh, Symbol.(name))
    rename!(table_prob, Symbol.(name))
    rename!(table_expi, Symbol.(name))
    rename!(table_obsi, Symbol.(name))
    rename!(table_dri, Symbol.(name))
    rename!(table_cls, Symbol.(name))
    rename!(table_ori, Symbol.(name))
    rename!(table_qua, Symbol.(name))
    
    table_int=hcat(table1, table_int)
    table_id=hcat(table1, table_id)
    table_sp=hcat(table1, table_sp)
    table_fr=hcat(table1, table_fr)
    table_cas=hcat(table1, table_cas)
    table_in=hcat(table1, table_in)
    table_intr=hcat(table1, table_intr)
    table_sim=hcat(table1, table_sim)
    table_high=hcat(table1, table_high)
    table_basm=hcat(table1, table_basm)
    table_actm=hcat(table1, table_actm)
    table_fwhh=hcat(table1, table_fwhh)
    table_prob=hcat(table1, table_prob)
    table_expi=hcat(table1, table_expi)
    table_obsi=hcat(table1, table_obsi)
    table_dri=hcat(table1, table_dri)
    table_cls=hcat(table1, table_cls)
    table_ori=hcat(table1, table_ori)
    table_qua=hcat(table1, table_qua)

    sort!(table_int, [:AveRt1D, :AveRt2D])
    sort!(table_id, [:AveRt1D, :AveRt2D])
    sort!(table_sp, [:AveRt1D, :AveRt2D])
    sort!(table_fr, [:AveRt1D, :AveRt2D])
    sort!(table_cas, [:AveRt1D, :AveRt2D])
    sort!(table_in, [:AveRt1D, :AveRt2D])
    sort!(table_intr, [:AveRt1D, :AveRt2D])
    sort!(table_sim, [:AveRt1D, :AveRt2D])
    sort!(table_high, [:AveRt1D, :AveRt2D])
    sort!(table_basm, [:AveRt1D, :AveRt2D])
    sort!(table_actm, [:AveRt1D, :AveRt2D])
    sort!(table_fwhh, [:AveRt1D, :AveRt2D])
    sort!(table_prob, [:AveRt1D, :AveRt2D])
    sort!(table_expi, [:AveRt1D, :AveRt2D])
    sort!(table_obsi, [:AveRt1D, :AveRt2D])
    sort!(table_dri, [:AveRt1D, :AveRt2D])
    sort!(table_cls, [:AveRt1D, :AveRt2D])
    sort!(table_ori, [:AveRt1D, :AveRt2D])
    sort!(table_qua, [:AveRt1D, :AveRt2D])

    table_int.Nr = 1:size(feature, 1)-1
    table_id.Nr = 1:size(feature, 1)-1
    table_sp.Nr = 1:size(feature, 1)-1
    table_fr.Nr = 1:size(feature, 1)-1
    table_cas.Nr = 1:size(feature, 1)-1
    table_in.Nr = 1:size(feature, 1)-1
    table_intr.Nr = 1:size(feature, 1)-1
    table_high.Nr = 1:size(feature, 1)-1
    table_basm.Nr = 1:size(feature, 1)-1
    table_actm.Nr = 1:size(feature, 1)-1
    table_fwhh.Nr = 1:size(feature, 1)-1
    table_prob.Nr = 1:size(feature, 1)-1
    table_expi.Nr = 1:size(feature, 1)-1
    table_obsi.Nr = 1:size(feature, 1)-1
    table_dri.Nr = 1:size(feature, 1)-1
    table_cls.Nr = 1:size(feature, 1)-1
    table_ori.Nr = 1:size(feature, 1)-1
    table_qua.Nr = 1:size(feature, 1)-1

    datetime = Dates.format(Dates.now(), "yymmdd")

    output = joinpath(path2files, datetime*"_Aligned_Area.csv")
    CSV.write(output, table_int)
    output = joinpath(path2files, datetime*"_Aligned_ID.csv")
    CSV.write(output, table_id)
    output = joinpath(path2files, datetime*"_Aligned_Spec.csv")
    CSV.write(output, table_sp)
    output = joinpath(path2files, datetime*"_Aligned_Formula.csv")
    CSV.write(output, table_fr)
    output = joinpath(path2files, datetime*"_Aligned_CAS.csv")
    CSV.write(output, table_cas)
    output = joinpath(path2files, datetime*"_Aligned_InChIKey.csv")
    CSV.write(output, table_in)
    output = joinpath(path2files, datetime*"_Aligned_InChIKeyTrunc.csv")
    CSV.write(output, table_intr)
    output = joinpath(path2files, datetime*"_Aligned_Similarity.csv")
    CSV.write(output, table_sim)
    output = joinpath(path2files, datetime*"_Aligned_Height.csv")
    CSV.write(output, table_high)
    output = joinpath(path2files, datetime*"_Aligned_BaseMass.csv")
    CSV.write(output, table_basm)
    output = joinpath(path2files, datetime*"_Aligned_ActualMass.csv")
    CSV.write(output, table_actm)
    output = joinpath(path2files, datetime*"_Aligned_FWHH.csv")
    CSV.write(output, table_fwhh)
    output = joinpath(path2files, datetime*"_Aligned_Probability.csv")
    CSV.write(output, table_prob)
    output = joinpath(path2files, datetime*"_Aligned_ExperimentalIon.csv")
    CSV.write(output, table_expi)
    output = joinpath(path2files, datetime*"_Aligned_ObservedIon.csv")
    CSV.write(output, table_obsi)
    output = joinpath(path2files, datetime*"_Aligned_DeltaRI.csv")
    CSV.write(output, table_dri)
    output = joinpath(path2files, datetime*"_Aligned_Classification.csv")
    CSV.write(output, table_cls)
    output = joinpath(path2files, datetime*"_Aligned_Origin.csv")
    CSV.write(output, table_ori)
    output = joinpath(path2files, datetime*"_Aligned_QuantMasses.csv")
    CSV.write(output, table_qua)

    # table_basm, table_actm, table_fwhh, table_expi, table_obsi, table_cls, table_ori, table_qua
    overview = finalOverviewTable(table_int, table_id, table_sp, table_fr, table_cas, table_in, table_intr, table_sim, table_high, table_prob, table_dri,
                                  table_basm, table_actm, table_cls, table_ori, table_qua, overviewMerging, weights)

    output=joinpath(path2files, datetime*"_Aligned.csv")
    CSV.write(output, overview)

    println("\n", "Done aligning ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .csv files.")
    println("Start time: ", time_b4_align)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)

    return
end



# A function to merge features of the same origin in the sample, e.g., the main peak & tailing
function internal_feature_merging(path2files, tr1_dev, tr2_dev, sim_thresh, mz_tol)
    nn = readdir(path2files)
    time_b4_merge = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Merging features of the same chemical origin... It will take a while... ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"), "\n")

    for i = 1:size(nn, 1)
        file = split(nn[i, :][1], ".")
        if length(file[1]) == 0 || file[end] != "csv" || any(contains.(file, "Aligned")) || any(contains.(file, "Original")) || any(contains.(file, "REMOVED")) || any(contains.(file, "merged")) || any(contains.(file, "IntMerged")) || any(contains.(file, "_IDresults")) || any(contains.(file, "_IDbest")) || any(contains.(file, "_IDconsens"))
            continue
        end 
        f_n = joinpath(path2files, nn[i, :][1])
        println("$i : " * f_n)
        df = CSV.read(f_n, DataFrame)

        # Dataframe for tracking what is being merged
        df_m = DataFrame()

        # Copy original version
        if !isfile(f_n[1:end-4] * "_Original.csv")
            cp(f_n, f_n[1:end-4] * "_Original.csv")
        end

        # Remove missin area entries
        df = df[ismissing.(df[!, "Area"]) .== 0, :]

        # Sort from high to low
        df = df[sortperm(df[!, "Height"], rev = true), :]

        # Extract retention times
        tr1, tr2 = extractRTs(df)

        mNum = 1

        for c = 1:size(df, 1)
            if tr1 == -Inf
                continue
            end
            # println(c)

            # Get entries RT1
            ind = findall((tr1[c]-tr1_dev[1] .<= tr1 .<= tr1[c]+tr1_dev[2]) .& (tr1 .!= -Inf))
            if (length(ind) == 1) && (ind[1] == c)
                # Only parent entry remains
                continue
            end

            # Get entries RT2
            ind = ind[(tr2[c]-tr2_dev[1] .<= tr2[ind] .<= tr2[c]+tr2_dev[2]) .& (tr2[ind] .!= -Inf)]
            if (length(ind) == 1) && (ind[1] == c)
                # Only parent entry remains
                continue
            end

            # Remove parent from ind
            ind = ind[ind .!= c]

            # Check spectral similarity
            vals = zeros(length(ind))
            numfrags = 150 # default number of fragments to use
            pfs = getSpec(df[c, "Spectrum"])
            for v = 1:length(ind)
                spec = getSpec(df[ind[v], "Spectrum"])
                master_vect_f = master_mass_gen(pfs[:, 1], pfs[:, 2], spec[:, 1], spec[:, 2], mz_tol)
                if size(master_vect_f, 1) .> numfrags
                    master_vect_fu = master_vect_f[sortperm(vec(sum(master_vect_f[:, 3], dims = 2)), rev = true), :][1:numfrags, :]
                    master_vect_fr = master_vect_f[sortperm(vec(sum(master_vect_f[:, 2], dims = 2)), rev = true), :][1:numfrags, :]
                    master_vect_f = [master_vect_fu; master_vect_fr]
                    master_vect_f = unique(master_vect_f, dims = 1)
                end
                # Mean center each vector
                master_vect_f[:, 2:3] = master_vect_f[:, 2:3] .- mean(master_vect_f[:, 2:3], dims = 1)
                vals[v] = dot(master_vect_f[:, 2], master_vect_f[:, 3]) / (norm(master_vect_f[:, 2] * norm(master_vect_f[:, 3])))
            end

            if any(vals .>= sim_thresh)
                # Group inter matches
                ind = ind[vals .>= sim_thresh]

                # Merge to parent (highest intensity spectrum) + sum area
                df[c, "Area"] = df[c, "Area"] + sum(df[ind, "Area"])

                # Save what was merged for data review
                if isempty(df_m)
                    df_m = [DataFrame(pairNum = mNum .* ones(size(df[[c; ind], :], 1))) df[[c; ind], :]]
                    mNum += 1
                else
                    df_m = append!(df_m, [DataFrame(pairNum = mNum .* ones(size(df[[c; ind], :], 1))) df[[c; ind], :]])
                    mNum += 1
                end

                # Set reference RT to weird values
                tr1[ind] .= -Inf
                tr2[ind] .= -Inf
            end
        end
        # Save df
        CSV.write(f_n, df[tr1 .!= -Inf, :])
        # Save merged pairs
        CSV.write(f_n[1:end-4] * "_IntMerged.csv", df_m)
    end

    println("\n", "Done internal feature merging for ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .csv files.")
    println("Start time: ", time_b4_merge)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to classify the features based on thier base mass and RT1/RT2 region
function filter_signal_external(path2files, filtType, base, mz, int, mz_tol, tr1_range, tr2_range)
    nn = readdir(path2files)
    time_b4_org = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Classifying features as bleed, toluene, DCM, alkanes, and acids... ", "\n")

    # Check if there is a RT range
    tr1_check = true
    if tr1_range == [0, 0]
        tr1_check = false
    end

    tr2_check = true
    if tr2_range == [0, 0]
        tr2_check = false
    end

    # Check if base is in correct direction
    lo = size(base)
    if lo[1] > 1
        base = base'
    end

    for i = 1:size(nn,1)
        file = split(nn[i, :][1], ".")
        if length(file[1]) == 0 || file[end] != "csv" || any(contains.(file, "_Aligned")) || any(contains.(file, "Original")) || any(contains.(file, "REMOVED")) || any(contains.(file, "merged")) || any(contains.(file, "IntMerged")) || any(contains.(file, "_IDresults")) || any(contains.(file, "_IDbest")) || any(contains.(file, "_IDconsens"))
            continue
        end 
        f_n = joinpath(path2files, nn[i, :][1])
        println("$i : " * f_n)
        df = CSV.read(f_n, DataFrame)
        if !any(names(df) .== "Origin")
            df[!, "Origin"] = fill("NA", size(df, 1))
        end

        # Copy original version
        if !isfile(f_n[1:end-4] * "_Original.csv")
            cp(f_n, f_n[1:end-4] * "_Original.csv")
        end

        # Get RTs
        if tr1_check || tr2_check
            rt1, rt2 = extractRTs(df)
        end

        # Filter file
        ind = findall(vec(any(abs.(round.(df[!, "Base Mass"]).-base) .< mz_tol, dims = 2)))
        for j = 1:length(ind)
            # Check if index fall within specified range
            if tr1_check && tr2_check
                if !((tr2_range[1] <= rt2[ind[j]] <= tr2_range[2]) && (tr1_range[1] <= rt1[ind[j]] <= tr1_range[2]))
                    continue
                end
            elseif tr1_check
                if !(tr1_range[1] <= rt1[ind[j]] <= tr1_range[2])
                    continue
                end
            elseif tr2_check
                if !(tr2_range[1] <= rt2[ind[j]] <= tr2_range[2])
                    continue
                end
            end

            # Get spectrum
            spec = round.(getSpec(df[ind[j], "Spectrum"]))

            # Determine reference intensity
            if any(int .== 1)
                # Use m/z with int 1 as reference
                mz_f = mz[findfirst(int .== 1)]
                if any(abs.(spec[:, 1] .- mz_f) .<= mz_tol)
                    refmz = maximum(spec[findall(abs.(spec[:, 1] .- mz_f) .<= mz_tol), 2])
                else
                    # This fragment is already missing -> go to next index
                    continue
                end
            else
                # Use base peak intensity
                refmz = maximum(spec[findall(vec(any(abs.(spec[:, 1] .- base) .<= mz_tol, dims = 2))), 2])
            end

            # For each candidate check spectrum for m/z:relative intensity pairs
            c = 0   # Counter for true reference m/z values
            for m = 1:length(mz)
                # Check presence
                if any(abs.(spec[:, 1] .- mz[m]) .<= mz_tol)
                    # If present check intensity
                    if maximum(spec[findall(abs.(spec[:, 1] .- mz[m]) .<= mz_tol), 2]) / refmz .>= int[m]
                        c += 1
                    else
                        # In case a fragment was not found at specified intensity it stays (i.e. goes to check the next index)
                        break
                    end
                else
                    # In case a fragment was not found it stays (i.e. goes to check the next index)
                    break
                end
            end
            # Check if all m/z's were correct
            if c >= length(mz)
                # If not set base mass to the filterType (i.e. entry will be removed)
                # println("Removed $(nn[i]):")      # lines to print what was removed
                # println(df[ind[j], :])
                df[ind[j], "Origin"] = filtType
            end
        end
        
        # Save file
        CSV.write(f_n, df)
    end

    println("\n", "Done classifying features as bleed, toluene, DCM, alkanes, and acids in ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4_org)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to reset previous filters (classification) if any
function reset_filters(path2files)
    nn = readdir(path2files)
    indf = findall(contains.(nn, "Original"))
    time_b4_reset = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Resetting previous filters (classification) if any... ")

    for f in indf
        # Remove filtered file
        if isfile(joinpath(path2files, nn[f])[1:end-13] * ".csv")
            rm(joinpath(path2files, nn[f])[1:end-13] * ".csv")
        end

        # Remove filtered compound file
        # if isfile(joinpath(path2files, nn[f])[1:end-13] * "_REMOVED.csv")
        #     rm(joinpath(path2files, nn[f])[1:end-13] * "_REMOVED.csv")
        # end

        # Remove merged features file
        if isfile(joinpath(path2files, nn[f])[1:end-13] * "_merged.csv")
            rm(joinpath(path2files, nn[f])[1:end-13] * "_merged.csv")
        end

        # mv original
        mv(joinpath(path2files, nn[f]), joinpath(path2files, nn[f])[1:end-13] * ".csv")
    end

    println("Done resetting ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4_reset)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to convert input TXT files to CSV format
function convert_txt2csv(path2files)
    nn = readdir(path2files)
    time_b4_conv = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Converting files to .csv format... ")

    for f in nn
        if contains.(f, ".txt")
            tmp = CSV.read(joinpath(path2files,f), DataFrame)
            CSV.write(joinpath(path2files, f)[1:end-3] * "csv",tmp)
        end
    end
    println("Done converting ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .txt files. ")
    println("Start time: ", time_b4_conv)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to remove m/z:intensity pairs with intensity values containing '0.' from spectra
function remove_zero_intensity_pairs(path2files)
    nn = readdir(path2files)
    time_b4_0 = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Removing m/z:intensity pairs with intensity values containing '0.' from spectra... ")

    # Find and process each CSV file in the directory
    for f in nn
        if contains.(f, ".csv")
            df = CSV.read(joinpath(path2files, f), DataFrame)

            # Check if 'Spectrum' column exists
            if "Spectrum" in names(df)
                # Replace missing values in 'Spectrum' with "NA"
                df.Spectrum .= coalesce.(df.Spectrum, "NA")

                # Remove m/z:intensity pairs with intensity values containing '0.' from spectra
                df.Spectrum = map(spectrum -> begin
                    pairs = split(spectrum, ' ')
                    filtered_pairs = filter(pair -> !occursin(":0.", pair), pairs)
                    # println("Filtered Pairs 0 Int m/z: \n", filtered_pairs)
                    
                    join(filtered_pairs, ' ')
                end, df.Spectrum)

                # Save the filtered DataFrame back as a CSV file
                CSV.write(joinpath(path2files, f), df)
            else
                println("Warning: 'Spectrum' column not found in file $(f). Skipping filtration.")
            end
        end
    end

    println("Done removing m/z:intensity pairs with intensity values containing '0.' from spectra of ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4_0)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to remove m/z:intensity pairs with intensity values <= 3% of the highest intensity values for each spectrum
function remove_low_intensity_pairs(path2files)
    nn = readdir(path2files)
    time_b4_low = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Removing m/z:intensity pairs with intensity values <= 3% of the highest intensity values from spectra... ")

    # Find and process each CSV file in the directory
    for f in nn
        if contains(f, ".csv")
            df = CSV.read(joinpath(path2files, f), DataFrame)

            # Check if 'Spectrum' column exists
            if "Spectrum" in names(df)
                # Replace missing values in 'Spectrum' with "NA"
                df.Spectrum .= coalesce.(df.Spectrum, "NA")

                # Remove m/z:intensity pairs with intensity values <= 3% of the highest intensity values from spectra
                df.Spectrum = map(spectrum -> begin
                    pairs = split(spectrum, ' ')
                    intensities = [parse(Float64, split(pair, ':')[2]) for pair in pairs]
                    max_intensity = maximum(intensities)
                    threshold = 0.03 * max_intensity
                    filtered_pairs = filter(pair -> parse(Float64, split(pair, ':')[2]) > threshold, pairs)
                    # println("Filtered Pairs Low Int m/z: \n", filtered_pairs)
                    
                    join(filtered_pairs, ' ')
                end, df.Spectrum)

                # Save the filtered DataFrame back as a CSV file
                CSV.write(joinpath(path2files, f), df)
            else
                println("Warning: 'Spectrum' column not found in file $(f). Skipping filtration.")
            end
        end
    end
    println("Done removing m/z:intensity pairs with intensity values <= 3% of the highest intensity values from spectra of ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4_low)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to filter out rows where "Classifications" from ChromaTOF contains unwanted values
function filter_classifications(path2files)
    nn = readdir(path2files)
    time_b4 = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Filtering away features classified by ChromaTOF as 'Bleed', 'Toluene', or 'DCM'... ")

    # Find and process each CSV file in the directory
    for f in nn
        if contains.(f, ".csv")
            df = CSV.read(joinpath(path2files,f), DataFrame)

            # Check if 'Classifications' column exists
            if "Classifications" in names(df)
                # Replace missing values in 'Classifications' with "NA"
                df.Classifications .= coalesce.(df.Classifications, "NA")

                # Delete features (rows) where 'Classifications' contains "Bleed", "Toluene", or "DCM"
                df_filtered = filter(row -> 
                    !(occursin(r"(?i)Bleed", row.Classifications) || 
                    occursin(r"(?i)Toluene", row.Classifications) || 
                    occursin(r"(?i)DCM", row.Classifications)), df)

                # Save the filtered DataFrame back as a CSV file
                CSV.write(joinpath(path2files, f), df_filtered)
            else
                println("Warning: 'Classifications' column not found in file $(f). Skipping filtration.")
            end
        end
    end

    println("Done filtering ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to filter out rows where "Name" contains unwanted values, e.g., siloxanes & other bleed
function filter_names(path2files, keywords::Vector{String})
    # Pre-escape keywords for regex safety
    # Ensures that keywords with special characters ([, ], -, ,) are correctly interpreted by the regex engine
    escaped_keywords = [Regex("(?i)" * escape_regex(k)) for k in keywords]
    nn = readdir(path2files)
    time_b4_filter = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Filtering away bleed features by name... ")
    println(" ")
    println("-"^150)

    # Find and process each CSV file in the directory
    for f in nn
        if occursin(".csv", f)
            println("Processing file: $f.")
            df = CSV.read(joinpath(path2files, f), DataFrame)

            # Check if 'Name' column exists
            if "Name" in names(df)
                # Replace missing values in 'Name' with "NA"
                df.Name .= coalesce.(df.Name, "NA")
                df.Name .= strip.(df.Name)  # Remove leading/trailing spaces

                # Print the original number of rows for debugging
                original_rows = nrow(df)
                println("Original number of rows: $original_rows.")
                # Find rows matching keywords for debugging
                matched_rows = filter(row -> any(regex -> occursin(regex, row.Name), escaped_keywords), df)
                println("Rows matching keywords: ", nrow(matched_rows))
                println("Matched rows preview (first 5 entries):")
                println(matched_rows[!, "Name"][1:min(end, 5)])  # Print the first 5 matched names for debugging
                println()

                # Delete features (rows) where 'Name' contains bleed compound name parts
                df_filtered = filter(row -> !any(regex -> occursin(regex, row.Name), escaped_keywords), df)
                
                # Print the number of rows after filtering
                filtered_rows = nrow(df_filtered)
                println("Number of rows after filtering: $filtered_rows.")

                # Save the filtered DataFrame back as a CSV file
                CSV.write(joinpath(path2files, f), df_filtered)

                # Print the difference for debugging
                removed_rows = original_rows - filtered_rows
                println("Rows removed: $removed_rows.")
                println("-"^150)
            else
                println("Warning: 'Name' column not found in file $(f). Skipping filtration.")
            end
        end
    end

    println("Done filtering ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4_filter)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to escape special regex characters
function escape_regex(str::String)
    # regex special characters that need to be escaped
    special_chars = r"([\\.+*?^=!:${}()|\[\]\/&])"  # the correct pattern
    return replace(str, special_chars => s"\\\1")
end



# A function to filter out rows where "Name" contains, e.g., compounds with spectra similar to toluene (m/z 91 or 92),
# if their spectra do not contains other higher m/z values, like molecular ions
function filter_names_if_missing_mz(path2files, filters::Vector{Tuple{String, Int}})
    # Pre-escape keywords for regex safety
    # Ensures that keywords with special characters ([, ], -, ,) are correctly interpreted by the regex engine
    escaped_filters = [(Regex("(?i)" * escape_regex(f[1])), f[2]) for f in filters]
    nn = readdir(path2files)
    time_b4_filtr = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Filtering away bleed features by name and specific m/z conditions... ")
    println(" ")
    println("-"^150)
    
    # Find and process each CSV file in the directory
    for f in nn
        if occursin(".csv", f)
            println("Processing file: $f.")
            df = CSV.read(joinpath(path2files, f), DataFrame)

            # Check if 'Name' and 'Spectrum' columns exist
            if "Name" in names(df) && "Spectrum" in names(df)
                # Replace missing values in 'Name' with "NA"
                df.Name .= coalesce.(df.Name, "NA")
                df.Name .= strip.(df.Name)  # Remove leading/trailing spaces

                # Print the original number of rows for debugging
                original_rows = nrow(df)
                println("Original number of rows: $original_rows.")

                # Find rows matching keywords for debugging
                matched_rows = filter(row -> any(regex -> occursin(regex[1], row.Name), escaped_filters), df)
                println("Rows matching keywords: ", nrow(matched_rows))
                # println("Matched rows preview:")
                # println(matched_rows[!, "Name"][1:min(end, 5)])  # Print the first 5 matched names for debugging
                # println(matched_rows[!, "Name"][1:end])  # Print all matched names for debugging
                println()

                # Print string (names) and integer (m/z's) keywords found in the input files
                for row in eachrow(df)
                    for regex in escaped_filters
                        if occursin(regex[1], row.Name)
                            mz_values = floor.(parse.(Float64, split(row.Spectrum, r"[: ]"))[1:2:end])
                            println("Found keyword: $(regex[1]) in 'Name': $(row.Name)")
                            if any(mz -> mz == regex[2], mz_values)
                                println("Target m/z $(regex[2]) found in 'Spectrum'")
                            end
                        end
                    end
                end
                println()

                # Delete features (rows) where 'Name' contains bleed compound name parts & DO NOT contain specific m/z's
                df_filtered = filter(row -> 
                    !any(regex -> occursin(regex[1], row.Name) && 
                    !any(floor.(parse.(Float64, split(row.Spectrum, r"[: ]"))[1:2:end]) .== regex[2]), escaped_filters), df)

                # Print the number of rows after filtering
                filtered_rows = nrow(df_filtered)
                println("Number of rows after filtering: $filtered_rows.")

                # Save the filtered DataFrame back as a CSV file
                CSV.write(joinpath(path2files, f), df_filtered)

                # Print the difference for debugging
                removed_rows = original_rows - filtered_rows
                println("Rows removed: $removed_rows.")
                println("-"^150)
            else
                println("Warning: 'Name' or 'Spectrum' column not found in file $(f). Skipping filtration.")
            end
        end
    end

    println("Done filtering ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4_filtr)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to filter out rows based on 'Name' containing "Peak" and specific m/z values for the highest intensity pair in 'Spectrum'
function filter_peak_names_with_bleed_mz(path2files, bleed_filters::Vector{Tuple{String, Int}})
    # Pre-escape keywords for regex safety
    escaped_filters = [(Regex("(?i)" * escape_regex(f[1])), f[2]) for f in bleed_filters]
    nn = readdir(path2files)
    time_b4_flt = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Filtering features by 'Name' containing 'Peak' and specific m/z conditions... ")
    println(" ")
    println("-"^150)

    # Find and process each CSV file in the directory
    for f in nn
        if occursin(".csv", f)
            println("Processing file: $f.")
            df = CSV.read(joinpath(path2files, f), DataFrame)

            # Check if 'Name' and 'Spectrum' columns exist
            if "Name" in names(df) && "Spectrum" in names(df)
                # Replace missing values in 'Name' with "NA"
                df.Name .= coalesce.(df.Name, "NA")
                df.Name .= strip.(df.Name)  # Remove leading/trailing spaces

                # Print the original number of rows for debugging
                original_rows = nrow(df)
                println("Original number of rows: $original_rows.")

                # Find rows matching keywords for debugging
                matched_rows = filter(row -> any(regex -> occursin(regex[1], row.Name), escaped_filters), df)
                println("Rows matching keywords: ", nrow(matched_rows))
                println()

                # Print string (names) and integer (m/z's) keywords found in the input files
                for row in eachrow(df)
                    for regex in escaped_filters
                        if occursin(regex[1], row.Name)
                            pairs = split(row.Spectrum, ' ')
                            intensities = [parse(Float64, split(pair, ':')[2]) for pair in pairs]
                            max_intensity_index = argmax(intensities)
                            max_mz = floor(parse(Float64, split(pairs[max_intensity_index], ':')[1]))
                            if max_mz == regex[2]
                                println("Removing row with 'Name': $(row.Name) and base m/z: $(max_mz)")
                            end
                        end
                    end
                end
                println()

                # Delete features (rows) where 'Name' contains "Peak" and the highest intensity m/z is one of the specified values
                df_filtered = filter(row -> 
                    !any(regex -> occursin(regex[1], row.Name) && 
                    floor(parse(Float64, split(split(row.Spectrum, ' ')[argmax([parse(Float64, split(pair, ':')[2]) for pair in split(row.Spectrum, ' ')])], ':')[1])) == regex[2], escaped_filters), df)

                # Print the number of rows after filtering
                filtered_rows = nrow(df_filtered)
                println("Number of rows after filtering: $filtered_rows.")

                # Save the filtered DataFrame back as a CSV file
                CSV.write(joinpath(path2files, f), df_filtered)

                # Print the difference for debugging
                removed_rows = original_rows - filtered_rows
                println("Rows removed: $removed_rows.")
                println("-"^150)
            else
                println("Warning: 'Name' or 'Spectrum' column not found in file $(f). Skipping filtration.")
            end
        end
    end

    println("Done filtering ", length(filter(f -> endswith(f, ".csv") && !occursin(r"Original|merged|IntMerged|Aligned", f), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4_flt)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to filter out rows after alignment; similar to filter_names_if_missing_mz
function filter_names_mz_postalign(path2files, filters::Vector{Tuple{String, Int}})
    # Pre-escape keywords for regex safety
    # Ensures that keywords with special characters ([, ], -, ,) are correctly interpreted by the regex engine
    escaped_filters = [(Regex("(?i)" * escape_regex(f[1])), f[2]) for f in filters]
    nn = readdir(path2files)
    time_b4_filt = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Filtering away bleed features by name and specific m/z conditions... ")
    println(" ")
    println("-"^150)

    # Find and process each CSV file in the directory
    for f in nn
        if occursin("Aligned.csv", f)
            println("Processing file: $f.")
            df = CSV.read(joinpath(path2files, f), DataFrame)

            # Check if 'Name' and 'Spectrum' columns exist
            if "Name" in names(df) && "Spectrum_Consensus" in names(df)
                # Replace missing values in 'Name' with "NA"
                df.Name .= coalesce.(df.Name, "NA")
                df.Name .= strip.(df.Name)  # Remove leading/trailing spaces

                # Print the original number of rows for debugging
                original_rows = nrow(df)
                println("Original number of rows: $original_rows.")

                # Find rows matching keywords for debugging
                matched_rows = filter(row -> any(regex -> occursin(regex[1], row.Name), escaped_filters), df)
                println("Rows matching keywords: ", nrow(matched_rows))
                # println("Matched rows preview:")
                # println(matched_rows[!, "Name"][1:min(end, 5)])  # Print the first 5 matched names for debugging
                # println(matched_rows[!, "Name"][1:end])  # Print all matched names for debugging
                println()

                # Print string (names) and integer (m/z's) keywords found in the input files
                for row in eachrow(df)
                    for regex in escaped_filters
                        if occursin(regex[1], row.Name)
                            mz_values = floor.(parse.(Float64, split(row.Spectrum_Consensus, r"[: ]"))[1:2:end])
                            println("Found keyword: $(regex[1]) in 'Name': $(row.Name)")
                            if any(mz -> mz == regex[2], mz_values)
                                println("Target m/z $(regex[2]) found in 'Spectrum_Consensus'")
                            end
                        end
                    end
                end
                println()

                # Delete features (rows) where 'Name' contains bleed compound name parts & DO NOT contain specific m/z's
                df_filtered = filter(row -> 
                    !any(regex -> occursin(regex[1], row.Name) && 
                    !any(floor.(parse.(Float64, split(row.Spectrum_Consensus, r"[: ]"))[1:2:end]) .== regex[2]), escaped_filters), df)

                # Print the number of rows after filtering
                filtered_rows = nrow(df_filtered)
                println("Number of rows after filtering: $filtered_rows.")

                # Save the filtered DataFrame back as a CSV file
                filtered_filename = f[1:end-4] * "_Fltrd.csv"
                CSV.write(joinpath(path2files, filtered_filename), df_filtered)
                println("Filtered file saved as: $filtered_filename")

                # Print the difference for debugging
                removed_rows = original_rows - filtered_rows
                println("Rows removed: $removed_rows.")
                println("-"^150)
            else
                println("Warning: 'Name' or 'Spectrum_Consensus' column not found in file $(f). Skipping filtration.")
            end
        end
    end
    
    println("Done filtering ", length(filter(f -> endswith(f, "_Fltrd.csv"), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4_filt)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to merge features with the same name, base mass and dRI < threshold after alignment
function merge_postalign(path2files::String, threshold::Float64, weights::Vector{Float64}=[1.0, 1.0, 1.0])
    nn = readdir(path2files)
    time_b4_mrg = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Merging features with the same name, base mass and dRI < N... ")
    println(" ")
    println("-"^150)

    for f in nn
        if occursin("_Fltrd.csv", f)
            println("-"^150)
            println("Processing file: $f.")
            println()
            df = CSV.read(joinpath(path2files, f), DataFrame; stringtype=String)

            # Fill empty cells with 'NA'
            df = coalesce.(df, "NA")

            # Delete specified columns
            select!(df, Not([:QuantMasses, :ActualMasses, :All_QuantMasses, :All_ActualMasses]))

            # Sort the DataFrame by 'Name'
            sort!(df, :Name)

            initial_rows = nrow(df)
            println("Initial number of rows: $initial_rows")
            # println()
            # println("Input DataFrame: ")
            # println(df)
            # println()
            println("-"^150)

            merged_rows = DataFrame()
            initial_features = DataFrame()

            i = 1
            while i <= nrow(df)
                row = df[i, :]
                j = i + 1
                while j <= nrow(df)
                    match = df[j, :]
                    if strip(row.Name) == strip(match.Name) && floor(Int, row.BaseMass) == floor(Int, match.BaseMass) && abs(row.MedRI - match.MedRI) < threshold
                        println("Match found for row $i with row $j: $(row.Name)")

                        initial_features = vcat(initial_features, DataFrame(row), DataFrame(match))

                        # Merge rows
                        row.MinRt1D = min(row.MinRt1D, match.MinRt1D)
                        row.MaxRt1D = max(row.MaxRt1D, match.MaxRt1D)
                        row.AveRt1D = (row.MinRt1D + row.MaxRt1D) / 2

                        row.MinRt2D = min(row.MinRt2D, match.MinRt2D)
                        row.MaxRt2D = max(row.MaxRt2D, match.MaxRt2D)
                        row.AveRt2D = (row.MinRt2D + row.MaxRt2D) / 2

                        row.MinRI = min(row.MinRI, match.MinRI)
                        row.MaxRI = max(row.MaxRI, match.MaxRI)
                        row.MedRI = median([row.MinRI, row.MaxRI])

                        row.MinRIlib = coalesce(min(row.MinRIlib, match.MinRIlib), row.MinRIlib, match.MinRIlib)
                        row.MaxRIlib = coalesce(max(row.MaxRIlib, match.MaxRIlib), row.MaxRIlib, match.MaxRIlib)
                        row.MedRIlib = median([row.MinRIlib, row.MaxRIlib])

                        row.Area += match.Area
                        row.TotalArea += match.TotalArea

                        row.BaseMass = median([row.BaseMass, match.BaseMass])
                        row.Similarity = max(row.Similarity, match.Similarity)
                        row.Probability = max(row.Probability, match.Probability)
                        row.dRI = max(row.dRI, match.dRI)

                        row.Score = (row.Similarity * weights[1] + row.Probability * weights[2] + row.dRI * weights[3]) / sum(weights)

                        row.Spectrum_BestMatch = ifelse(row.Probability > match.Probability, row.Spectrum_BestMatch, match.Spectrum_BestMatch)
                        row.Spectrum_Consensus = ifelse(row.Probability > match.Probability, row.Spectrum_Consensus, match.Spectrum_Consensus)

                        row.All_Names *= " &&& " * match.All_Names
                        row.All_Formulae *= " &&& " * match.All_Formulae
                        row.All_InChIKey *= " &&& " * match.All_InChIKey
                        row.All_TotalArea *= " &&& " * match.All_TotalArea
                        row.All_BaseMass *= " &&& " * match.All_BaseMass
                        row.All_Origin *= " &&& " * match.All_Origin
                        row.All_Class *= " &&& " * match.All_Class
                        row.All_CAS *= " &&& " * match.All_CAS
                        row.All_Similarity *= " &&& " * match.All_Similarity

                        area_columns = names(df)[findall(x -> x == "Spectrum_Consensus", names(df))[1]+1:end]
                        for col in area_columns
                            row[col] = coalesce(row[col], 0) + coalesce(match[col], 0)
                        end

                        # row.DetectionNum = count(!isnan, row[area_columns])
                        row.DetectionNum = count(x -> !isnan(x) && x > 0, row[area_columns])
                        row.DetectionFreq = round(row.DetectionNum / length(area_columns) * 100)

                        # Remove the matched row
                        df = df[Not(j), :]
                    else
                        j += 1
                    end
                end
                merged_rows = vcat(merged_rows, DataFrame(row))
                i += 1
            end

            merged_rows.Nr = 1:nrow(merged_rows)

            println()
            println("-"^150)
            # println("Merged DataFrame: ")
            # println(merged_rows)
            # println("*"^150)            

            merged_rows_count = nrow(initial_features) // 2
            final_rows = nrow(merged_rows)
            println("Initial number of rows: $initial_rows")
            # println("Number of rows merged: $merged_rows_count")
            println("Final number of rows: $final_rows")
            row_difference = initial_rows - final_rows
            println("Number of rows difference: $row_difference")
            println("-"^150)

            if nrow(initial_features) > 0
                CSV.write(joinpath(path2files, replace(f, "_Fltrd" => "_Fltrd_MrgdFtrs")), initial_features)                
            else
                println("No features were merged.")
                println("!"^150)
            end

            CSV.write(joinpath(path2files, replace(f, "_Fltrd" => "_Fltrd_Mrgd")), merged_rows)            
            println("Input file: $f.")
            println("File with merged features saved as: $(replace(f, "_Fltrd" => "_Fltrd_MrgdFtrs")).")
            println("Processed file saved as: $(replace(f, "_Fltrd" => "_Fltrd_Mrgd")).")
            println("-"^150)
        end
    end

    println("Done merging ", length(filter(f -> endswith(f, "_Mrgd.csv"), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4_mrg)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to retrieve InChIKey, SMILES & CAS# values from PubChem based on "InChIKey" and "Name" column values
# A helper function to get CID from InChIKey or Name
function get_pubchem_cid(identifier::Union{String, Nothing}, type::Symbol)
    if identifier === nothing || identifier in ["NA", "NoID", ""]
        return nothing
    end

    identifier = strip(identifier)  # remove leading/trailing spaces
    baseurl = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"

    url = type == :inchikey ?
        "$baseurl/inchikey/$identifier/cids/JSON" :
        "$baseurl/name/$identifier/cids/JSON"

    try
        resp = HTTP.get(url)
        if resp.status != 200
            return nothing
        end

        data = JSON3.read(String(resp.body))
        if haskey(data, :IdentifierList) && haskey(data[:IdentifierList], :CID)
            cids = data[:IdentifierList][:CID]
            return cids[1]  # take first CID
        else
            return nothing
        end

    catch e
        return nothing
    end
end



# Safe wrapper: try InChIKey first, fallback to Name
function safe_get_cid(inchikey::Union{String, Nothing}, name::Union{String, Nothing})
    cid = get_pubchem_cid(inchikey, :inchikey)

    if cid === nothing
        cid = get_pubchem_cid(name, :name)
    end

    return cid
end



# A CAS# retrieval function
function get_main_cas(cid::Int)
    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/$cid/JSON"
    try
        resp = HTTP.get(url; readtimeout=15)
        if resp.status != 200
            return "NA"
        end
        data = JSON3.read(String(resp.body))

        # Recursive function to search sections for "CAS" heading
        function extract_cas(sections)
            for sec in sections
                # check if SectionTitle contains "CAS"
                if haskey(sec, :TOCHeading) && occursin(r"CAS", sec[:TOCHeading])
                    # look in Information list
                    if haskey(sec, :Information)
                        for info in sec[:Information]
                            if haskey(info, :Value) && haskey(info[:Value], :StringWithMarkup)
                                for s in info[:Value][:StringWithMarkup]
                                    text = s[:String]
                                    m = match(r"\d{2,7}-\d{2}-\d", text)
                                    if m !== nothing
                                        return m.match  # first match
                                    end
                                end
                            end
                        end
                    end
                end
                # recurse into subsections
                if haskey(sec, :Section)
                    cas = extract_cas(sec[:Section])
                    if cas !== nothing
                        return cas
                    end
                end
            end
            return nothing
        end

        sections = get(data, :Record, Dict())[:Section]  # top-level sections
        cas = extract_cas(sections)
        return cas === nothing ? "NA" : cas

    catch
        return "NA"
    end
end



# The main function to retrieve data from PubChem
function PubChemRetriever(path2files::String)
    files = readdir(path2files)
    time_b4_rtrv = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Retrieving InChIKey, SMILES, CAS#, DTXSID, Uses, and Use Classification values from PubChem... ")
    println("-"^150)
    println(" ")

    for f in files
        if occursin("_Mrgd.csv", f)
            println("-"^150)
            println("Processing file: $f.")
            println()

            df = CSV.read(joinpath(path2files, f), DataFrame; stringtype=String)

            # Add empty columns for PubChem data
            df.InChIKey_PubChem = fill("NA", nrow(df))
            df.SMILES_PubChem = fill("NA", nrow(df))
            df.CAS_PubChem = fill("NA", nrow(df))
            df.InChIKey_Consensus = fill("NA", nrow(df))
            df.CAS_Consensus = fill("NA", nrow(df))

            println("Retrieving data from PubChem. This will take a while... \n")
            println("-"^150)

            for (i, row) in enumerate(eachrow(df))
                name = row[:Name]
                inchikey = row[:InChIKey]
                
                # Skip if "Peak" in the compound name
                if occursin("Peak", name) || occursin("Feature", name)
                    println("$i: $name --> Skipped: non-hit ('Name' contains 'Peak' or 'Feature')")
                    println("-"^150)
                    continue
                end

                try
                    cid = safe_get_cid(inchikey, name)

                    # If CID is still nothing, report and skip
                    if cid === nothing
                        println("$i: $name --> Not found in PubChem by Name or InChIKey")
                        continue
                    end

                    if cid != nothing
                        # Retrieve InChIKey and SMILES
                        properties_json = JSON3.read(get_for_cids(cid, properties="SMILES,InChIKey", output="JSON"))                        
                        if haskey(properties_json, :PropertyTable) && haskey(properties_json[:PropertyTable], :Properties)
                            props = properties_json[:PropertyTable][:Properties][1]
                            # Check both old CanonicalSMILES and new SMILES
                            smiles_val = haskey(props, :SMILES) ? props[:SMILES] :
                                         haskey(props, :CanonicalSMILES) ? props[:CanonicalSMILES] : "NA"
                            df.SMILES_PubChem[i] = smiles_val
                            df.InChIKey_PubChem[i] = get(props, :InChIKey, "NA")
                        else
                            df.SMILES_PubChem[i] = "NA"
                            df.InChIKey_PubChem[i] = "NA"
                        end

                        # Retrieve CAS numbers
                        df.CAS_PubChem[i] = get_main_cas(cid)
                                                
                        println("$i: $name --> InChIKey: $(df.InChIKey_PubChem[i]) // SMILES: $(df.SMILES_PubChem[i]) // CAS#: $(df.CAS_PubChem[i])")
                        println("-"^150)
                    else
                        println("$i: $name --> Not found in PubChem")
                    end
                catch e
                    println("Error processing $name: $e")
                    println("!"^150)
                end
            end

            # Copy all InChIKey_PubChem and CAS_PubChem values to InChIKey_Consensus and CAS_Consensus
            df.InChIKey_Consensus .= df.InChIKey_PubChem
            df.CAS_Consensus .= df.CAS_PubChem

            # Update InChIKey_Consensus based on conditions: if "NA" in 'InChIKey_PubChem', copy from 'InChIKey'
            for i in 1:nrow(df)
                if df.InChIKey_Consensus[i] == "NA" && !(occursin("Peak", df.Name[i]) || occursin("Feature", df.Name[i])) && !(df.InChIKey[i] in ["NA", "NoID"])
                    df.InChIKey_Consensus[i] = df.InChIKey[i]
                end
            end

            # Update 'CAS_Consensus' based on conditions: if "NA" in 'CAS_PubChem', copy from 'CAS'
            for i in 1:nrow(df)
                if df.CAS_Consensus[i] == "NA" && !(occursin("Peak", df.Name[i]) || occursin("Feature", df.Name[i])) && !(df.CAS[i] in ["NA", "NoID"])
                    df.CAS_Consensus[i] = df.CAS[i]
                end
            end

            # Sort by InChIKey_Consensus in descending order, NAs at the bottom
            sort!(df, :InChIKey_Consensus, rev=true, by=x->x == "NA" ? "" : x)

            output_file = replace(f, ".csv" => "_PubChem.csv")
            CSV.write(joinpath(path2files, output_file), df)

            println("-"^150)
            println("Processing complete. Saved to $output_file.")
            println("-"^150)
        end
    end

    println("Done retrieving data for ", length(filter(f -> endswith(f, "_PubChem.csv"), readdir(path2files))), " .csv files. ")
    println("Start time: ", time_b4_rtrv)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# # A function to retrieve functional use data from EPA FUse, NORMAN dust trial & in-house harmonized use databases
# # The databases should be in TXT (tab-delimited) format
# function ClassFinder(path_file, path_harm, path_fuse, path_norm)
#     time_b4_cls = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
#     println("Retrieving functional use classification from EPA FUse, NORMAN dust trial & in-house harmonized use databases... ")
#     println("-"^150)
#     println(" ")
#     println("Reading input file...")
#     println("-"^150)
#     # Reading the input file
#     input_df = CSV.read(path_file, header=1, DataFrame)
#     println("Input file read. Rows: ", nrow(input_df), " Columns: ", ncol(input_df))
#     println("-"^150)

#     # Strip any extra spaces from column names
#     for col in names(input_df)
#         rename!(input_df, Symbol(col) => Symbol(strip(col)))
#     end
#     println("Input DF columns after stripping spaces: \n", names(input_df))
#     println("-"^150)

#     # Adding empty classification columns
#     input_df[:, :Class_Harmonized] = fill("", nrow(input_df))
#     input_df[:, :Class_FUse] = fill("", nrow(input_df))
#     input_df[:, :Class_NORMAN] = fill("", nrow(input_df))

#     # println("Input file: \n", input_df)
#     # println("-"^150)

#     # Ensure Name column exists in input_df and rename it to ensure column is correctly identified
#     if :Name in names(input_df)
#         println("Found column :Name")
#     else
#         println("Column :Name not found, renaming...")
#         rename!(input_df, Symbol("Name") => :Name)
#     end

#     # Ensure InChIKey_PubChem column exists in input_df and rename it to ensure column is correctly identified
#     if :InChIKey_Consensus in names(input_df)
#         println("Found column :InChIKey_Consensus")
#     else
#         println("Column :InChIKey_Consensus not found, renaming...")
#         rename!(input_df, Symbol("InChIKey_Consensus") => :InChIKey_Consensus)
#     end

#     # Ensure CAS_PubChem column exists in input_df and rename it to ensure column is correctly identified
#     if :CAS_PubChem in names(input_df)
#         println("Found column :CAS_Consensus")
#     else
#         println("Column :CAS_Consensus not found, renaming...")
#         rename!(input_df, Symbol("CAS_Consensus") => :CAS_Consensus)
#     end

#     println("-"^150)
#     println("Input DF columns: \n", names(input_df))
#     println("-"^150)

#     # Reading classification files
#     function read_classification_file(path)
#         println("Reading classification file: ", path)
#         println("-"^150)
        
#         # Check file extension and read accordingly
#         if any(ext -> endswith(lowercase(path), ext), [".csv", ".txt"])
#             df = CSV.read(path, DataFrame; delim='\t')
#         elseif any(ext -> endswith(lowercase(path), ext), [".xls", ".xlsx"])
#             df = DataFrame(XLSX.readtable(path, 1)[])
#         else
#             error("Unsupported file format: $path")
#             println("!"^150)
#         end
        
#         # Strip any extra spaces from column names
#         for col in names(df)
#             rename!(df, Symbol(col) => Symbol(strip(col)))
#         end
        
#         println("Classification DF columns after stripping spaces: ", names(df))
#         println("-"^150)
        
#         # Remove quotes from the 'Name' column if present
#         if :Name in names(df)
#             println("Stripping quotes from the 'Name' column...")
#             df[!, :Name] .= replace.(df[!, :Name], r"^\"(.*)\"$" => s -> s.match[1])
#         end
        
#         # Rename columns to match the input file if needed
#         if :InChIKey in names(df)
#             rename!(df, :InChIKey => :InChIKey_Consensus)
#         end
        
#         if :CAS in names(df)
#             rename!(df, :CAS => :CAS_Consensus)
#         end
        
#         return df
#     end

#     class_harm = read_classification_file(path_harm)
#     class_fuse = read_classification_file(path_fuse)
#     class_norm = read_classification_file(path_norm)

#     # Rename columns in classification files for consistency
#     rename!(class_harm, :InChIKey => :InChIKey_Consensus, :CAS => :CAS_Consensus)
#     rename!(class_fuse, :InChIKey => :InChIKey_Consensus, :CAS => :CAS_Consensus)
#     rename!(class_norm, :InChIKey => :InChIKey_Consensus, :CAS => :CAS_Consensus)

#     println("Classification files read. Harmonized rows: ", nrow(class_harm), ", FUse rows: ", nrow(class_fuse), ", NORMAN rows: ", nrow(class_norm))
#     println("-"^150)
#     println("Harmonized cols: ", names(class_harm), 
#         "\nFUse cols: ", names(class_fuse), 
#         "\nNORMAN cols: ", names(class_norm))
#     println("-"^150)

#     # Classification function
#     function update_class!(input_df, class_df, key_col, class_col, target_col)
#         class_map = Dict(class_df[!, key_col] .=> class_df[!, class_col])
#         println("Updating classes for key column: ", key_col, " and target column: ", target_col)
#         println("-"^150)
    
#         for row in eachrow(input_df)
#             key = row[key_col]
#             # Check if the key is missing or empty
#             if key == "" || ismissing(key)
#                 println("Skipping row with missing or empty key: ", row)
#                 println("-"^150)
#                 continue
#             end

#             # Skip rows where :Name contains "Peak" or "Feature"
#             if key_col == "Name" && occursin.(r"(Peak|Feature)", [row[:Name]])[1]
#                 println("Skipping entry with name containing 'Peak' or 'Feature': ", row[:Name])
#                 continue
#             end
    
#             if row[target_col] == "" || ismissing(row[target_col])
#                 if haskey(class_map, key)
#                     # If key exists in the class_map, update the target column
#                     row[target_col] = class_map[key]
#                 else
#                     println("No match found for key: ", key, " in ", target_col)
#                     println("-"^150)
#                     # Update with "NA" if no match found
#                     row[target_col] = "NA"
#                 end
#             end
#         end
#     end

#     # Update classes based on Name, InChIKey_Consensus, and CAS_Consensus
#     println("Updating classes based on Name, InChIKey_Consensus, and CAS_Consensus...")
#     println("-"^150)
#     for (key_col, target_col, class_col, class_df) in [
#         ("Name", "Class_Harmonized", "Class", class_harm),
#         ("Name", "Class_FUse", "harmonized_function", class_fuse),
#         ("Name", "Class_NORMAN", "Class", class_norm),
#         ("InChIKey_Consensus", "Class_Harmonized", "Class", class_harm),
#         ("InChIKey_Consensus", "Class_FUse", "harmonized_function", class_fuse),
#         ("InChIKey_Consensus", "Class_NORMAN", "Class", class_norm),
#         ("CAS_Consensus", "Class_Harmonized", "Class", class_harm),
#         ("CAS_Consensus", "Class_FUse", "harmonized_function", class_fuse),
#         ("CAS_Consensus", "Class_NORMAN", "Class", class_norm)
#     ]
#         update_class!(input_df, class_df, key_col, class_col, target_col)
#     end

#     # Fill empty cells with 'NA'
#     println("Filling empty cells with 'NA'...")
#     println("-"^150)
#     for col in names(input_df)
#         if eltype(input_df[!, col]) <: AbstractString
#             replace!(input_df[!, col], "" => "NA")
#         end
#     end

#     # Save the output file
#     output_path = replace(path_file, ".csv" => "_Class.csv")
#     CSV.write(output_path, input_df)
#     println("Classification completed successfully for ", length(filter(f -> endswith(f, "_CFB.csv"), readdir(path2files))), " .csv files. ")
#     println("File path: $output_path")
#     println("-"^150)
#     println("Start time: ", time_b4_cls)
#     println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
#     println("*"^150)
# end





# A helper function to find first matching column name and rename to a standard name
function standardize_column!(df::DataFrame, target::Symbol, variants::Vector{String})
    current_names = String.(names(df))

    # Already exists
    if String(target) in current_names
        println("Found standard column: $target")
        return
    end

    # Search variants (case-insensitive)
    lookup = Dict(lowercase(nm) => nm for nm in current_names)

    for var in variants
        if haskey(lookup, lowercase(var))
            old_name = Symbol(lookup[lowercase(var)])
            rename!(df, old_name => target)
            println("Renamed $(old_name) -> $(target)")
            return
        end
    end
    println("WARNING: Could not find column for $(target)")
end



# A function to retrieve functional use data from EPA FUse, NORMAN dust trial & in-house harmonized use databases
# The databases should be in TXT (tab-delimited) format
function ClassFinder(path_file, path_harm, path_fuse, path_norm)
    time_b4_cls = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Retrieving functional use classification from EPA FUse, NORMAN dust trial & in-house harmonized use databases... ")
    println("-"^150)
    println(" ")
    println("Reading input file...")
    println("-"^150)
    # Reading the input file
    input_df = CSV.read(path_file, header=1, DataFrame)
    println("Input file read. Rows: ", nrow(input_df), " Columns: ", ncol(input_df))
    println("-"^150)

    # Strip any extra spaces from column names
    for col in names(input_df)
        rename!(input_df, Symbol(col) => Symbol(strip(col)))
    end
    println("Input DF columns after stripping spaces: \n", names(input_df))
    println("-"^150)

    # Adding empty classification columns
    input_df[:, :Class_Harmonized] = fill("", nrow(input_df))
    input_df[:, :Class_FUse] = fill("", nrow(input_df))
    input_df[:, :Class_NORMAN] = fill("", nrow(input_df))

    # println("Input file: \n", input_df)
    # println("-"^150)

    # Ensure Name column exists in input_df and rename it to ensure column is correctly identified
    if :Name in names(input_df)
        println("Found column :Name")
    else
        println("Column :Name not found, renaming...")
        rename!(input_df, Symbol("Name") => :Name)
    end
    
    # Standardize identifier columns
    standardize_column!(input_df, :Name,
        ["Name", "Metabolite name", "Metabolite Name", "chemical_name",
            "Feature name", "feature_name", "Feature Name", "Compound name", "Compound Name"])

    standardize_column!(input_df, :CAS_Consensus,
        ["CAS_Consensus", "CAS", "CASRN", "CAS#"])

    standardize_column!(input_df, :InChIKey_Consensus,
        ["InChIKey_Consensus", "InChIKey", "inchikey"])

    println("-"^150)
    println("Input DF columns: \n", names(input_df))
    println("-"^150)

    # Reading classification files
    function read_classification_file(path)
        println("Reading classification file: ", path)
        println("-"^150)
        
        # Check file extension and read accordingly
        if any(ext -> endswith(lowercase(path), ext), [".csv", ".txt"])
            df = CSV.read(path, DataFrame; delim='\t')
        elseif any(ext -> endswith(lowercase(path), ext), [".xls", ".xlsx"])
            df = DataFrame(XLSX.readtable(path, 1)[])
        else
            error("Unsupported file format: $path")
            println("!"^150)
        end
        
        # Strip any extra spaces from column names
        for col in names(df)
            rename!(df, Symbol(col) => Symbol(strip(col)))
        end
        
        println("Classification DF columns after stripping spaces: ", names(df))
        println("-"^150)
        
        # Remove quotes from the 'Name' column if present
        if :Name in names(df)
            println("Stripping quotes from the 'Name' column...")
            df[!, :Name] .= replace.(df[!, :Name], r"^\"(.*)\"$" => s -> s.match[1])
        end
        
        # Rename columns to match the input file if needed
        standardize_column!(df, :InChIKey_Consensus,
            ["InChIKey_Consensus", "InChIKey", "inchikey"])

        standardize_column!(df, :CAS_Consensus,
            ["CAS_Consensus", "CAS", "CASRN", "CAS#"])
    end
    
    class_harm = read_classification_file(path_harm)
    class_fuse = read_classification_file(path_fuse)
    class_norm = read_classification_file(path_norm)

    println("Classification files read. Harmonized rows: ", nrow(class_harm), ", FUse rows: ", nrow(class_fuse), ", NORMAN rows: ", nrow(class_norm))
    println("-"^150)
    println("Harmonized cols: ", names(class_harm), 
        "\nFUse cols: ", names(class_fuse), 
        "\nNORMAN cols: ", names(class_norm))
    println("-"^150)

    # Classification function
    function update_class!(input_df, class_df, key_col, class_col, target_col)
        class_map = Dict(class_df[!, key_col] .=> class_df[!, class_col])
        println("Updating classes for key column: ", key_col, " and target column: ", target_col)
        println("-"^150)
    
        for row in eachrow(input_df)
            key = row[key_col]
            # Check if the key is missing or empty
            if key == "" || ismissing(key)
                println("Skipping row with missing or empty key: ", row)
                println("-"^150)
                continue
            end

            # Skip rows where :Name contains "Peak" or "Feature"
            if key_col == "Name" && occursin.(r"(Peak|Feature)", [row[:Name]])[1]
                println("Skipping entry with name containing 'Peak' or 'Feature': ", row[:Name])
                continue
            end
    
            if row[target_col] == "" || ismissing(row[target_col])
                if haskey(class_map, key)
                    # If key exists in the class_map, update the target column
                    row[target_col] = class_map[key]
                else
                    println("No match found for key: ", key, " in ", target_col)
                    println("-"^150)
                    # Update with "NA" if no match found
                    row[target_col] = "NA"
                end
            end
        end
    end

    # Update classes based on Name, InChIKey_Consensus, and CAS_Consensus
    println("Updating classes based on Name, InChIKey_Consensus, and CAS_Consensus...")
    println("-"^150)
    for (key_col, target_col, class_col, class_df) in [
        ("Name", "Class_Harmonized", "Class", class_harm),
        ("Name", "Class_FUse", "harmonized_function", class_fuse),
        ("Name", "Class_NORMAN", "Class", class_norm),
        ("InChIKey_Consensus", "Class_Harmonized", "Class", class_harm),
        ("InChIKey_Consensus", "Class_FUse", "harmonized_function", class_fuse),
        ("InChIKey_Consensus", "Class_NORMAN", "Class", class_norm),
        ("CAS_Consensus", "Class_Harmonized", "Class", class_harm),
        ("CAS_Consensus", "Class_FUse", "harmonized_function", class_fuse),
        ("CAS_Consensus", "Class_NORMAN", "Class", class_norm)
    ]
        update_class!(input_df, class_df, key_col, class_col, target_col)
    end

    # Fill empty cells with 'NA'
    println("Filling empty cells with 'NA'...")
    println("-"^150)
    for col in names(input_df)
        if eltype(input_df[!, col]) <: AbstractString
            replace!(input_df[!, col], "" => "NA")
        end
    end

    # Save the output file
    output_path = replace(path_file, ".csv" => "_Class.csv")
    CSV.write(output_path, input_df)
    println("Classification completed successfully for ", length(filter(f -> endswith(f, "_CFB.csv"), readdir(path2files))), " .csv files. ")
    println("File path: $output_path")
    println("-"^150)
    println("Start time: ", time_b4_cls)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to search library/database postalignment to potentially identify some unknowns and verify previous identifications
function librarySearch(pathDB, pathFile, addRIwin, specID, mz_tol, numfrags = 15, index = "all", similarity_method::String = "DISCO")
    println("*"^150)
    println("Searching library/database... It will take a while... ")
    time_b4_search = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("-"^150)

    # Load database and overview files
    xf = XLSX.readxlsx(pathDB)
    sn = XLSX.sheetnames(xf)
    db = DataFrame(XLSX.readtable(pathDB, sn[1], infer_eltypes=true))
    # Create decided RI columns
    RIval, RItol = mergeToLeft(Matrix(db[:, 3:12])[:, 1:3:end], addRIwin)
    db[!,"RImatching"] = RIval
    db[!,"RImatchingTol"] = RItol
    spec = []

    # Iterate candidate spectra
    df = CSV.read(pathFile, DataFrame)
    df[!, :Class] .= coalesce.(df[!, :Class], "NA")
    df[!, :All_Class] .= coalesce.(df[!, :All_Class], "NA")
    df[!, :Origin] .= coalesce.(df[!, :Origin], "NA")
    df[!, :All_Origin] .= coalesce.(df[!, :All_Origin], "NA")
    if specID == "BM"
        println("Getting best match spectra")
        println("-"^150)
        spec = df[:, "Spectrum_BestMatch"]
        namekey = "Name"
        BaseMassKey = "BaseMass"
        OriginKey = "Origin"
        ClassKey = "Class"
    elseif specID == "CS"
        println("Getting consensus spectra")
        println("-"^150)
        spec = df[:, "Spectrum_Consensus"]
        namekey = "All_Names"
        BaseMassKey = "All_BaseMass"
        OriginKey = "All_Origin"
        ClassKey = "All_Class"
        for i = 1:size(df, 1)
            if contains(df[i, namekey], " &&& ")
                df[i, namekey] = split(df[i, namekey], " &&& ")[1]
            end
            
            if contains(df[i, BaseMassKey], " &&& ")
                df[i, BaseMassKey] = split(df[i, BaseMassKey], " &&& ")[1]
            end
            
            if contains(df[i, OriginKey], " &&& ")
                df[i, OriginKey] = split(df[i, OriginKey], " &&& ")[1]
            end
            
            if contains(df[i, ClassKey], " &&& ")
                df[i, ClassKey] = split(df[i, ClassKey], " &&& ")[1]
            end
        end
    else
        println("Provide a valid specID argument: \"BM\" or \"CS\"")
        println("!"^150)
        return
    end

    # Setup empty dataframe
    identAll = DataFrame(Align_Num = 0, DB_Num = 0, minRI = 0., maxRI = 0., DB_RI = 0., dot = 0., BaseMass = 0., Origin = "", Class = "",  Align_Name = "", DB_Name = "", DB_InChIKey = "", DB_CAS = "", DB_Formula = "", Spec = "", DB_spec = "", )

    if typeof(index) == String
        index = 1:size(df, 1)
    end

    for i in index
        println(i)
        sp = getSpec(spec[i])

        # Add/remove tolerance on a vector below
        # Select candidates
        if isnan(df[i, "MinRIlib"])
            minRI = fill(df[i, "MinRI"], size(db, 1)) - db[!, "RImatchingTol"]
        else
            minRI = fill(minimum([df[i, "MinRI"], df[i, "MinRIlib"]]), size(db, 1)) - db[!, "RImatchingTol"] # window is added to library info
        end
        if isnan(df[i, "MaxRIlib"])
            maxRI = fill(df[i, "MaxRI"], size(db, 1)) + db[!, "RImatchingTol"]
        else
            maxRI = fill(maximum([df[i, "MaxRI"], df[i, "MaxRIlib"]]), size(db, 1)) + db[!, "RImatchingTol"]
        end

        ind_s = findall(minRI .<= db[!, "RImatching"] .<= maxRI)
        if isempty(ind_s)
            continue
        end
        ident = DataFrame(Align_Num = 0, DB_Num = 0, minRI = 0., maxRI = 0., DB_RI = 0., dot = 0., BaseMass = 0., Origin = "", Class = "",  Align_Name = "", DB_Name = "", DB_InChIKey = "", DB_CAS = "", DB_Formula = "", Spec = "", DB_spec = "", )

        for s in ind_s
            # Calculate dot product
            spdb = getSpec(db[s,"Spectrum"])
            master_vect_f = master_mass_gen(sp[:, 1], sp[:, 2], spdb[:, 1], spdb[:, 2], mz_tol)
            # direct_m, revers_m = dot_prod(master_vect_f)
            # vals[v] = direct_m
            if size(master_vect_f, 1) .> numfrags
                master_vect_fu = master_vect_f[sortperm(vec(sum(master_vect_f[:, 3], dims = 2)), rev = true), :][1:numfrags, :]
                master_vect_fr = master_vect_f[sortperm(vec(sum(master_vect_f[:, 2], dims = 2)), rev = true), :][1:numfrags, :]
                master_vect_f = [master_vect_fu; master_vect_fr]
                master_vect_f = unique(master_vect_f, dims = 1)
            end
            # Mean center each vector or not based on method
            if similarity_method == "DISCO"
                # Mean center intensities
                master_vect_f[:, 2:3] = master_vect_f[:, 2:3] .- mean(master_vect_f[:, 2:3], dims = 1)
                p = abs(dot(master_vect_f[:, 2], master_vect_f[:, 3]) / (norm(master_vect_f[:, 2] * norm(master_vect_f[:, 3]))))
            elseif similarity_method == "NDP"
                # Do not mean-center; just normalize directly
                p = dot(master_vect_f[:, 2], master_vect_f[:, 3]) / (norm(master_vect_f[:, 2]) * norm(master_vect_f[:, 3]))
            else
                error("Unknown similarity method: $similarity_method. Use \"DISCO\" or \"NDP\".")
            end
            # master_vect_f[:, 2:3] = master_vect_f[:, 2:3] .- mean(master_vect_f[:, 2:3], dims = 1)
            # p = abs(dot(master_vect_f[:, 2], master_vect_f[:, 3]) / (norm(master_vect_f[:, 2] * norm(master_vect_f[:, 3]))))

            ident = push!(ident,(df[i, "Nr"], db[s, "#"], minRI[s], maxRI[s], db[s, "RImatching"], p, df[i, BaseMassKey], df[i, OriginKey], df[i, ClassKey], df[i, namekey], db[s, "Name"], db[s, "InChIKey"], string(db[s, "CAS Number"]), db[s, "Formula"], spec[i], db[s, "Spectrum"]), promote = true)
        end
        
        # Only keep top 3 matches
        ident = ident[sortperm(ident[:, "dot"], rev = true), :]
        if size(ident, 1) > 3
            ident = ident[1:3, :]
        end
        # Add to all identifications
        identAll = [identAll; ident]
    end

    identAll = identAll[2:end, :]

    if specID == "BM"
        CSV.write(pathFile[1:end-4] * "_LibSearch_BM.csv", identAll)
    elseif specID == "CS"
        CSV.write(pathFile[1:end-4] * "_LibSearch_CS.csv", identAll)
    end
    println("-"^150)
    println("Done library/database searching.")
    println("Start time: ", time_b4_search)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to merge RI values to the leftmost non-missing value 
# and pair them with corresponding RI tolerance values
function mergeToLeft(ri, addRIwin)
    RIval = fill(NaN, size(ri, 1))
    RItol = fill(NaN, size(ri, 1))
    for r = 1:size(ri, 1)
        if ismissing(ri[r, 1])
            if !ismissing(ri[r, 2])
                RIval[r] = ri[r, 2]
                RItol[r] = addRIwin[2]
            elseif !ismissing(ri[r, 3])
                RIval[r] = ri[r, 3]
                RItol[r] = addRIwin[3]
            elseif !ismissing(ri[r, 4])
                RIval[r] = ri[r, 4]
                RItol[r] = addRIwin[4]
            end
        else
            RIval[r] = ri[r, 1]
            RItol[r] = addRIwin[1]
        end
    end
    RIval = float.(RIval)
    RIval[ismissing.(RIval)] .= NaN

    return RIval, RItol
end



# A function to create head-to-tail graphs for user vs. library spectra
function libraryVisualization(pathFile, index)
    time_b4_vis = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Creating graphs... It will take a while... ")
   
    df = CSV.read(pathFile, DataFrame)
    df[!,"Location"] = fill("NA", size(df, 1))

    # Debugging: delete rows where 'Spec' column contains missing values
    dropmissing!(df, :Spec)

    if typeof(index) == String
        index = 1:size(df, 1)
    end

    # Create a directory for the graphs
    pathout = pathFile[1:end-4]
    if !isdir(pathout)
        mkdir(pathout)
    end

    # Keep only printable ASCII (32..126), drop ALL control chars like \r \n \t
    printable_ascii(c) = (' ' <= c <= '~')
    function ascii_printable_only(s::AbstractString)
        io = IOBuffer()
        for c in String(s)
            if printable_ascii(c)
                write(io, c)
            end
        end
        return String(take!(io))
    end

    # Targeted replacements for common unicode science symbols BEFORE ASCII-only filter
    function replace_common_unicode(s::AbstractString)
        str = String(s)
        # common mappings
        pairs = [
            ("β","beta"), ("α","alpha"), ("γ","gamma"), ("δ","delta"),
            ("µ","mu"),   ("ß","beta"), ("a-","alpha-"), ("-a-","-alpha-"),
            ("–","-"), ("—","-"),
            ("’","'"), ("‘","'"), ("“","\""), ("”","\"")
        ]
        for (sym, repl) in pairs
            str = replace(str, string(sym) => repl)
        end
        return str
    end

    # Final sanitizer: map known unicodes, strip control chars, drop non-ASCII
    sanitize_for_json(s::AbstractString) = ascii_printable_only(replace_common_unicode(s))

    # Helper to wrap long text every n characters (works safely after ASCII sanitizing)
    function wrap_text(s::String, n::Int=20)
        s = sanitize_for_json(s)
        return isempty(s) ? s : join([s[i:min(i+n-1,end)] for i in 1:n:length(s)], "<br>")
    end

    for i in index
        us = getSpec(df[i, "Spec"])
        db = getSpec(df[i, "DB_spec"])

        # Normalize intensities; guard against zero/NaN maxima (avoid NaNs in JSON)
        usmax = maximum(us[:,2])
        dbmax = maximum(db[:,2])
        us_y = (isfinite(usmax) && usmax > 0) ? (us[:, 2] ./ usmax) .* 100 : zeros(size(us, 1))
        db_y = (isfinite(dbmax) && dbmax > 0) ? (db[:, 2] ./ dbmax) .* -100 : zeros(size(db, 1))

        # Base head-to-tail spectrum
        trace_user = bar(
            x=us[:, 1], y=us_y, name="User spectrum",
            marker_color="rgb(54,96,146)", width=1.0, offset=0.0
        )

        trace_db = bar(
            x=db[:, 1], y=db_y, name="Database spectrum",
            marker_color="rgb(150,54,52)", width=1.0, offset=0.0
        )

        # Annotate peaks: user spectrum
        ann = []
        bit = ones(size(us, 1))
        while any(bit .== 1)
            ind = argmax(us[:, 2] .* bit)
            push!(ann, attr(x=us[ind, 1], y=us_y[ind]+8,
                            text=string(Int(us[ind, 1])), showarrow=false, font=attr(size=9)))
            bit[us[ind, 1]-5 .< us[:, 1] .<= us[ind, 1]+5] .= 0
        end

        # Annotate peaks: database spectrum
        bit = ones(size(db, 1))
        while any(bit .== 1)
            ind = argmax(db[:, 2] .* bit)
            push!(ann, attr(x=db[ind, 1], y=db_y[ind]-8,
                            text=string(Int(db[ind, 1])), showarrow=false, font=attr(size=9)))
            bit[db[ind, 1]-5 .< db[:, 1] .<= db[ind, 1]+5] .= 0
        end

        # Sanitize feature names
        user_name = wrap_text(string(df[i, "Align_Name"]), 20)
        db_name = wrap_text(string(df[i, "DB_Name"]),   20)

        # Feature name annotations outside the plot area, relative to the paper
        push!(ann, attr(
            x=1.005, # slightly right of plot area (in paper coordinates)
            y=0.75, # top spectrum, relative to 0-1 paper coordinates
            xref="paper",
            yref="paper",
            text="User: "*user_name,
            showarrow=false,
            xanchor="left",
            yanchor="middle",
            textangle=-90,
            font=attr(size=14, color="rgb(54,96,146)")
        ))
        push!(ann, attr(
            x=1.005, # slightly right of plot area
            y=0.25, # bottom spectrum, relative to 0-1 paper coordinates
            xref="paper",
            yref="paper",
            text="DB: "*db_name,
            showarrow=false,
            xanchor="left",
            yanchor="middle",
            textangle=-90,
            font=attr(size=14, color="rgb(150,54,52)")
        ))

        # Layout with secondary y-axis for feature names
        layout = Layout(
            title=attr(text="Head-to-tail spectrum comparison", x=0.5), # Centered title
            xaxis=attr(title="m/z"),
            yaxis=attr(title="Relative intensity", range=[-110,110], zeroline=true, zerolinecolor="black"),
            bargap=0.05,
            showlegend=false,
            annotations=ann,
            margin=attr(r=100)  # Increase r for more white space on the right
        )

        fig = Plot([trace_user, trace_db], layout)

        # Save as PNG
        fname = "$(i)_User_$(df[i,"Align_Num"])=DB_$(df[i,"DB_Num"]).png"
        fpath = joinpath(pathout, fname)
        try
            savefig(fig, fpath)  # PNG via Kaleido
            df[i,"Location"] = "file:///" * fpath
        catch err
            # Retry with annotations removed (if text still problematic)
            try
                layout_no_ann = deepcopy(layout); layout_no_ann[:annotations] = Any[]
                fig2 = Plot([trace_user, trace_db], layout_no_ann)
                savefig(fig2, fpath)
                df[i,"Location"] = "file:///" * fpath
                @warn "Saved without name annotations for entry $(i) due to text encoding."
            catch err2
                # Final fallback: save HTML (always works) so you at least keep a visual
                htmlpath = replace(fpath, ".png" => ".html")
                PlotlyJS.savehtml(fig, htmlpath)
                df[i,"Location"] = "file:///" * htmlpath
                @warn "PNG export failed for entry $(i). Saved HTML instead at $(htmlpath)."
            end
        end
    end

    filtered_pathFile = pathFile[1:end-4] * "_NoNA.csv"
    CSV.write(filtered_pathFile, df)
    println("\n", "Done creating graphs.")
    println("Start time: ", time_b4_vis)
    println("End time: ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)
end



# A function to resort columns in a specific order and delete specified columns from the aligned table(s)
function resort_delete_cols(path2files::AbstractString,
                     fullfilepath::AbstractString,
                     col_names::Vector{String},
                     meta_order::Vector{String};
                     suffix::String = "_Final")

    time_b4_drop = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    println("Deleting specified columns and resorting... ")
    println(" ")

    # Validate input file exists
    if !isfile(fullfilepath)
        error("File not found: $fullfilepath")
    end

    # Load input CSV
    df = CSV.read(fullfilepath, DataFrame)

    println("Columns in file:")
    println(names(df))

    # Apply metadata column ordering (based only on columns that exist in df)
    df_names = String.(names(df))

    # Keep only metadata cols that exist in the dataframe
    existing_meta = [m for m in meta_order if m in df_names]

    # Sample columns (everything else)
    sample_cols = [c for c in df_names if !(c in meta_order)]

    # Final column order
    final_order = vcat(existing_meta, sample_cols)

    # Apply safe reordering
    df = df[:, Symbol.(final_order)]


    # Drop requested columns (only exact matches)
    println("\nSpecified columns to drop:")
    println(col_names)

    # Find which exist
    matched_strings = [nm for nm in df_names if nm in col_names]

    println("\nMatched & will be dropped (exact matches):")
    println(matched_strings)

    existing_drop = Symbol.(matched_strings)

    df_out = select(df, Not(existing_drop))

    # Write output file with suffix
    filename_only = basename(fullfilepath)

    if endswith(filename_only, ".csv")
        out_file = joinpath(
            path2files,
            replace(filename_only, ".csv" => string(suffix, ".csv"))
        )
    else
        out_file = joinpath(path2files, filename_only * suffix * ".csv")
    end

    CSV.write(out_file, df_out)
    
    println("\nDone. Saved processed file: $out_file")
    println("Start time: ", time_b4_drop)
    println("End time:   ", Dates.format(now(), "yyyy-mm-dd HH:MM:SS"))
    println("*"^150)

    return out_file
end

# End of module
end
