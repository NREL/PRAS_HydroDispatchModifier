function loadhydro!(
    regions::Regions, gens::PRAS.Generators, region_gen_idxs::Vector{UnitRange{Int64}},
    timesteps::StepRange{ZonedDateTime,Hour}, hydrodata::DataFrame, hydroflex::Symbol, meta::NamedTuple)

    hydroflex == :low &&
    return loadhydro_lowflex!(regions, gens, region_gen_idxs,
                              timesteps, hydrodata, meta)

    n_regions = length(regions)

    hydronames = Set(hydrodata[!, "Generator Name"])
    hydplant_idx = findall(gens.names .∈ Ref(hydronames))
    hydplant_regidx = [findall(id .∈ region_gen_idxs)[1] for id in hydplant_idx]
    hydrogens = DataFrame([gens.names[hydplant_idx]],["Name"])
    hydrogens[!,:RegionIdx] .= hydplant_regidx

    hydrotype = [hydrodata[findfirst(name .== (hydrodata[!,"Generator Name"])),
                            "Resolution"]
                for name in hydrogens.Name]
    hydrogens[!, :HydroType] .= hydrotype

    sort!(hydrogens, :RegionIdx)
    n_gens = nrow(hydrogens)
    hydronames = Set(hydrogens.Name)

    weekenergy, weekmingen, weekmaxgen =
        loadhydrotime(hydrodata, "Weekly")
    monthenergy, monthmingen, monthmaxgen =
        loadhydrotime(hydrodata, "Monthly")    

    # Remove hydro generators from the current generators Array
    newgens,new_reggenidx = remove_hydro_generators!(gens,region_gen_idxs,hydronames,meta)

    fixedhydro_generation = DataFrame(zeros((meta.N,length(hydronames))),[name for name in hydronames])

    region_genstors = searchsorted.(Ref(hydrogens.RegionIdx), 1:n_regions)

    allzero = zeros(Float64, n_gens, meta.N)
    allone = ones(Float64, n_gens, meta.N)

    inflow = zeros(Int, n_gens, meta.N)
    energy_capacity = zeros(Int, n_gens, meta.N)
    discharge_capacity = zeros(Int, n_gens, meta.N)

    for (t, dt) in enumerate(timesteps)
        for (i, row) in enumerate(eachrow(hydrogens))

            # TODO: Prorate inflow and mingen_e when not starting
            #       on first day of week/month

            if row.HydroType == "Weekly"

                w = week(dt)
                w_prev = week(dt - Hour(1))

                # According to PNNL sometimes these values are reversed?
                # ... seems suspicious
                w_mingen, w_maxgen =
                    minmax(weekmingen[w, row.Name], weekmaxgen[w, row.Name])

                w_energy = round(Int, weekenergy[w, row.Name])
                w_capacity = round(Int, w_maxgen)

                if hydroflex == :mid

                    w_capacity -= round(Int, w_mingen)

                    # When mingen values are infeasible, we relax mingen to the
                    # lowest value consistent with the energy budget
                    # Note this leaves the unit with no PRAS-dispatchable energy

                    w_minenergy = round(Int, w_mingen * 24 * 7)

                    if w_energy < w_minenergy
                        fixedhydro_generation[t, row.Name] = w_energy / (24 * 7)
                        w_energy = 0
                    else
                        fixedhydro_generation[t, row.Name] = w_mingen
                        w_energy -= w_minenergy
                    end

                    (w_capacity < 0 || w_energy < 0) && (println(row.Name, "\t", w_mingen, "\t", w); error("Bad hydro data (weekly)"))

                end

                inflow[i,t] = (t == 1 || w != w_prev) ? w_energy : 0
                energy_capacity[i,t] = w_energy
                discharge_capacity[i,t] = w_capacity

            else # Monthly

                m = month(dt)
                m_prev = month(dt - Hour(1))

                m_mingen, m_maxgen =
                    minmax(monthmingen[m, row.Name], monthmaxgen[m, row.Name])

                m_energy = round(Int, monthenergy[m, row.Name])
                m_capacity = round(Int, m_maxgen)

                if hydroflex == :mid

                    # When mingen values are infeasible, we relax mingen to the
                    # lowest value consistent with the energy budget
                    # Note this leaves the unit with no PRAS-dispatchable energy

                    m_energy_fixed = round(Int, m_mingen * 24 * monthlength[m])
                    m_capacity -= round(Int, m_mingen)

                    if m_energy < m_energy_fixed
                        fixedhydro_generation[t, row.Name] = m_energy / (24 * monthlength[m])
                        m_energy = 0
                    else
                        fixedhydro_generation[t, row.Name] = m_mingen
                        m_energy -= m_energy_fixed
                    end

                    (m_capacity < 0 || m_energy < 0) && (println(row.Name, "\t", m_mingen, "\t", m, "\t",m_capacity, "\t", m_energy); error("Bad hydro data (monthly)"))
                end


                inflow[i,t] = (t == 1 || m != m_prev) ? m_energy : 0
                energy_capacity[i,t] = m_energy
                discharge_capacity[i,t] = m_capacity

            end

        end
    end

    hydroflex == :mid && netloadadjustment!(regions, hydrogens, fixedhydro_generation)

    select!(fixedhydro_generation, hydronames...)

    genstors = GeneratorStorages{meta.N, meta.L, meta.T, meta.P, meta.E}(
        hydrogens.Name, hydrogens.HydroType .* "Hydro",
        inflow, discharge_capacity, energy_capacity,
        allone, allone, allone,
        inflow, Int.(allzero), discharge_capacity, allzero, allone)

    return newgens, new_reggenidx, genstors, region_genstors
           #fixedhydro_generation, fixedhydro_inflow, hydro_mingen, hydro_maxgen

end

function loadhydro_lowflex!(
    regions::Regions, gens::PRAS.Generators, region_gen_idxs::Vector{UnitRange{Int64}},
    timesteps::StepRange{ZonedDateTime,Hour}, hydrodata::DataFrame,
    meta)

    n_regions = length(regions)

    hydronames = Set(hydrodata[!, "Generator Name"])
    hydplant_idx = findall(gens.names .∈ Ref(hydronames))
    hydplant_regidx = [findall(id .∈ region_gen_idxs)[1] for id in hydplant_idx]
    hydrogens = DataFrame([gens.names[hydplant_idx]],["Name"])
    hydrogens[!,:RegionIdx] .= hydplant_regidx

    hydrotype = [hydrodata[findfirst(name .== (hydrodata[!,"Generator Name"])),
                            "Resolution"]
                for name in hydrogens.Name]
    hydrogens[!, :HydroType] .= hydrotype

    sort!(hydrogens, :RegionIdx)
    n_gens = nrow(hydrogens)
    hydronames = Set(hydrogens.Name)

    weekenergy, weekmingen, weekmaxgen =
        loadhydrotime(hydrodata, "Weekly")
    monthenergy, monthmingen, monthmaxgen =
        loadhydrotime(hydrodata, "Monthly")    

    # Remove hydro generators from the current generators Array
    newgens,new_reggenidx = remove_hydro_generators!(gens,region_gen_idxs,hydronames,meta)

    fixedhydro_generation = DataFrame(zeros((meta.N,length(hydronames))),[name for name in hydronames])

    for (t, dt) in enumerate(timesteps)
        for (i, row) in enumerate(eachrow(hydrogens))

            # TODO: Prorate inflow and mingen_e when not starting
            #       on first day of week/month

            if row.HydroType == "Weekly"

                w = week(dt)
                w_prev = week(dt - Hour(1))

                # According to PNNL sometimes these values are reversed?
                # ... seems suspicious
                w_mingen, w_maxgen =
                    minmax(weekmingen[w, row.Name], weekmaxgen[w, row.Name])

                w_energy = weekenergy[w, row.Name]
                w_fixedgen = min(w_energy / (24 * 7), w_maxgen)

                fixedhydro_generation[t, row.Name] = w_fixedgen

            else # Monthly

                m = month(dt)
                m_prev = month(dt - Hour(1))

                m_mingen, m_maxgen =
                    minmax(monthmingen[m, row.Name], monthmaxgen[m, row.Name])

                m_energy = monthenergy[m, row.Name]
                m_fixedgen = min(m_energy / (24 * monthlength[m]), m_maxgen)

                fixedhydro_generation[t, row.Name] = m_fixedgen

            end

        end
    end

    netloadadjustment!(regions, hydrogens, fixedhydro_generation)


    empty_int = zeros(Int, 0, meta.N)
    empty_float = zeros(Float64, 0, meta.N)
    genstors = GeneratorStorages{meta.N, meta.L, meta.T, meta.P, meta.E}(
        String[], String[],
        empty_int, empty_int, empty_int,
        empty_float, empty_float, empty_float,
        empty_int, empty_int, empty_int, empty_float, empty_float)

    return newgens, new_reggenidx, genstors, fill(1:0, n_regions)
           #fixedhydro_generation, fixedhydro_inflow, hydro_mingen, hydro_maxgen

end

function loadhydrotime(hydrodata::DataFrame,hydrotype::String)

    hydro_typesubset = subset(hydrodata, "Resolution" => ByRow(n -> n == hydrotype))

    energy = timecolumns(hydro_typesubset, hydrotype * "Energy")
    mingen = timecolumns(hydro_typesubset, "MinGen")
    maxgen = timecolumns(hydro_typesubset, ["MaxGen", "MaxCap"])

    return energy, mingen, maxgen

end

function remove_hydro_generators!(gens,reggenidx_init,hydronames,meta)
    
    genidxremove = findall(gens.names .∈ Ref(hydronames))
    adjreggenidx = cumsum([sum(genidxremove .∈ Ref(sublist)) for sublist in reggenidx_init])

    numregions = length(reggenidx_init)
    region_gen_idxs = Array{UnitRange{Int64},1}(undef,numregions); 

    for (idx,reg_range_idx) in enumerate(reggenidx_init)
        if idx == 1 
            region_gen_idxs[idx] = 1:last(reg_range_idx)-adjreggenidx[idx]    
        else
            region_gen_idxs[idx] = range(first(reg_range_idx)-adjreggenidx[idx-1],
                                        last(reg_range_idx)-adjreggenidx[idx])
        end                                    
    end

    nothyd_genidx = findall(gens.names .∉ Ref(hydronames))

    return (Generators{meta.N, meta.L, meta.T, meta.P}(
                gens.names[nothyd_genidx], 
                gens.categories[nothyd_genidx], 
                gens.capacity[nothyd_genidx,:],
                gens.λ[nothyd_genidx,:],
                gens.μ[nothyd_genidx,:]), 
            region_gen_idxs)
end

function netloadadjustment!(
    regions::Regions, gens::DataFrame, generation::DataFrame,
    mod_categories::Vector{String}=String[])

    modifiers = if length(mod_categories) > 0
        subset(gens, :Category => ByRow(c -> c in mod_categories))
    else
        gens
    end
    modifiers = select(modifiers, :Name, :RegionIdx)

    missingunits = 0
    for mod in eachrow(modifiers)

        if !(mod.Name in names(generation))
            #println(mod.Name, " (", mod.Category, ") has no generation data")
            missingunits +=1
            continue
        end

        regions.load[mod.RegionIdx, :] .-=
            round.(Int, generation[!, mod.Name])

    end
    #println(missingunits, " load modifiers with missing data")

end

# Note: Using GridView week definition, not ISO-8601
# TODO provide flags to use different week definitions
week(dt::ZonedDateTime) = div(dayofyear(dt) - 1, 7) + 1

timestack(df::DataFrame) =
    stack(df, Not(["Generator Name", "DataTypeName"]), variable_name="Time")

timecolumns(df::DataFrame, datatypename::String) =
    timecolumns(df, [datatypename])

function timecolumns(df::DataFrame, datatypenames::Vector{String})

    result = df[df.DataTypeName .∈ Ref(datatypenames), :]
    result = unstack(result, "Time", "Generator Name", "value")
    disallowmissing!(result)

    return result

end

monthlength = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

