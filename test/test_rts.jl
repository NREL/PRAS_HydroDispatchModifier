using PRAS_HydroDispatchModifier
using DataFrames,CSV
using PRAS

all_hydros = DataFrame(CSV.File("rts_hdgens.csv"))
rts_sys = PRAS.rts_gmlc()
hd_gens_rts = rts_sys.generators.names[findall(rts_sys.generators.categories .== "Hydro")]

@testset "Medium flexibility dispatch" begin
    dispatch_modifier = :mid
    new_sys = modify_hydro(rts_sys,all_hydros,dispatch_modifier)
    
    hd_gens_new = new_sys.generators.names[findall(new_sys.generators.categories .== "Hydro")]
    num_gens_dataset = length(unique(all_hydros[:,"Generator Name"]))
    
    
    # Test number of generators removed
    @test length(hd_gens_new) == length(hd_gens_rts) - num_gens_dataset 

    mingen = sum(unique(filter(x->x["DataTypeName"] ∈ ["MinGen"] ,all_hydros)[:,["Generator Name","value"]]).value)
    maxgen = sum(unique(filter(x->x["DataTypeName"] ∈ ["MaxGen"] ,all_hydros)[:,["Generator Name","value"]]).value)

    inflow_dset = sum(filter(x->x["DataTypeName"] ∈ ["WeeklyEnergy","MonthlyEnergy"],all_hydros).value)
    
    discharge_cap_newsys = sum(unique(new_sys.generatorstorages.discharge_capacity,dims=2))
    inflow_newsys = sum(new_sys.generatorstorages.inflow)
    
    # Test load adjustment
    @test sum(rts_sys.regions.load) - sum(new_sys.regions.load) == mingen*length(new_sys.timestamps) 
    
    # Test discharge_capacity = maxgen - mingen
    @test discharge_cap_newsys == maxgen - mingen 
    
    # Test inflow = available_energy - load adjustment
    # Test broken for now due to imperfect inflow adjustment (53 weeks =/= 366 days)
    @test_broken inflow_newsys == inflow_dset - mingen*length(new_sys.timestamps) 

end

@testset "High flexibility dispatch" begin
    dispatch_modifier = :high
    new_sys = modify_hydro(rts_sys,all_hydros,dispatch_modifier)
    
    hd_gens_new = new_sys.generators.names[findall(new_sys.generators.categories .== "Hydro")]
    num_gens_dataset = length(unique(all_hydros[:,"Generator Name"]))
    
    
    # Test number of generators removed
    # @test length(hd_gens_new) == length(hd_gens_rts) - num_gens_dataset

    maxgen = sum(unique(filter(x->x["DataTypeName"] ∈ ["MaxGen"] ,all_hydros)[:,["Generator Name","value"]]).value)
    inflow_dset = sum(filter(x->x["DataTypeName"] ∈ ["WeeklyEnergy","MonthlyEnergy"],all_hydros).value)
    
    discharge_cap_newsys = sum(unique(new_sys.generatorstorages.discharge_capacity,dims=2))
    inflow_newsys = sum(new_sys.generatorstorages.inflow)

    
    # Test load unadjusted
    @test sum(rts_sys.regions.load) == sum(new_sys.regions.load) 
    
    # Test discharge_capacity = maxgen - mingen
    @test discharge_cap_newsys == maxgen 
    
    # Test Total inflow = total load
    @test inflow_newsys == inflow_dset 

end
