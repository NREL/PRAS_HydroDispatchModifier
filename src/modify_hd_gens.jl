using DataFrames,CSV
using PRAS
using Revise

using Dates
using DelimitedFiles
using Statistics
using TimeZones

include("utils.jl")

"""
`hydro_flex` can be one of either:
  - :high (represent reservoir units as GeneratorStorages)
  - :mid (same as above, but enforce min gen via capacity, energy,
          and regional net load adjustment)
  - :low (dispatch hydro at a fixed level consistent with energy and power constraints)

"""
function modify_hydro(sys::PRAS.SystemModel{N,L,T,P,E},
                        hydro_data::DataFrame,
                        hydro_flex) where {N,L,T,P,E}
    meta = (N = N, L = L, T = T, P = P, E = E)
    
    regions = deepcopy(sys.regions)
    generators = deepcopy(sys.generators)
    region_gen_idxs = deepcopy(sys.region_gen_idxs)
    
    basegens, new_region_gen_idxs, generatorstorages, region_genstor_idxs = loadhydro!(
        regions, generators, region_gen_idxs, sys.timestamps, hydro_data, hydro_flex, meta)
    
    new_sys = SystemModel(regions, sys.interfaces, basegens, new_region_gen_idxs,
            sys.storages, sys.region_stor_idxs, generatorstorages, region_genstor_idxs,
            sys.lines, sys.interface_line_idxs, sys.timestamps)

    return new_sys
end