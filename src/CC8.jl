#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022

__precompile__()

module CC8

    include("types.jl")
        export Bounds, xySta, BaseParams
    
    include("utils.jl")
        export r2p, bm2, refsta
    
    include("hdfutils.jl")
        export save_main_hdf, save_sup_hdf

    include("cc8mre.jl")
        export cc8mre_run


end