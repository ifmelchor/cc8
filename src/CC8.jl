#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022

__precompile__()

module CC8

    export run

    include("types.jl")

    include("utils.jl")

    include("cc8mre.jl")

end