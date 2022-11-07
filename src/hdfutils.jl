#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022

using HDF5

function save_main_hdf(fname::String, ip::Integer; kwargs...)

    # open hdf5 and read the step process bin (stb)
    fid = h5open(fname, "r+")
    hdr = fid["header"]
    stb = read_attribute(hdr, "last_time_bin")
    
    # save tseries
    tseries = fid[string(ip-1)]
    tseries["slow"][stb+2] = kwargs[:slow]
    tseries["baz"][stb+2] = kwargs[:baz]
    tseries["cc_max"][stb+2] = kwargs[:ccmax]

    # save bounds
    bds = tseries["bounds"]
    bds["baz_min"][stb+2] = kwargs[:bounds].azimin
    bds["baz_max"][stb+2] = kwargs[:bounds].azimax
    bds["slo_min"][stb+2] = kwargs[:bounds].slomin
    bds["slo_max"][stb+2] = kwargs[:bounds].slomax

    # and add +1 to step process bin
    attrs(hdr)["last_time_bin"] = stb + 1

    # close hdf5
    close(fid)

end

function save_sup_hdf(sfname::String, ip::Integer; kwargs...)

    # open hdf5
    fname = string(sfname, ".", ip-1)
    fid = h5open(fname, "r+")

    # read the step process bin (stb)
    dset = fid["sumap"]
    stb = read_attribute(dset, "last_time_bin")

    # save tseries
    dset[stb+2,:,:] = kwargs[:sumap]

    # and add +1 to step process bin
    attrs(dset)["last_time_bin"] = stb+1

    # close hdf5
    close(fid)
    
end