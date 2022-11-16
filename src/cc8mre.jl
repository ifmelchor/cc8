#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022

using LinearAlgebra

function run(data::Array{Float64}, xStaUTM::Vector{Float64}, yStaUTM::Vector{Float64}, pmax::Vector{Float64}, pinc::Vector{Float64}, fsem::Integer, lwin::Integer, nwin::Integer, padv::Float64, ccerr::Float64, nini::Integer)

    # define baseparams
    nsta = length(xStaUTM)
    xysta = [xySta(xStaUTM[s], yStaUTM[s]) for s in 1:nsta]
    refxy = refsta(xysta)

    cciter = []
    for ii in 1:nsta-1
        for jj in ii+1:nsta
            push!(cciter, (ii, jj))
        end
    end

    nadv = floor(Int64, padv*lwin)
    base = BaseParams(data, xysta, refxy[1], refxy[2], ccerr, fsem, lwin, nwin, nini, padv, nadv, pmax, pinc, cciter)

    res = _run(base)

    return res
end


function _idx(ip::Integer, nite::Integer, pxy0::Tuple{Float64, Float64}, base::BaseParams)
    pxx = Array{Float64}(undef, nite)
    pyy = Array{Float64}(undef, nite)
    for i in 1:nite
    # for x
    px = pxy0[1] - base.pmax[ip] + base.pinc[ip]*(i-1)
    pxx[i] = base.pinc[ip] * floor(Int64, px/base.pinc[ip])
    # for y
    py = pxy0[2] - base.pmax[ip] + base.pinc[ip]*(i-1)
    pyy[i] = base.pinc[ip] * floor(Int64, py/base.pinc[ip])
    end

    return [[px, py] for px in pxx for py in pyy]
end


function _run(base::BaseParams)
    # initialize variables

    function _pccorr(nkk, pxy)
        # compute delat times for each station
        dtimes = [pxy[1]*(sta.x-base.xref) + pxy[2]*(sta.y-base.yref) for sta in base.stalist]
        nsta = length(base.stalist)

        # build cc matrix
        cc = zeros(Float64, nsta, nsta)
        for ii in 1:nsta
            mii = nkk + floor(Int64, base.fsem * dtimes[ii])
            dii = @view base.data[ii, 1+mii:base.lwin+mii]
            for jj in ii:nsta
                mjj = nkk + floor(Int64, base.fsem * dtimes[jj])
                djj = @view base.data[jj, 1+mjj:base.lwin+mjj]
                cc[ii,jj] += dot(dii,djj)
            end
        end

        # computes crosscorr coefficient
        suma = sum([cc[ii,jj]/sqrt(cc[ii,ii]*cc[jj,jj]) for (ii, jj) in base.citer])
        return (2*suma + nsta) / nsta^2
    end

    pxy0 = (0.,0.)
    procdatalist = []
    for ip in 1:length(base.pmax)  
        nite = 2*floor(Int64, base.pmax[ip]/base.pinc[ip])+1
        pxylist = _idx(ip, nite, pxy0, base)

        # data to store
        maac = Array{Float64}(undef, base.nwin)
        slow = Array{Float64}(undef, base.nwin)
        bazm = Array{Float64}(undef, base.nwin)
        sumap = Array{Float64}(undef, base.nwin, nite, nite)
        
        best_maac::Float64 = -1.
        for nk in 1:base.nwin
            nkk = base.nini + base.nadv*(nk-1)
            ccmap = map(pxyl->_pccorr(nkk, pxyl), pxylist)
            sumap[nk,:,:] = reshape(ccmap, nite, nite)

            # find max
            find_ccmax = findmax(ccmap)
            ccmax = find_ccmax[1]
            maac[nk] = ccmax
            px = pxylist[find_ccmax[2]][1]
            py = pxylist[find_ccmax[2]][2]
            slow[nk], bazm[nk] = r2p(-px, -py)

            # get the best pxy
            if ccmax > best_maac
                best_maac = ccmax
                best_pxy = (px, py)
            end
        end

        # data to save
        push!(procdatalist, ProcessData(ip, slow, bazm, maac, sumap))

        # change px0 py0 and keep computing with next slowness domain
        if best_maac > base.ccerr
          pxy0 = best_pxy
        else
          break
        end
    end

    return procdatalist
end