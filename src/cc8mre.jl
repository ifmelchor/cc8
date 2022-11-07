#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022

using Distributed

function cc8mre_run(data::Array{Float64}, xStaUTM::Vector{Float64}, yStaUTM::Vector{Float64}, pmax::Vector{Float64}, pinc::Vector{Float64}, fsem::Integer, lwin::Integer, nwin::Integer, nadv::Float64, ccerr::Float64, nini::Integer, fname::String, supfile::Union{String,Nothing})

  # define baseparams
  base = BaseParams(ccerr, fsem, lwin, nwin, nini, nadv, pmax, pinc)

  # define xysta vector
  nsta = length(xStaUTM)
  xysta = [xySta(xStaUTM[s], yStaUTM[s]) for s in 1:nsta]

  # compute
  _run(xysta, data, base, fname, supfile)  

end


function _idx(ip::Integer, nite::Integer, px0::Float64, py0::Float64, base::BaseParams)
  pxx = Array{Float64}(undef, nite)
  pyy = Array{Float64}(undef, nite)
  
  for i in 1:nite
    # for x
    px = px0 - base.pmax[ip] + base.pinc[ip]*(i-1)
    pxx[i] = base.pinc[ip] * floor(Int64, px/base.pinc[ip])
    # for y
    py = py0 - base.pmax[ip] + base.pinc[ip]*(i-1)
    pyy[i] = base.pinc[ip] * floor(Int64, py/base.pinc[ip])
  end

  return [(px, py) for px in pxx for py in pyy]
end


function _shift(ninikk::Integer, data::Array{Float64}, t::Vector{Float64}, base::BaseParams)
  nsta = size(data, 1)
  shiftedtraces = Array{Float64}(undef, nsta, base.lwin)
  
  for n in 1:nsta
    m = ninikk + floor(Int64, base.fsem * t[n])
    shiftedtraces[n, :] = @view data[n, 1+m:base.lwin+m]
  end

  return shiftedtraces
end


function _ccorr(cal::Array{Float64}, base::BaseParams)
  nsta = size(cal, 1)
  cc = Array{Float64}(undef, nsta, nsta)

  for ii in 1:nsta, jj in ii:nsta
    cc[ii,jj] = 0.
    for l in 1:base.lwin
        cc[ii,jj] = cc[ii,jj] + cal[ii,l]*cal[jj,l]
    end
  end

  return cc
end


function _ccorrcoef(cc::Array{Float64})
  nsta = size(cc, 1)
  suma::Float64 = 0

  for ii in 1:nsta-1, jj in ii+1:nsta
    suma = suma + cc[ii,jj] / sqrt(cc[ii,ii]*cc[jj,jj])
  end

  return (2*suma + nsta) / nsta^2
end


function _reshape(mat::Vector{Float64})
  nite = convert(Int16, sqrt(size(mat, 1)))
  smap = Array{Float64}(undef, nite, nite)
  n::Int = 1
  
  for i in 1:nite, j in 1:nite
    smap[i,j] = mat[n]
    n = n + 1
  end
  
  return smap
end
 

function _run(station_list::Vector{xySta}, data::Array{Float64}, base::BaseParams, fname::String, supfile::Union{String,Nothing})

  # compute center location of the array
  xref, yref = refsta(station_list)
  
  # define time_delay and croscorr partial funcions
  time_delay((px, py)) = [px*(sta.x-xref) + py*(sta.y-yref) for sta in station_list]
  zlcc = cal -> _ccorr(cal, base)
  
  # iterate over time
  for n in 1:base.nwin
    @time begin
    nini_adv = floor(Int64, base.nadv*base.lwin) * (n-1)
    ninikk = base.nini + nini_adv
    
    # define partial function for time delay
    shift_delaymap = tdelay -> _shift(ninikk, data, tdelay, base)

    # initialise variables
    px0::Float64 = 0.
    py0::Float64 = 0.

    # iterate over slowness domain
    for ip in 1:length(base.pmax)
      # compute the number of iterations in slowness
      nite = 2*floor(Int64, base.pmax[ip]/base.pinc[ip])+1
      
      # define the slowness intervals
      idx = _idx(ip, nite, px0, py0, base)

      # time delay between stations and ref for all slowness vector (px, py)
      delay_map = pmap(time_delay, idx)
      
      # shift the traces to align them in time
      shiftmap = pmap(shift_delaymap, delay_map)
      # shiftmap is a nite^2 vector of arrays with seismic data. Each vector contain an array of dim (nsta x lwin)

      # calculate the zero-lag cross-correlation of the array
      ccmap = pmap(zlcc, shiftmap)

      # compute the array-averaged cross-correlation coeficient
      sumap = pmap(_ccorrcoef, ccmap)

      # determine MACC (Maximum Array-averaged Cross-Correlation)
      find_ccmax = findmax(sumap)
      ccmax = find_ccmax[1]
      px0 = idx[find_ccmax[2]][1]
      py0 = idx[find_ccmax[2]][2]
      slow, azm = r2p(-px0, -py0)

      # reshape sumap into a nite x nite matrix
      # sumap = _reshape(sumap)
      sumap = reshape(sumap, (nite, nite))

      # get error bounds
      bds = bm2(sumap, nite, base.pmax[ip], base.pinc[ip], ccmax, base.ccerr)

      # save data to hdf5 file
      # save_main_hdf(fname, ip; slow=slow, baz=azm, ccmax=ccmax, bounds=bds)

      if typeof(supfile) <: String
        # save sumap to hdf5 file
        save_sup_hdf(supfile, ip; sumap=sumap)
      end
    end
    end
  end

end
