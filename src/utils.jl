#!/usr/local/bin julia
# coding=utf-8

# Utility functions for cc8mre.jl

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022

# include("types.jl")

"""
  r2p(x, y)
    
    Get slowness and back-azimuth angles
    
"""

function r2p(x::Float64, y::Float64)

  slowness = sqrt(x^2+y^2)
  azimuth = 180.
  
  rad = 57.2957795 
  
  if y < 0

    azimuth = 180. + rad * atan(x/y)
  
  elseif y > 0

    azimuth = rad * atan(x/y)

    if x < 0
      azimuth += 360.
    end
  
  else # y == 0

    if x > 0
      azimuth = 90.
    
    elseif x < 0
      azimuth = 270.
    
    else # x == 0
      azimuth = 666.
    
    end
  
  end
  
  return (slowness, azimuth)
end


function bm2(msum::Array{Float64}, pmax::Float64, pinc::Float64, ccmax::Float64, ccerr::Float64)

  nite = size(msum, 1)
  bnd = Bounds(666., -1., 666., -1.)
  q = Array{Float64}(undef, nite, nite)

  ccth = ccmax - ccerr

  for i in 1:nite
    px = pinc * (i-1) - pmax  
    
    for j in 1:nite 
      py = pinc * (j-1) - pmax 
      
      ( (px == 0) && (py == 0) ) && continue # skip for

      if msum[i,j] >= ccth
        q[i,j] = 1
        
        for x in (-px+pinc, -px-pinc)
          for y in (-py+pinc, -py-pinc)
            s, a = r2p(x, y)
            if ( s > bnd.slomax ) ; bnd.slomax = s end
            if ( s < bnd.slomin ) ; bnd.slomin = s end
            if ( a > bnd.azimax ) ; bnd.azimax = a end
            if ( a < bnd.azimin ) ; bnd.azimin = a end
          end
        end
      
      else
        q[i,j] = 0
      
      end

    end

  end

  if (bnd.azimax > 355) && (bnd.azimin < 5)
    bnd.azimin = 666.
    bnd.azimax = 1.
    
    for i in 1:nite
      px = pinc * (i-1) - pmax 

      for j in 1:nite
        py = pinc * (j-1) - pmax

        ( (px == 0) && (py == 0) ) && continue # skip for
        
        if convert(Bool, q[i,j])
          for x in [-px+pinc, -px-pinc]
            for y in [-py+pinc, -py-pinc]
              s, a = r2p(x, y)
              if ( x > 0 ) && ( a > bnd.azimax ) ; bnd.azimax = a end
              if ( x < 0 ) && ( a < bnd.azimin ) ; bnd.azimin = a end
            end
          end
        end

      end
    end

  end

  return bnd
end


function refsta(station_list::Vector{xySta})
  nsta = length(station_list)
  xref = sum([sta.x for sta in station_list])
  yref = sum([sta.y for sta in station_list])
  
  return (xref/nsta, yref/nsta)
end


