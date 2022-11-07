#!/usr/local/bin julia
# coding=utf-8

# GNU GPL v2 licenced to I. Melchor and J. Almendros 08/2022


mutable struct Bounds
  azimin :: Float64
  azimax :: Float64
  slomin :: Float64
  slomax :: Float64
end


struct xySta
  x :: Float64
  y :: Float64
end


struct BaseParams
  ccerr :: Float64           #  --> correlation level
  fsem  :: Int64             #  --> sampling rate
  lwin  :: Int64             #  --> time window length
  nwin  :: Int64             #  --> number of time windows
  nini  :: Int64             #  --> first sample
  nadv  :: Float64           #  --> percentage of time advance (0--1)
  pmax  :: Vector{Float64}   #  --> maximum slownes for the grid
  pinc  :: Vector{Float64}   #  --> slownes interval for the grid
end
