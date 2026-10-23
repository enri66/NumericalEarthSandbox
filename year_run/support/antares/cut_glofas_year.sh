#!/bin/bash
J=$HOME/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/bin/julia
export JULIA_SSL_CA_ROOTS_PATH=$HOME/.julia/juliaup/julia-1.12.7+0.x64.linux.gnu/share/julia/cert.pem
cd ~/cut_env
CACHE_WHAT=glofas CACHE_LON=-76,-64 CACHE_LAT=34,42 CACHE_START=2019-03-25 CACHE_END=2020-04-05 CACHE_OUT=$HOME/Data/mab_year_cut $J --project=. ~/regional_cache_from_archives.jl
echo CUTDONE
