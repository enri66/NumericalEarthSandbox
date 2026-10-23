#!/bin/bash
. /t0/workdir/enrique/year_run/env/ib_env.sh
export JULIA_PKG_GC_AUTO=false
julia --project=/t0/workdir/enrique/year_run/env -e "using Pkg; Pkg.precompile(); using Oceananigans, NumericalEarth; println(\"LOADED \", pkgversion(Oceananigans), \" \", pkgversion(NumericalEarth))"
