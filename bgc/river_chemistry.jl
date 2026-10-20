# River chemistry for the BGC run's river loads, from RC4USCoast (Gomez et al. 2023, Earth System Science Data 15, 2223,
# doi:10.5194/essd-15-2223-2023; NOAA NCEI accession 0260455, doi:10.25921/9jfw-ph50): monthly climatologies (1990-2022)
# of USGS river chemistry and discharge for 140 US rivers, files mclim_19902022_chem.nc and mclim_19902022_disc.nc in
# MAB_RC4_DIR. Included by mab_bgc.jl.
#
# The river cells of the model are in classes (the nearest hand-placed mouth, see river_class in mab_bgc.jl). Each class
# takes the discharge-weighted monthly mean of the chemistry of the RC4USCoast rivers it stands for:
#   Delaware   the Delaware (USGS station at Trenton) and the Schuylkill: the GloFAS cell of the Delaware lies below their
#              confluence, so its discharge includes the Schuylkill
#   Hudson     no station in RC4USCoast: a PROXY, the mean of the neighbouring Connecticut, Housatonic, Delaware and
#              Schuylkill (the published ">50 µM" nitrate of the Hudson is consistent with it)
#   other      all the other rivers with a mouth inside the model box
# Model tracer ← RC4USCoast variable (mmol/m³; Alk meq/m³): NO₃ ← no3, NH₄ ← nh4, DOM ← don (dissolved organic nitrogen),
# Alk ← alk, DIC ← dic (calculated from Alk, pH and temperature), O₂ ← do. P, Z and the particles are not in the data (0).
using Statistics

const RC4_NCD = Base.loaded_modules[only(id for id in keys(Base.loaded_modules) if id.name == "NCDatasets")]
const RC4_VARIABLES = (NO₃ = "no3", NH₄ = "nh4", DOM = "don", Alk = "alk", DIC = "dic", O₂ = "do")

# table[tracer][class] is the 12 monthly concentrations of that class. `box` = (λ_min, λ_max, φ_min, φ_max).
function river_chemistry_table(dir; box, classes = (:Delaware, :Hudson, :other))
    chem = RC4_NCD.NCDataset(joinpath(dir, "mclim_19902022_chem.nc"))
    disc = RC4_NCD.NCDataset(joinpath(dir, "mclim_19902022_disc.nc"))
    names = String.(chem["river_name"][:]); mlon = Float64.(chem["mouth_lon"][:]); mlat = Float64.(chem["mouth_lat"][:])
    grab(ds, v) = [ismissing(x) || !isfinite(x) ? NaN : Float64(x) for x in ds[v][:, :]]     # (river, month)
    Q = grab(disc, "disc")
    rivers = Dict(Symbol(:Delaware) => ["Delaware", "Schuylkill"],
                  Symbol(:Hudson)   => ["Connecticut", "Housatonic", "Delaware", "Schuylkill"],
                  Symbol(:other)    => [names[i] for i in eachindex(names)
                                        if box[1] < mlon[i] < box[2] && box[3] < mlat[i] < box[4] && !(names[i] in ("Delaware", "Schuylkill"))])
    index(list) = [findfirst(==(r), names) for r in list]
    table = Dict(tracer => Dict(class => zeros(12) for class in classes) for tracer in keys(RC4_VARIABLES))
    for (tracer, variable) in pairs(RC4_VARIABLES)
        C = grab(chem, variable)
        for class in classes
            idx = index(rivers[class]); any(isnothing, idx) && error("RC4USCoast has no river in $(rivers[class])")
            for m in 1:12
                ok = [i for i in idx if isfinite(C[i, m])]
                isempty(ok) && error("RC4USCoast: no $variable for $class in month $m")
                w = [isfinite(Q[i, m]) ? max(Q[i, m], 0.0) : 0.0 for i in ok]
                table[tracer][class][m] = sum(w) > 0 ? sum(w .* C[ok, m]) / sum(w) : mean(C[ok, m])
            end
        end
    end
    close(chem); close(disc)
    return table, rivers
end
