# using Pkg
# Pkg.instantiate()

using MAGEMinCalc
using Plots

bulk = [50, 15, 10, 10, 10, 0.5, 1.5, 1.5, 8.5, 0.2, 0.3]

T_Liq = MAGEMinCalc.findliq(bulk = bulk, P_kbar = 7.0, T_start_C = 1400.0,
                fo2_buffer = "qfm", fo2_offset = -2.0, Model = "ig")

T_C = collect(T_Liq:-2.0:1000.0)
P_kbar = fill(7.0, length(T_C))

out = MAGEMinCalc.path(bulk = bulk, T_C = T_C, P_kbar = P_kbar, frac_xtal = true, 
                    Model = "ig", fo2_buffer = "qfm", fo2_offset = -2.0)

plot(out["liq1"]["MgO"].*100, out["liq1"]["Al2O3"].*100, xlabel = "MgO (wt%)", ylabel = "Al2O3 (wt%)")