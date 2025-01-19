using Plots
using MAGEMinCalc

comp = Dict("SiO2_Liq" => 47.5, "Al2O3_Liq" => 16.4, "CaO_Liq" => 11.6, "MgO_Liq" => 9.38,
            "FeOt_Liq" => 9.16, "K2O_Liq" => 0.329, "Na2O_Liq" => 2.25, "TiO2_Liq" => 2.29, 
            "Fe3Fet_Liq" => 0.15, "Cr2O3_Liq" => 0.0, "H2O_Liq" => 0.68)

Results = MAGEMinCalc.path(comp = comp, T_end_C = 1100.0, dt_C = 2.0, 
            P_bar = 1000.0, frac_xtal = true, 
            Model = "ig",
            find_liquidus = true)

# plot(out["liq1"]["MgO"].*100, out["liq1"]["Al2O3"].*100, xlabel = "MgO (wt%)", ylabel = "Al2O3 (wt%)")