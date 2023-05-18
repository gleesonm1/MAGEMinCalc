using MAGEMin_C

function satPhase(P_kbar, T_C, bulk)

		gv, DB = init_MAGEMin();
		sys_in = "wt";

		gv.verbose = -1;

		new_bulk = bulk/sum(bulk);
		
		out = point_wise_minimization(P_kbar, T_C, new_bulk, gv, DB, sys_in);
		Phase = out.ph
		Oxides = out.oxides
		if "liq" in Phase
			Liq_index = findfirst(x -> occursin("liq", x), out.ph)
			Liq_Comp = out.SS_vec[Liq_index].Comp_wt
			Liq_Frac = out.ph_frac_wt[Liq_index]
		
			finalize_MAGEMin(gv,DB);
			Ret = Dict("Phase" => Phase, "Oxides" => Oxides, "Liq_Comp" => Liq_Comp, "Liq_Frac" => Liq_Frac)
			return Ret
		else
			finalize_MAGEMin(gv, DB);
			Ret = Dict("Phase" => Phase, "Oxides" => Oxides)
			return Ret
		end
	end

function PT_minimisation(P_kbar, T_C, bulk)

		gv, z_b, DB, splx_data = init_MAGEMin();
		sys_in = "wt";

		gv.verbose = -1;

		new_bulk = bulk/sum(bulk);

		gv = define_bulk_rock(gv, new_bulk);		

		out = point_wise_minimization(P_kbar, T_C, gv, z_b, DB, splx_data, sys_in);
		Phase = out.ph;
		Oxides = out.oxides;
		Type = out.ph_type;

		unique_phases = unique(Phase)
		
		Ret = Dict();
		if length(Phase) > 0		
			Ret["sys"] = Dict("Phase" => Phase, "Oxides" => Oxides, "Comp" => out.bulk, "Entropy" => out.entropy);
			
			i = 0
			j = 0
			for index in 1:len(Phase)
				Frac = out.ph_frac_wt[index];
				if Type[index] == 0
					i = i + 1
					Comp = out.PP_vec[i].Comp_wt;
					Ret[Phase[index]] = Dict("Frac" => Frac, "Comp" => Comp);
				else
					j = j +1
					Comp = out.SS_vec[j].Comp_wt;
					Ret[Phase[index]] = Dict("Frac" => Frac, "Comp" => Comp);
				end
			end
		end

		finalize_MAGEMin(gv, DB);
		return Ret
	end