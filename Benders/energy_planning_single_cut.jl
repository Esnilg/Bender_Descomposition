using JuMP
using Gurobi
using CPLEX
using DelimitedFiles
using DataFrames
using CSV

GRB_ENV = Gurobi.Env()

df_simulations = DataFrame(Scenario=[],LB=[],UB=[],Time=[])

nScenarios = 100

include("../DATOS/Conjuntos_ses_main.jl");
include("../DATOS/datos_ses_main.jl");
include("../DATOS/parametros_inciertos.jl");

Mob = union(TECHNOLOGIES_OF_END_USES_CATEGORY["MOBILITY_PASSENGER"],TECHNOLOGIES_OF_END_USES_CATEGORY["MOBILITY_FREIGHT"])
TECHNOLOGIES = setdiff(TECHNOLOGIES,Mob)

Mob2 = ["MOB_PUBLIC","MOB_PRIVATE","MOB_FREIGHT_RAIL","MOB_FREIGHT_ROAD"]
LAYERS = setdiff(LAYERS,Mob2)
END_USES_TYPES = setdiff(END_USES_TYPES,Mob2)

probability = ones(NS) / NS; # probabilities


# FIRST-STAGE MODEL


# CREATE STOCHASTIC MODEL
master = Model(CPLEX.Optimizer);
set_optimizer_attribute(master, "CPX_PARAM_MIPDISPLAY", 2);
set_optimizer_attribute(master, "CPX_PARAM_SCRIND", 0); 

# FIRST-STAGE MODEL

# var I stage

# feasibility primal
@variable(master,F_Mult[i in TECHNOLOGIES] >=0);
@variable(master,Share_Heat_Dhn >=0);
@variable(master,θ >= -1e10);

#sigmaLB
@constraint(master, row_sigmaLB[i=TECHNOLOGIES],  F_Mult[i] - f_min[i] >= 0.0 )
#sigmaUB
@constraint(master, row_sigmaUB[i=TECHNOLOGIES], -F_Mult[i] + f_max[i] >= 0.0 )
#thetaLB
@constraint(master, row_thetaLB, Share_Heat_Dhn - share_heat_dhn_min >=  0.0 )
#thetaUB
@constraint(master, row_thetaUB,-Share_Heat_Dhn + share_heat_dhn_max >=  0.0)
#gamma5
@constraint(master, row_gamma5,F_Mult["GRID"] - (1 +(9400/c_inv["GRID"])*(F_Mult["WIND"]+F_Mult["PV"])/(f_max["WIND"]+f_max["PV"]))== 0.0)
#gamma4
@constraint(master, row_gamma4,-F_Mult["POWER2GAS_1"] + F_Mult["POWER2GAS_3"] >= 0)
#gamma3
@constraint(master, row_gamma3,-F_Mult["POWER2GAS_2"] + F_Mult["POWER2GAS_3"] >= 0)
#gamma2
@constraint(master, row_gamma2,F_Mult["EFFICIENCY"] - 1/(1 + i_rate) == 0.0)
#gamma1
@constraint(master, row_gamma1,F_Mult["NUCLEAR"] == 0.0)
#gamma8
@constraint(master, row_gamma8,F_Mult["DHN"] - sum(F_Mult[j] for j in TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DHN"]) == 0.0)
#gamma10
@constraint(master, row_gamma10,f_max["PUMPED_HYDRO"] * (F_Mult["NEW_HYDRO_DAM"] - f_min["NEW_HYDRO_DAM"])/(f_max["NEW_HYDRO_DAM"] - f_min["NEW_HYDRO_DAM"]) - F_Mult["PUMPED_HYDRO"] >= 0)

@objective(master, Min, sum((tau[i] * c_inv[i] + c_maint[i]) * F_Mult[i] for i=TECHNOLOGIES)  + θ)

UB = 1e10;
LB = -1e10;

global iter = 0

t0 = time();

while(true)

	global iter = iter + 1 

	JuMP.optimize!(master)
	feasible = JuMP.primal_status(master) == MOI.FEASIBLE_POINT

	variable_value = Dict{JuMP.VariableRef, Float64}()

	if feasible
		for vref in values(F_Mult)
			variable_value[vref] = JuMP.value(vref)
		end
		for vref in values(Share_Heat_Dhn)
			variable_value[vref] = JuMP.value(Share_Heat_Dhn)
		end
		for θ in values(θ)
			variable_value[θ] = JuMP.value(θ)
		end
	end

	# SECOND-STAGE MODELS

	global objective_y =  Dict()
	global dual_F = Dict()
	global dual_SH_DHN = Dict()

	for s in 1:NS
		# stochastic block
		submodel = Model(() -> Gurobi.Optimizer(GRB_ENV))
		set_optimizer_attribute(submodel, "TimeLimit", 10000)
		set_optimizer_attribute(submodel, "Presolve", 0)
		set_optimizer_attribute(submodel, "OutputFlag", 0)

		# second-stage variables
		@variable(submodel,F_Mult_t[union(RESOURCES,TECHNOLOGIES),PERIODS] >=0);
		@variable(submodel,Storage_In[i=STORAGE_TECH,l=LAYERS,PERIODS; storage_eff_in[i,l]>0] >=0);
		@variable(submodel,Storage_Out[i=STORAGE_TECH,l=LAYERS,PERIODS; storage_eff_out[i,l]>0] >=0);
		@variable(submodel,Max_Heat_Demand_DHN >= 0);
		@variable(submodel,InSample >= 0);
		@variable(submodel,w[LAYERS,PERIODS] >=0);
		@variable(submodel,F[i=TECHNOLOGIES] >= 0);
		@variable(submodel,SH_DHN >= 0);
		

		# AUXLIAR VARIABLES

		@expression(submodel, Losses[i=END_USES_TYPES,t=PERIODS],sum(layers_in_out[j,i]*F_Mult_t[j,t] for j in setdiff(union(RESOURCES,TECHNOLOGIES),STORAGE_TECH) if layers_in_out[j, i] > 0)*M_loss_coeff[i][s]);

		@expression(submodel, End_Uses[l=LAYERS, t=PERIODS],   
			if l=="ELECTRICITY"
				M_end_uses_input[l][s]/total_time + M_end_uses_input["LIGHTING"][s]*lighting_month[t]/t_op[t]+Losses[l,t]
			elseif l=="HEAT_LOW_T_DHN"
				(M_end_uses_input["HEAT_LOW_T_HW"][s]/total_time + M_end_uses_input["HEAT_LOW_T_SH"][s]*heating_month[t]/t_op[t])*SH_DHN
			elseif l=="HEAT_LOW_T_DECEN"
				(M_end_uses_input["HEAT_LOW_T_HW"][s]/total_time + M_end_uses_input["HEAT_LOW_T_SH"][s]*heating_month[t]/t_op[t])*(1 - SH_DHN)
			elseif l=="HEAT_HIGH_T"
				M_end_uses_input[l][s] / total_time
			else 0.0 end);


		### CONSTRAINTS ###

		#xi
		@constraint(submodel, row_xi[i=TECHNOLOGIES,t=PERIODS],F[i]*M_c_p_t[i,t][s] - F_Mult_t[i,t] >=0 )
		#delta
		@constraint(submodel, row_delta[i=TECHNOLOGIES],F[i]*M_c_p[i][s]*total_time - sum(F_Mult_t[i,t]*t_op[t] for t in PERIODS) >= 0)
		#beta
		EXP = [1]
		@constraint(submodel, row_beta[i=STORAGE_TECH,t=PERIODS],F_Mult_t[i,t] - ((in(t,EXP) ? F_Mult_t[i,length(PERIODS)] : F_Mult_t[i,t-1])  + (sum(Storage_In[i, l, t]*storage_eff_in[i, l] for l in LAYERS if  storage_eff_in[i,l] > 0) - sum(Storage_Out[i, l, t]/storage_eff_out[i,l] for l in LAYERS if storage_eff_out[i,l] > 0))*t_op[t])== 0.0)
		#gamma12
		@constraint(submodel, row_gamma12[i=END_USES_TYPES,j=TECHNOLOGIES_OF_END_USES_TYPE[i]],-fmin_perc[j]*sum(F_Mult_t[j2,t2]*t_op[t2] for j2 in TECHNOLOGIES_OF_END_USES_TYPE[i], t2 in PERIODS) + sum(F_Mult_t[j, t]*t_op[t] for t in PERIODS) >= 0)
		#gamma11
		@constraint(submodel, row_gamma11[i=END_USES_TYPES,j=TECHNOLOGIES_OF_END_USES_TYPE[i]],fmax_perc[j]*sum(F_Mult_t[j2,t2]*t_op[t2] for j2 in TECHNOLOGIES_OF_END_USES_TYPE[i], t2 in PERIODS) - sum(F_Mult_t[j, t]*t_op[t] for t in PERIODS) >= 0)
		#gamma9
		@constraint(submodel, row_gamma9[t=PERIODS],(F_Mult_t["HYDRO_DAM",t] + F_Mult_t["NEW_HYDRO_DAM",t]) - Storage_In["PUMPED_HYDRO","ELECTRICITY",t] >= 0)
		#gamma7
		@constraint(submodel, row_gamma7[t=PERIODS],-(End_Uses["HEAT_LOW_T_DHN",t] + Losses["HEAT_LOW_T_DHN",t] - w["HEAT_LOW_T_DHN",t]) + Max_Heat_Demand_DHN >= 0)
		#gamma6
		@constraint(submodel, row_gamma6,-peak_dhn_factor*Max_Heat_Demand_DHN + F["DHN"] >= 0)
		#lambda
		EXP = ["HEAT_LOW_T_DHN"]
		@constraint(submodel, row_lambda[l=LAYERS,t=PERIODS], sum(layers_in_out[i,l]*F_Mult_t[i,t] for i in setdiff(union(RESOURCES,TECHNOLOGIES),STORAGE_TECH)) + sum(Storage_Out[j,l,t] - Storage_In[j,l,t] for j in STORAGE_TECH if storage_eff_in[j,l] + storage_eff_out[j,l] > 0) - End_Uses[l,t] - (in(l,EXP) ? Losses[l,t] : 0.0) + w[l,t] == 0.0)
		#mu
		@constraint(submodel, row_mu[i=RESOURCES],M_avail[i][s] - sum(F_Mult_t[i,t]*t_op[t] for t in PERIODS) >= 0)
		@constraint(submodel, InSample == sum(M_c_op[i,t][s] * F_Mult_t[i,t] * t_op[t] for t=PERIODS for i=RESOURCES) + 5000*sum(w[i,p] for i=LAYERS for p = PERIODS))
		
		@constraint(submodel, [l=RESO1,t=PERIODS], w[l,t] == 0.0)
		
		@constraint(submodel, dual_T[i=TECHNOLOGIES], F[i] == variable_value[F_Mult[i]])
		@constraint(submodel, dual_DHN, SH_DHN == variable_value[Share_Heat_Dhn])
		
		@objective(submodel, Min, InSample);
		
		optimize!(submodel)
		
		feasible = JuMP.primal_status(submodel) == MOI.FEASIBLE_POINT
	
		#println(feasible)
		
		objective_y[s] = objective_value(submodel)
		dual_F[s] = JuMP.dual.(dual_T)
		dual_SH_DHN[s] = JuMP.dual(dual_DHN)
		#println(JuMP.dual.(dual))
		#println(objective_value(submodel))
	end

	obj = sum((tau[i] * c_inv[i] + c_maint[i]) * variable_value[F_Mult[i]] for i=TECHNOLOGIES)

	UB = obj + sum(probability[s]*objective_y[s] for s=1:NS)

	LB = obj + variable_value[θ]

	if abs(UB-LB) < 0.01
		time_end = time() - t0
		println("\n================== END ==================");
		println("Lower Bound: $(LB)")
		println("Upper Bound: $(UB)")
		println("Gap = ", abs(UB-LB))
		println("Elapsed Time: ", time_end, " seconds")
		break
	else
		@constraint(master, θ >= sum(probability[s]*objective_y[s] for s=1:NS) + sum(probability[s]*dual_F[s][i]*(F_Mult[i] - variable_value[F_Mult[i]]) for i=TECHNOLOGIES for s=1:NS) + sum(probability[s]*dual_SH_DHN[s]*(Share_Heat_Dhn - variable_value[Share_Heat_Dhn]) for s=1:NS))
	end

	time_end = time() - t0

	println("\n================== Iteration $(iter) ==================");
	println("Lower Bound: $(LB)")
	println("Upper Bound: $(UB)")
	println("Gap = ", abs(UB-LB))
	println("Elapsed Time: ", time_end, " seconds")
	
	push!(df_simulations,[nScenarios,LB,UB,time_end])

	#CSV.write("df_Bender_single_cut_$(nScenarios).csv",df_simulations,delim=';', decimal='.');
	
end


