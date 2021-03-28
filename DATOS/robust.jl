
for i = PERIODS

# max value
c_op["ELECTRICITY",i]=c_op["ELECTRICITY",i]*(1+0.8992)
c_op["GASOLINE",i]   =c_op["GASOLINE",i]*(1+0.8992)
c_op["DIESEL",i]     =c_op["DIESEL",i]*(1+0.8992)
c_op["LFO",i]        =c_op["LFO",i]*(1+0.8992)
c_op["NG",i]         =c_op["NG",i]*(1+0.8992)
c_op["WOOD",i]       =0.0959
c_op["COAL",i]       =c_op["COAL",i]*(1+0.8992)
c_op["URANIUM",i]    =c_op["URANIUM",i]*(1+0.8992)
c_op["NG_CCS",i]     =c_op["NG",i]
c_op["SNG",i]        =c_op["NG",i]
c_op["COAL",i]       =c_op["COAL_CCS",i]
c_op["BIODIESEL",i]  =c_op["DIESEL",i]
c_op["BIOETHANOL",i] =c_op["GASOLINE",i]

# min value
c_p_t["PV",i]              =0.8893*c_p_t["PV",i]
c_p_t["WIND",i]            =0.8893*c_p_t["WIND",i]
c_p_t["HYDRO_DAM",i]       =0.8893*c_p_t["HYDRO_DAM",i]
c_p_t["HYDRO_RIVER",i]     =0.8893*c_p_t["HYDRO_RIVER",i]
c_p_t["NEW_HYDRO_DAM",i]   =0.8893*c_p_t["NEW_HYDRO_DAM",i]
c_p_t["NEW_HYDRO_RIVER",i] =0.8893*c_p_t["NEW_HYDRO_RIVER",i]
c_p_t["DEC_SOLAR",i]       =0.8893*c_p_t["DEC_SOLAR",i]

end

# max value
end_uses_demand_year[("ELECTRICITY","HOUSEHOLDS")]           =11317.4668168738
end_uses_demand_year[("ELECTRICITY","SERVICES")]             =15637.1646491666
end_uses_demand_year[("ELECTRICITY","INDUSTRY")]             =11064.8764607173
end_uses_demand_year[("LIGHTING","HOUSEHOLDS")]              =443.4928829798
end_uses_demand_year[("LIGHTING","SERVICES")]                =3959.8402104954
end_uses_demand_year[("LIGHTING","INDUSTRY")]                =1338.9946733427
end_uses_demand_year[("HEAT_HIGH_T","INDUSTRY")]             =20153.2577773288
end_uses_demand_year[("HEAT_LOW_T_SH","HOUSEHOLDS")]         =30765.1148547813
end_uses_demand_year[("HEAT_LOW_T_SH","SERVICES")]           =15115.0759721968
end_uses_demand_year[("HEAT_LOW_T_SH","INDUSTRY")]           =5241.8706649494
end_uses_demand_year[("HEAT_LOW_T_HW","HOUSEHOLDS")]         =7863.9394338392
end_uses_demand_year[("HEAT_LOW_T_HW","SERVICES")]           =3388.3211724411
end_uses_demand_year[("HEAT_LOW_T_HW","INDUSTRY")]           =1358.0656530232


end_uses_input = Dict(i => sum([end_uses_demand_year[i,s] for s in SECTORS]) for i in END_USES_INPUT);

# min value
c_p["NUCLEAR"]         =0.8282
c_p["CCGT"]            =0.8292
c_p["CCGT_CCS"]        =0.8292
c_p["COAL_US"]         =0.8467
c_p["COAL_IGCC"]       =0.835
c_p["COAL_US_CCS"]     =0.8467
c_p["COAL_IGCC_CCS"]   =0.835
c_p["GEOTHERMAL"]      =0.8389
c_p["IND_COGEN_GAS"]   =0.8292
c_p["IND_COGEN_WOOD"]  =0.8292
c_p["IND_COGEN_WASTE"] =0.8292
c_p["IND_BOILER_GAS"]  =0.9267
c_p["IND_BOILER_WOOD"] =0.878
c_p["IND_BOILER_OIL"]  =0.9267
c_p["IND_BOILER_COAL"] =0.878
c_p["IND_BOILER_WASTE"]=0.878
c_p["IND_DIRECT_ELEC"] =0.9267
c_p["DHN_HP_ELEC"]     =0.9267
c_p["DHN_COGEN_GAS"]   =0.8292
c_p["DHN_COGEN_WOOD"]  =0.8292
c_p["DHN_COGEN_WASTE"] =0.8292
c_p["DHN_BOILER_GAS"]  =0.9267
c_p["DHN_BOILER_WOOD"] =0.878
c_p["DHN_BOILER_OIL"]  =0.9267
c_p["DHN_DEEP_GEO"]    =0.8292
c_p["POWER2GAS_1"]     =0.878
c_p["POWER2GAS_2"]    =0.878
c_p["H2_ELECTROLYSIS"] =0.878
c_p["H2_NG"]           =0.8389
c_p["H2_BIOMASS"]      =0.84
c_p["GASIFICATION_SNG"]=0.83
c_p["PYROLYSIS"]       =0.83

# min value
avail["WOOD"]  =8340.45
avail["WASTE"] =7568.15

# max value
loss_coeff["ELECTRICITY"]   =0.0714
loss_coeff["HEAT_LOW_T_DHN"]=0.051
