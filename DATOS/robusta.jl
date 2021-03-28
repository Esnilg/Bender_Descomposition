

Lista_cpt = ["PV","WIND","HYDRO_DAM","HYDRO_RIVER","DEC_SOLAR","NEW_HYDRO_DAM","NEW_HYDRO_RIVER"]

for i in Lista_cpt, for t = 1:12
c_p_t[i,t] = 0.8893*c_p_t[i,t]
end

"NUCLEAR","CCGT","CCGT_CCS","COAL_US","COAL_IGCC",
             "COAL_US_CCS","COAL_IGCC_CCS","GEOTHERMAL","IND_COGEN_GAS",
             "IND_COGEN_WOOD","IND_COGEN_WASTE","IND_BOILER_GAS",
             "IND_BOILER_WOOD","IND_BOILER_OIL","IND_BOILER_COAL",
             "IND_BOILER_WASTE","IND_DIRECT_ELEC","DHN_HP_ELEC",
             "DHN_COGEN_GAS","DHN_COGEN_WOOD","DHN_COGEN_WASTE",
             "DHN_BOILER_GAS","DHN_BOILER_WOOD","DHN_BOILER_OIL",
             "DHN_DEEP_GEO","POWER2GAS_1","POWER2GAS_2","H2_ELECTROLYSIS",
             "H2_NG","H2_BIOMASS","GASIFICATION_SNG","PYROLYSIS"
c_p[]

c_inv["NUCLEAR"]=11346.7
c_inv["CCGT"]=1030.49
c_inv["CCGT_CCS"]=1777.9
c_inv["COAL_US"]=3359.4
c_inv["COAL_IGCC"]=4839.99
c_inv["COAL_US_CCS"]=6042.32
c_inv["COAL_IGCC_CCS"]=8440.56
c_inv["PV"]=1396.34
c_inv["WIND"]=1801.55
c_inv["HYDRO_DAM"]=8393.05
c_inv["NEW_HYDRO_DAM"]=5974.64
c_inv["HYDRO_RIVER"]=6550.56
c_inv["NEW_HYDRO_RIVER"]=7197.07
c_inv["GEOTHERMAL"]=18580.33
c_inv["IND_COGEN_GAS"]=1828.26
c_inv["IND_COGEN_WOOD"]=1403.35
c_inv["IND_COGEN_WASTE"]=3801.68
c_inv["IND_BOILER_GAS"]=76.47
c_inv["IND_BOILER_WOOD"]=149.55
c_inv["IND_BOILER_OIL"]=71.23
c_inv["IND_BOILER_COAL"]=149.55
c_inv["IND_BOILER_WASTE"]=149.55
c_inv["IND_DIRECT_ELEC"]=431.53
c_inv["DHN_HP_ELEC"]=447.65
c_inv["DHN_COGEN_GAS"]=1628.89
c_inv["DHN_COGEN_WOOD"]=1403.35
c_inv["DHN_COGEN_WASTE"]=3801.68
c_inv["DHN_BOILER_GAS"]=76.47
c_inv["DHN_BOILER_WOOD"]=149.55
c_inv["DHN_BOILER_OIL"]=71.23
c_inv["DHN_DEEP_GEO"]=2625.8
c_inv["DEC_HP_ELEC"]=638.89
c_inv["DEC_THHP_GAS"]=409.9
c_inv["DEC_COGEN_GAS"]=1828.26
c_inv["DEC_COGEN_OIL"]=1695.2
c_inv["DEC_ADVCOGEN_GAS"]=10799.33
c_inv["DEC_ADVCOGEN_H2"]=10799.33
c_inv["DEC_BOILER_GAS"]=205.85
c_inv["DEC_BOILER_WOOD"]=600.45
c_inv["DEC_BOILER_OIL"]=184.88
c_inv["DEC_SOLAR"]=933.68
c_inv["DEC_DIRECT_ELEC"]=51.91
c_inv["POWER2GAS"]=0.58
c_inv["POWER2GAS_3"]=4354.36
c_inv["EFFICIENCY"]=2585.23
c_inv["DHN"]=1228.48
c_inv["GRID"]=85106.34
c_inv["H2_ELECTROLYSIS"]=458.66
c_inv["H2_NG"]=1015.84
c_inv["H2_BIOMASS"]=3765.73
c_inv["GASIFICATION_SNG"]=4091.43
c_inv["PYROLYSIS"]=2004.43
c_maint["NUCLEAR"]=149.16
c_maint["CCGT"]=28.59
c_maint["CCGT_CCS"]=41.02
c_maint["COAL_US"]=43.04
c_maint["COAL_IGCC"]=70.93
c_maint["COAL_US_CCS"]=91.71
c_maint["COAL_IGCC_CCS"]=100.23
c_maint["PV"]=21.55
c_maint["WIND"]=31.08
c_maint["HYDRO_DAM"]=32.76
c_maint["NEW_HYDRO_DAM"]=3.92
c_maint["HYDRO_RIVER"]=73.1
c_maint["NEW_HYDRO_RIVER"]=103.51
c_maint["GEOTHERMAL"]=631.06
c_maint["IND_COGEN_GAS"]=134.21
c_maint["IND_COGEN_WOOD"]=58.68
c_maint["IND_COGEN_WASTE"]=161.32
c_maint["IND_BOILER_GAS"]=1.71
c_maint["IND_BOILER_WOOD"]=3.34
c_maint["IND_BOILER_OIL"]=1.71
c_maint["IND_BOILER_COAL"]=3.34
c_maint["IND_BOILER_WASTE"]=3.34
c_maint["IND_DIRECT_ELEC"]=2.18
c_maint["DHN_HP_ELEC"]=17.38
c_maint["DHN_COGEN_GAS"]=54.39
c_maint["DHN_COGEN_WOOD"]=58.68
c_maint["DHN_COGEN_WASTE"]=161.32
c_maint["DHN_BOILER_GAS"]=1.71
c_maint["DHN_BOILER_WOOD"]=3.34
c_maint["DHN_BOILER_OIL"]=1.71
c_maint["DHN_DEEP_GEO"]=81.58
c_maint["DEC_HP_ELEC"]=30.51
c_maint["DEC_THHP_GAS"]=13.72
c_maint["DEC_COGEN_GAS"]=134.21
c_maint["DEC_COGEN_OIL"]=118.78
c_maint["DEC_ADVCOGEN_GAS"]=209.9
c_maint["DEC_ADVCOGEN_H2"]=209.9
c_maint["DEC_BOILER_GAS"]=6.89
c_maint["DEC_BOILER_WOOD"]=23.45
c_maint["DEC_BOILER_OIL"]=12.38
c_maint["DEC_SOLAR"]=11.72
c_maint["DEC_DIRECT_ELEC"]=0.26
c_maint["POWER2GAS_3"]=211.58
c_maint["H2_ELECTROLYSIS"]=44.58
c_maint["H2_NG"]=93.38
c_maint["H2_BIOMASS"]=283.6
c_maint["GASIFICATION_SNG"]=202.79
c_maint["PYROLYSIS"]=97.39
lifetime["NUCLEAR"]=44.12
lifetime["CCGT"]=18.382
lifetime["CCGT_CCS"]=18.382
lifetime["COAL_US"]=25.735
lifetime["COAL_IGCC"]=25.735
lifetime["COAL_US_CCS"]=25.735
lifetime["COAL_IGCC_CCS"]=25.735
lifetime["PV"]=18.382
lifetime["WIND"]=14.706
lifetime["HYDRO_DAM"]=29.412
lifetime["NEW_HYDRO_DAM"]=29.412
lifetime["HYDRO_RIVER"]=29.412
lifetime["NEW_HYDRO_RIVER"]=29.412
lifetime["GEOTHERMAL"]=22.059
lifetime["IND_COGEN_GAS"]=18.382
lifetime["IND_COGEN_WOOD"]=18.382
lifetime["IND_COGEN_WASTE"]=18.382
lifetime["IND_BOILER_GAS"]=12.5
lifetime["IND_BOILER_WOOD"]=12.5
lifetime["IND_BOILER_OIL"]=12.5
lifetime["IND_BOILER_COAL"]=12.5
lifetime["IND_BOILER_WASTE"]=12.5
lifetime["IND_DIRECT_ELEC"]=11.029
lifetime["DHN_HP_ELEC"]=18.382
lifetime["DHN_COGEN_GAS"]=18.382
lifetime["DHN_COGEN_WOOD"]=18.382
lifetime["DHN_COGEN_WASTE"]=18.382
lifetime["DHN_BOILER_GAS"]=12.5
lifetime["DHN_BOILER_WOOD"]=12.5
lifetime["DHN_BOILER_OIL"]=12.5
lifetime["DHN_DEEP_GEO"]=22.059
lifetime["DEC_HP_ELEC"]=13.235
lifetime["DEC_THHP_GAS"]=14.706
lifetime["DEC_COGEN_GAS"]=14.706
lifetime["DEC_COGEN_OIL"]=14.706
lifetime["DEC_ADVCOGEN_GAS"]=14.706
lifetime["DEC_ADVCOGEN_H2"]=14.706
lifetime["DEC_BOILER_GAS"]=12.5
lifetime["DEC_BOILER_WOOD"]=12.5
lifetime["DEC_BOILER_OIL"]=12.5
lifetime["DEC_SOLAR"]=14.706
lifetime["DEC_DIRECT_ELEC"]=11.029
lifetime["POWER2GAS_3"]=18.382
lifetime["DHN"]=44.118
lifetime["GRID"]=58.824
lifetime["H2_ELECTROLYSIS"]=11.029
lifetime["H2_NG"]=18.382
lifetime["H2_BIOMASS"]=18.382
lifetime["POWER2GAS"]=18.382
lifetime["GASIFICATION_SNG"]=18.382
lifetime["PYROLYSIS"]=18.382
end_uses_demand_year[("ELECTRICITY","HOUSEHOLDS")]=11317.47
end_uses_demand_year[("ELECTRICITY","SERVICES")]=15637.16
end_uses_demand_year[("ELECTRICITY","INDUSTRY")]=11064.88
end_uses_demand_year[("LIGHTING","HOUSEHOLDS")]=443.49
end_uses_demand_year[("LIGHTING","SERVICES")]=3959.84
end_uses_demand_year[("LIGHTING","INDUSTRY")]=1338.99
end_uses_demand_year[("HEAT_HIGH_T","INDUSTRY")]=20153.26
end_uses_demand_year[("HEAT_LOW_T_SH","HOUSEHOLDS")]=30765.11
end_uses_demand_year[("HEAT_LOW_T_SH","SERVICES")]=15115.08
end_uses_demand_year[("HEAT_LOW_T_SH","INDUSTRY")]=5241.87
end_uses_demand_year[("HEAT_LOW_T_HW","HOUSEHOLDS")]=7863.94
end_uses_demand_year[("HEAT_LOW_T_HW","SERVICES")]=3388.32
end_uses_demand_year[("HEAT_LOW_T_HW","INDUSTRY")]=1358.07
end_uses_demand_year[("MOBILITY_PASSENGER","TRANSPORTATION")]=150972.35
end_uses_demand_year[("MOBILITY_FREIGHT","TRANSPORTATION")]=41313.9
layers_in_out[("NUCLEAR","URANIUM")]=-2.865
layers_in_out[("CCGT","NG")]=-1.6826
layers_in_out[("CCGT_CCS","NG_CCS")]=-1.8597
layers_in_out[("COAL_US","COAL")]=-2.1633
layers_in_out[("COAL_IGCC","COAL")]=-1.9631
layers_in_out[("COAL_US_CCS","COAL_CCS")]=-2.5239
layers_in_out[("COAL_IGCC_CCS","COAL_CCS")]=-2.2084
layers_in_out[("IND_COGEN_GAS","ELECTRICITY")]=0.9023
layers_in_out[("IND_COGEN_WOOD","ELECTRICITY")]=0.3204
layers_in_out[("IND_COGEN_WASTE","ELECTRICITY")]=0.4192
layers_in_out[("IND_BOILER_GAS","NG")]=-1.1432
layers_in_out[("IND_BOILER_WOOD","WOOD")]=-1.2262
layers_in_out[("IND_BOILER_OIL","LFO")]=-1.2149
layers_in_out[("IND_BOILER_COAL","COAL")]=-1.2927
layers_in_out[("IND_BOILER_WASTE","WASTE")]=-1.2927
layers_in_out[("IND_DIRECT_ELEC","ELECTRICITY")]=-1.06
layers_in_out[("DHN_HP_ELEC","ELECTRICITY")]=-0.265
layers_in_out[("DHN_COGEN_GAS","ELECTRICITY")]=1.1792
layers_in_out[("DHN_COGEN_WOOD","ELECTRICITY")]=0.3204
layers_in_out[("DHN_COGEN_WASTE","ELECTRICITY")]=0.353
layers_in_out[("DHN_BOILER_GAS","NG")]=-1.1432
layers_in_out[("DHN_BOILER_WOOD","WOOD")]=-1.2262
layers_in_out[("DHN_BOILER_OIL","LFO")]=-1.2149
layers_in_out[("DEC_HP_ELEC","ELECTRICITY")]=-0.3533
layers_in_out[("DEC_THHP_GAS","NG")]=-0.7067
layers_in_out[("DEC_COGEN_GAS","ELECTRICITY")]=0.9023
layers_in_out[("DEC_COGEN_OIL","ELECTRICITY")]=0.8556
layers_in_out[("DEC_ADVCOGEN_GAS","ELECTRICITY")]=2.0884
layers_in_out[("DEC_ADVCOGEN_H2","ELECTRICITY")]=2.0884
layers_in_out[("DEC_BOILER_GAS","NG")]=-1.1778
layers_in_out[("DEC_BOILER_WOOD","WOOD")]=-1.2471
layers_in_out[("DEC_BOILER_OIL","LFO")]=-1.2471
layers_in_out[("DEC_DIRECT_ELEC","ELECTRICITY")]=-1.06
layers_in_out[("TRAMWAY_TROLLEY","ELECTRICITY")]=-0.1752
layers_in_out[("BUS_COACH_DIESEL","DIESEL")]=-0.2814
layers_in_out[("BUS_COACH_HYDIESEL","DIESEL")]=-0.2308
layers_in_out[("BUS_COACH_CNG_STOICH","NG")]=-0.3246
layers_in_out[("BUS_COACH_FC_HYBRIDH2","H2")]=-0.2847
layers_in_out[("TRAIN_PUB","ELECTRICITY")]=-0.0972
layers_in_out[("CAR_GASOLINE","GASOLINE")]=-0.541
layers_in_out[("CAR_DIESEL","DIESEL")]=-0.487
layers_in_out[("CAR_NG","NG")]=-0.6076
layers_in_out[("CAR_HEV","GASOLINE")]=-0.3465
layers_in_out[("CAR_PHEV","GASOLINE")]=-0.2468
layers_in_out[("CAR_BEV","ELECTRICITY")]=-0.1495
layers_in_out[("CAR_FUEL_CELL","H2")]=-0.2516
layers_in_out[("TRAIN_FREIGHT","ELECTRICITY")]=-0.0724
layers_in_out[("TRUCK","DIESEL")]=-0.5434
layers_in_out[("POWER2GAS_1","LNG")]=0.7471
layers_in_out[("H2_ELECTROLYSIS","ELECTRICITY")]=-1.2471
layers_in_out[("H2_NG","NG")]=-1.4465
layers_in_out[("H2_BIOMASS","WOOD")]=-2.4509
