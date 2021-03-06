using DataFrames, Distributions, Random

Random.seed!(25061985);

if @isdefined nScenarios
    NS = nScenarios;
else
    NS = 1;
end

println("Number of scenarios is: $(NS)")

# leer los parametros inciertos

Parametros_GSA = Dict{Tuple{String,String},Array{Float64,1}}(
("1","i_rate") => [0.0322,0.0173,0.047],
("lifetime","NUCLEAR") => [60,44.12,75.88],
("lifetime","CCGT") => [25,18.382,31.618],
("lifetime","CCGT_CCS") => [25,18.382,31.618],
("lifetime","COAL_US") => [35,25.735,44.265],
("lifetime","COAL_IGCC") => [35,25.735,44.265],
("lifetime","COAL_US_CCS") => [35,25.735,44.265],
("lifetime","COAL_IGCC_CCS") => [35,25.735,44.265],
("lifetime","PV") => [25,18.382,31.618],
("lifetime","WIND") => [20,14.706,25.294],
("lifetime","HYDRO_DAM") => [40,29.412,50.588],
("lifetime","NEW_HYDRO_DAM") => [40,29.412,50.588],
("lifetime","HYDRO_RIVER") => [40,29.412,50.588],
("lifetime","NEW_HYDRO_RIVER") => [40,29.412,50.588],
("lifetime","GEOTHERMAL") => [30,22.059,37.941],
("lifetime","IND_COGEN_GAS") => [25,18.382,31.618],
("lifetime","IND_COGEN_WOOD") => [25,18.382,31.618],
("lifetime","IND_COGEN_WASTE") => [25,18.382,31.618],
("lifetime","IND_BOILER_GAS") => [17,12.5,21.5],
("lifetime","IND_BOILER_WOOD") => [17,12.5,21.5],
("lifetime","IND_BOILER_OIL") => [17,12.5,21.5],
("lifetime","IND_BOILER_COAL") => [17,12.5,21.5],
("lifetime","IND_BOILER_WASTE") => [17,12.5,21.5],
("lifetime","IND_DIRECT_ELEC") => [15,11.029,18.971],
("lifetime","DHN_HP_ELEC") => [25,18.382,31.618],
("lifetime","DHN_COGEN_GAS") => [25,18.382,31.618],
("lifetime","DHN_COGEN_WOOD") => [25,18.382,31.618],
("lifetime","DHN_COGEN_WASTE") => [25,18.382,31.618],
("lifetime","DHN_BOILER_GAS") => [17,12.5,21.5],
("lifetime","DHN_BOILER_WOOD") => [17,12.5,21.5],
("lifetime","DHN_BOILER_OIL") => [17,12.5,21.5],
("lifetime","DHN_DEEP_GEO") => [30,22.059,37.941],
("lifetime","DEC_HP_ELEC") => [18,13.235,22.765],
("lifetime","DEC_THHP_GAS")=> [20,14.706,25.294],
("lifetime","DEC_COGEN_GAS") => [20,14.706,25.294],
("lifetime","DEC_COGEN_OIL") => [20,14.706,25.294],
("lifetime","DEC_ADVCOGEN_GAS") => [20,14.706,25.294],
("lifetime","DEC_ADVCOGEN_H2") => [20,14.706,25.294],
("lifetime","DEC_BOILER_GAS") => [17,12.5,21.5],
("lifetime","DEC_BOILER_WOOD") => [17,12.5,21.5],
("lifetime","DEC_BOILER_OIL") => [17,12.5,21.5],
("lifetime","DEC_SOLAR") => [20,14.706,25.294],
("lifetime","DEC_DIRECT_ELEC") => [15,11.029,18.971],
("lifetime","POWER2GAS_3") => [25,18.382,31.618],
("lifetime","DHN") => [60,44.118,75.882],
("lifetime","GRID") => [80,58.824,101.176],
("lifetime","H2_ELECTROLYSIS") => [15,11.029,18.971],
("lifetime","H2_NG") => [25,18.382,31.618],
("lifetime","H2_BIOMASS") => [25,18.382,31.618],
("lifetime","POWER2GAS") => [25,18.382,31.618],
("lifetime","GASIFICATION_SNG") => [25,18.382,31.618],
("lifetime","PYROLYSIS") => [25,18.382,31.618],
("c_inv","NUCLEAR") => [ 5174.76,4057.59,11346.7],
("c_inv","CCGT") => [824.41,646.43,1030.49],
("c_inv","CCGT_CCS") => [1273.26,768.62,1777.9],				# new
("c_inv","COAL_US") => [2687.59,2107.37,3359.4],
("c_inv","COAL_IGCC") => [3466.2,2092.41,4839.99],				# new
("c_inv","COAL_US_CCS") => [4327.26,2612.2,6042.32],			# new
("c_inv","COAL_IGCC_CCS") => [6044.78,3649,8440.56],			# new
("c_inv","PV") => [1000,603.66,1396.34],
("c_inv","WIND") => [1465.62,1149.21,1801.55],
("c_inv","HYDRO_DAM") => [4828.39,3786,8393.05], 				# comentar out-of-sample
("c_inv","NEW_HYDRO_DAM") => [3437.12,2695.09,5974.64], 		# comentar out-of-sample
("c_inv","HYDRO_RIVER") => [5387.47,4224.38,6550.56], 			# mature
("c_inv","NEW_HYDRO_RIVER") => [5919.19,4641.31,7197.07],		# mature
("c_inv","GEOTHERMAL") => [11464.13,6910.49,18580.33],
("c_inv","IND_COGEN_GAS") => [1503.64,1179.02,1828.26],			# mature
("c_inv","IND_COGEN_WOOD") => [1154.18,905.01,1403.35],			# mature
("c_inv","IND_COGEN_WASTE") => [3126.67,2451.66,3801.68],		# mature
("c_inv","IND_BOILER_GAS") => [62.89,49.31,76.47],				# mature
("c_inv","IND_BOILER_WOOD") => [123,96.45,149.55],				# mature
("c_inv","IND_BOILER_OIL") => [58.58,45.93,71.23],				# mature
("c_inv","IND_BOILER_COAL") => [123,96.45,149.55],				# mature
("c_inv","IND_BOILER_WASTE") => [123,96.45,149.55],				# mature
("c_inv","IND_DIRECT_ELEC") => [354.91,278.29,431.53],			# mature
("c_inv","DHN_HP_ELEC") => [368.17,288.69,447.65],				# mature
("c_inv","DHN_COGEN_GAS") => [1339.67,1050.45,1628.89],			# mature
("c_inv","DHN_COGEN_WOOD") => [1154.18,905.01,1403.35],			# mature
("c_inv","DHN_COGEN_WASTE") => [3126.67,2451.66,3801.68],		# mature
("c_inv","DHN_BOILER_GAS") => [62.89,49.31,76.47],				# mature
("c_inv","DHN_BOILER_WOOD") => [123,96.45,149.55],				# mature
("c_inv","DHN_BOILER_OIL") => [58.58,45.93,71.23],				# mature
("c_inv","DHN_DEEP_GEO") => [1620.13,976.6,2625.8],
("c_inv","DEC_HP_ELEC") => [525.45,412.01,638.89],				# mature
("c_inv","DEC_THHP_GAS") => [337.12,264.34,409.9],				# mature
("c_inv","DEC_COGEN_GAS") => [1503.64,1179.02,1828.26],			# mature
("c_inv","DEC_COGEN_OIL") => [1394.21,1093.22,1695.2],			# mature
("c_inv","DEC_ADVCOGEN_GAS") => [7734.03,4668.73,10799.33],		# new
("c_inv","DEC_ADVCOGEN_H2") => [7734.03,4668.73,10799.33],		# new
("c_inv","DEC_BOILER_GAS") => [169.3,132.75,205.85],			# mature
("c_inv","DEC_BOILER_WOOD") => [493.84,387.23,600.45],			# mature
("c_inv","DEC_BOILER_OIL") => [152.05,119.22,184.88],			# mature
("c_inv","DEC_SOLAR") => [767.9,602.12,933.68],					# mature
("c_inv","DEC_DIRECT_ELEC") => [42.69,33.47,51.91],				# mature
("c_inv","POWER2GAS") => [0.4144,0.25,0.578643055],				# new
("c_inv","POWER2GAS_3") => [3118.41,1882.46,4354.36],			# new
("c_inv","EFFICIENCY") => [1856,1126.77,2585.23], 				# comentar out-of-sample
("c_inv","DHN") => [881.96,535.44,1228.48],
("c_inv","GRID") => [61100,37093.66,85106.34], 					# comentar out-of-sample
("c_inv","H2_ELECTROLYSIS") => [328.47,198.28,458.66],			# new
("c_inv","H2_NG") => [727.5,439.16,1015.84],					# new
("c_inv","H2_BIOMASS") => [2696.86,1627.99,3765.73],			# new
("c_inv","GASIFICATION_SNG") => [2930.11,1768.79,4091.43],
("c_inv","PYROLYSIS") => [1435.49,866.55,2004.43],
#("c_maint","perc") => [0.0295,0.01528,0.04],
("c_maint","NUCLEAR")=>[109.92,56.9803963257,149.1850376526],
("c_maint","CCGT")=>[21.07,10.9222793903,28.5965133128],
("c_maint","CCGT_CCS")=>[30.23,15.6706457508,41.0285997838],
("c_maint","COAL_US")=>[31.72,16.4430328553,43.0508496574],
("c_maint","COAL_IGCC")=>[52.27,27.0957543299,70.9416113365],
("c_maint","COAL_US_CCS")=>[67.58,35.0321614237,91.720568091],
("c_maint","COAL_IGCC_CCS")=>[73.86,38.2875916359,100.243876283],
("c_maint","PV")=>[15.88,8.2318840398,21.5525691223],
("c_maint","WIND")=>[22.9,11.8709159012,31.0802161776],
("c_maint","HYDRO_DAM")=>[24.14,12.5137078539,32.7631623811],
("c_maint","NEW_HYDRO_DAM")=>[2.89,1.4981199543,3.9223504259],
("c_maint","HYDRO_RIVER")=>[53.87,27.9251633012,73.1131548248],
("c_maint","NEW_HYDRO_RIVER")=>[76.28,39.5420727049,103.5283358092],
("c_maint","GEOTHERMAL")=>[465.04,241.0677174971,631.1591149015],
("c_maint","IND_COGEN_GAS")=>[98.9,51.2678420361,134.2285318763],
("c_maint","IND_COGEN_WOOD")=>[43.24,22.4147774483,58.6859627738],
("c_maint","IND_COGEN_WASTE")=>[118.88,61.6250865647,161.3456811876],
("c_maint","IND_BOILER_GAS")=>[1.26,0.6531595649,1.7100904971],
("c_maint","IND_BOILER_WOOD")=>[2.46,1.2752162933,3.3387481134],
("c_maint","IND_BOILER_OIL")=>[1.26,0.6531595649,1.7100904971],
("c_maint","IND_BOILER_COAL")=>[2.46,1.2752162933,3.3387481134],
("c_maint","IND_BOILER_WASTE")=>[2.46,1.2752162933,3.3387481134],
("c_maint","IND_DIRECT_ELEC")=>[1.61,0.8345927773,2.1851156352],
("c_maint","DHN_HP_ELEC")=>[12.81,6.6404555762,17.3859200539],
("c_maint","DHN_COGEN_GAS")=>[40.08,20.7766947301,54.3971643843],
("c_maint","DHN_COGEN_WOOD")=>[43.24,22.4147774483,58.6859627738],
("c_maint","DHN_COGEN_WASTE")=>[118.88,61.6250865647,161.3456811876],
("c_maint","DHN_BOILER_GAS")=>[1.26,0.6531595649,1.7100904971],
("c_maint","DHN_BOILER_WOOD")=>[2.46,1.2752162933,3.3387481134],
("c_maint","DHN_BOILER_OIL")=>[1.26,0.6531595649,1.7100904971],
("c_maint","DHN_DEEP_GEO")=>[60.12,31.1650420951,81.5957465764],
("c_maint","DEC_HP_ELEC")=>[22.48,11.6531960462,30.5101860119],
("c_maint","DEC_THHP_GAS")=>[10.11,5.2408279372,13.7214404173],
("c_maint","DEC_COGEN_GAS")=>[98.9,51.2678420361,134.2285318763],
("c_maint","DEC_COGEN_OIL")=>[87.53,45.3738545341,118.7970009619],
("c_maint","DEC_ADVCOGEN_GAS")=>[154.68,80.1831122967,209.9339667404],
("c_maint","DEC_ADVCOGEN_H2")=>[154.68,80.1831122967,209.9339667404],
("c_maint","DEC_BOILER_GAS")=>[5.08,2.6333734838,6.8946505756],
("c_maint","DEC_BOILER_WOOD")=>[17.28,8.9576168896,23.4526696746],
("c_maint","DEC_BOILER_OIL")=>[9.12,4.7276311362,12.3777978838],
("c_maint","DEC_SOLAR")=>[8.64,4.4788084448,11.7263348373],
("c_maint","DEC_DIRECT_ELEC")=>[0.19,0.0984923153,0.2578707892],
("c_maint","POWER2GAS_3")=>[155.92,80.8259042494,211.6169129439],
("c_maint","H2_ELECTROLYSIS")=>[32.85,17.0288029412,44.5845022461],
("c_maint","H2_NG")=>[68.81,35.6697695703,93.3899421477],
("c_maint","H2_BIOMASS")=>[208.99,108.3363630649,283.6442960246],
("c_maint","GASIFICATION_SNG")=>[149.44,77.4667979158,202.8221618159],
("c_maint","PYROLYSIS")=>[71.77,37.2041761671,97.4072976012],
("c_op","ELECTRICITY") => [0.0901,0.0901*(1-0.473),0.0901*(1+0.8992)],
("c_op","GASOLINE") => [0.0880,0.0880*(1-0.473),0.0880*(1+0.8992)],
("c_op","DIESEL") => [0.0852,0.0852*(1-0.473),0.0852*(1+0.8992)],
("c_op","LFO") => [0.0606,0.0606*(1-0.473),0.0606*(1+0.8992)],
("c_op","NG") => [0.0348,0.0348*(1-0.473),0.0348*(1+0.8992)],
("c_op","WOOD") => [0.0932,0.0905,0.0959],
("c_op","COAL") => [0.0302,0.0302*(1-0.473),0.0302*(1+0.8992)],
("c_op","URANIUM") => [0.0041,0.0041*(1-0.473),0.0041*(1+0.8992)],
("c_p_t","PV") => [1,0.8893,1.1107],
("c_p_t","WIND") => [1,0.8893,1.1107],
("c_p_t","HYDRO_DAM") => [1,0.8893,1.1107],
("c_p_t","HYDRO_RIVER") => [1,0.8893,1.1107],
("c_p_t","NEW_HYDRO_DAM") => [1,0.8893,1.1107],
("c_p_t","NEW_HYDRO_RIVER") => [1,0.8893,1.1107],
("c_p_t","DEC_SOLAR") => [1,0.8893,1.1107],
("1","peak_dhn_factor") => [2,1.6,2.4],
("f_max","PV") => [25,18.983,31.017],
("f_max","WIND") => [5.3,4.024,6.576],
("f_max","NEW_HYDRO_DAM") => [0.4,0.334,0.546],
("f_max","NEW_HYDRO_RIVER") => [0.9,0.645,1.055],
("f_max","GEOTHERMAL") => [0.7,0.532,0.868],
("fmax_perc","PV")=>[0.5,0.4,0.6],
("fmax_perc","WIND")=>[0.5,0.4,0.6],
("fmax_perc","IND_COGEN_GAS")=>[0.5,0.4,0.6],
("fmax_perc","IND_COGEN_WASTE")=>[0.5,0.4,0.6],
("fmax_perc","IND_BOILER_GAS")=>[0.6,0.48,0.72],
("fmax_perc","IND_BOILER_OIL")=>[0.5,0.4,0.6],
("fmax_perc","IND_BOILER_COAL")=>[0.5,0.4,0.6],
("fmax_perc","IND_DIRECT_ELEC")=>[0.2,0.16,0.24],
("fmax_perc","DHN_HP_ELEC")=>[0.5,0.4,0.6],
("fmax_perc","DHN_COGEN_GAS")=>[0.5,0.4,0.6],
("fmax_perc","DHN_COGEN_WASTE")=>[0.5,0.4,0.6],
("fmax_perc","DHN_BOILER_GAS")=>[0.8,0.64,0.96],
("fmax_perc","DHN_BOILER_OIL")=>[0.5,0.4,0.6],
("fmax_perc","DHN_DEEP_GEO")=>[0.5,0.4,0.6],
("fmax_perc","DEC_HP_ELEC")=>[0.5,0.4,0.6],
("fmax_perc","DEC_THHP_GAS")=>[0.2,0.16,0.24],
("fmax_perc","DEC_COGEN_GAS")=>[0.4,0.32,0.48],
("fmax_perc","DEC_COGEN_OIL")=>[0.4,0.32,0.48],
("fmax_perc","DEC_ADVCOGEN_GAS")=>[0.2,0.16,0.24],
("fmax_perc","DEC_ADVCOGEN_H2")=>[0.2,0.16,0.24],
("fmax_perc","DEC_BOILER_GAS")=>[0.8,0.64,0.96],
("fmax_perc","DEC_BOILER_OIL")=>[0.5,0.4,0.6],
("fmax_perc","DEC_SOLAR")=>[0.4,0.32,0.48],
("fmax_perc","DEC_DIRECT_ELEC")=>[0.2,0.16,0.24],
("fmax_perc","TRAMWAY_TROLLEY")=>[0.3,0.24,0.36],
("fmax_perc","BUS_COACH_DIESEL")=>[0.3,0.24,0.36],
("fmax_perc","BUS_COACH_HYDIESEL")=>[0.3,0.24,0.36],
("fmax_perc","BUS_COACH_CNG_STOICH")=>[0.3,0.24,0.36],
("fmax_perc","BUS_COACH_FC_HYBRIDH2")=>[0.2,0.16,0.24],
("fmax_perc","TRAIN_PUB")=>[0.8,0.64,0.96],
("fmax_perc","CAR_NG")=>[0.5,0.4,0.6],
("fmax_perc","CAR_HEV")=>[0.3,0.24,0.36],
("fmax_perc","CAR_PHEV")=>[0.3,0.24,0.36],
("fmax_perc","CAR_BEV")=>[0.3,0.24,0.36],
("fmax_perc","CAR_FUEL_CELL")=>[0.2,0.16,0.24],
("fmin_perc","DHN_BOILER_GAS")=>[0.2,0.16,0.24],
("fmin_perc","DEC_BOILER_GAS")=>[0.2,0.16,0.24],
("fmin_perc","DEC_BOILER_OIL")=>[0.1,0.08,0.12],
("fmin_perc","CAR_GASOLINE")=>[0.2,0.16,0.24],
("fmin_perc","CAR_DIESEL")=>[0.2,0.16,0.24],
("c_p","NUCLEAR")=>[0.85,0.8282,0.8693],
("c_p","CCGT")=>[0.85,0.8292,0.8703],
("c_p","CCGT_CCS")=>[0.85,0.8292,0.8703],
("c_p","COAL_US")=>[0.87,0.8467,0.8887],
("c_p","COAL_IGCC")=>[0.86,0.835,0.8764],
("c_p","COAL_US_CCS")=>[0.87,0.8467,0.8887],
("c_p","COAL_IGCC_CCS")=>[0.86,0.835,0.8764],
("c_p","GEOTHERMAL")=>[0.86,0.8389,0.8805],
("c_p","IND_COGEN_GAS")=>[0.85,0.8292,0.8703],
("c_p","IND_COGEN_WOOD")=>[0.85,0.8292,0.8703],
("c_p","IND_COGEN_WASTE")=>[0.85,0.8292,0.8703],
("c_p","IND_BOILER_GAS")=>[0.95,0.9267,0.9727],
("c_p","IND_BOILER_WOOD")=>[0.9,0.878,0.9215],
("c_p","IND_BOILER_OIL")=>[0.95,0.9267,0.9727],
("c_p","IND_BOILER_COAL")=>[0.9,0.878,0.9215],
("c_p","IND_BOILER_WASTE")=>[0.9,0.878,0.9215],
("c_p","IND_DIRECT_ELEC")=>[0.95,0.9267,0.9727],
("c_p","DHN_HP_ELEC")=>[0.95,0.9267,0.9727],
("c_p","DHN_COGEN_GAS")=>[0.85,0.8292,0.8703],
("c_p","DHN_COGEN_WOOD")=>[0.85,0.8292,0.8703],
("c_p","DHN_COGEN_WASTE")=>[0.85,0.8292,0.8703],
("c_p","DHN_BOILER_GAS")=>[0.95,0.9267,0.9727],
("c_p","DHN_BOILER_WOOD")=>[0.9,0.878,0.9215],
("c_p","DHN_BOILER_OIL")=>[0.95,0.9267,0.9727],
("c_p","DHN_DEEP_GEO")=>[0.85,0.8292,0.8703],
("c_p","POWER2GAS_1")=>[0.9,0.878,0.9215],
("c_p","POWER2GAS_2")=>[0.9,0.878,0.9215],
("c_p","H2_ELECTROLYSIS")=>[0.9,0.878,0.9215],
("c_p","H2_NG")=>[0.86,0.8389,0.8805],
("c_p","H2_BIOMASS")=>[0.86,0.84,0.88],
("c_p","GASIFICATION_SNG")=>[0.85,0.83,0.87],
("c_p","PYROLYSIS")=>[0.85,0.83,0.87],
("ref_size","NUCLEAR")=>[1,0.8,1.2],
("ref_size","CCGT")=>[0.5,0.4,0.6],
("ref_size","CCGT_CCS")=>[0.5,0.4,0.6],
("ref_size","COAL_US")=>[0.5,0.4,0.6],
("ref_size","COAL_IGCC")=>[0.5,0.4,0.6],
("ref_size","COAL_US_CCS")=>[0.5,0.4,0.6],
("ref_size","COAL_IGCC_CCS")=>[0.5,0.4,0.6],
("ref_size","PV")=>[0.000003,0.0000024,0.0000036],
("ref_size","WIND")=>[0.003,0.0024,0.0036],
("ref_size","NEW_HYDRO_DAM")=>[0.001,0.0008,0.0012],
("ref_size","NEW_HYDRO_RIVER")=>[0.0001,0.00008,0.00012],
("ref_size","GEOTHERMAL")=>[0.008,0.0064,0.0096],
("ref_size","IND_COGEN_GAS")=>[0.02,0.016,0.024],
("ref_size","IND_COGEN_WOOD")=>[0.02,0.016,0.024],
("ref_size","IND_COGEN_WASTE")=>[0.02,0.016,0.024],
("ref_size","IND_BOILER_GAS")=>[0.01,0.008,0.012],
("ref_size","IND_BOILER_WOOD")=>[0.01,0.008,0.012],
("ref_size","IND_BOILER_OIL")=>[0.01,0.008,0.012],
("ref_size","IND_BOILER_COAL")=>[0.001,0.0008,0.0012],
("ref_size","IND_BOILER_WASTE")=>[0.001,0.0008,0.0012],
("ref_size","IND_DIRECT_ELEC")=>[0.0001,0.00008,0.00012],
("ref_size","DHN_HP_ELEC")=>[0.001,0.0008,0.0012],
("ref_size","DHN_COGEN_GAS")=>[0.02,0.016,0.024],
("ref_size","DHN_COGEN_WOOD")=>[0.02,0.016,0.024],
("ref_size","DHN_COGEN_WASTE")=>[0.02,0.016,0.024],
("ref_size","DHN_BOILER_GAS")=>[0.01,0.008,0.012],
("ref_size","DHN_BOILER_WOOD")=>[0.01,0.008,0.012],
("ref_size","DHN_BOILER_OIL")=>[0.01,0.008,0.012],
("ref_size","DHN_DEEP_GEO")=>[0.023,0.0184,0.0276],
("ref_size","DEC_HP_ELEC")=>[0.00001,0.000008,0.000012],
("ref_size","DEC_THHP_GAS")=>[0.00001,0.000008,0.000012],
("ref_size","DEC_COGEN_GAS")=>[0.000005,0.000004,0.000006],
("ref_size","DEC_COGEN_OIL")=>[0.00001,0.000008,0.000012],
("ref_size","DEC_ADVCOGEN_GAS")=>[0.00001,0.000008,0.000012],
("ref_size","DEC_ADVCOGEN_H2")=>[0.00001,0.000008,0.000012],
("ref_size","DEC_BOILER_GAS")=>[0.00001,0.000008,0.000012],
("ref_size","DEC_BOILER_WOOD")=>[0.00001,0.000008,0.000012],
("ref_size","DEC_BOILER_OIL")=>[0.00001,0.000008,0.000012],
("ref_size","DEC_SOLAR")=>[0.00001,0.000008,0.000012],
("ref_size","DEC_DIRECT_ELEC")=>[0.00001,0.000008,0.000012],
("ref_size","TRAMWAY_TROLLEY")=>[0.001,0.0008,0.0012],
("ref_size","BUS_COACH_DIESEL")=>[0.001,0.0008,0.0012],
("ref_size","BUS_COACH_HYDIESEL")=>[0.001,0.0008,0.0012],
("ref_size","BUS_COACH_CNG_STOICH")=>[0.001,0.0008,0.0012],
("ref_size","BUS_COACH_FC_HYBRIDH2")=>[0.001,0.0008,0.0012],
("ref_size","TRAIN_PUB")=>[0.001,0.0008,0.0012],
("ref_size","CAR_GASOLINE")=>[0.001,0.0008,0.0012],
("ref_size","CAR_DIESEL")=>[0.001,0.0008,0.0012],
("ref_size","CAR_NG")=>[0.001,0.0008,0.0012],
("ref_size","CAR_HEV")=>[0.001,0.0008,0.0012],
("ref_size","CAR_PHEV")=>[0.001,0.0008,0.0012],
("ref_size","CAR_BEV")=>[0.001,0.0008,0.0012],
("ref_size","CAR_FUEL_CELL")=>[0.001,0.0008,0.0012],
("ref_size","TRAIN_FREIGHT")=>[0.001,0.0008,0.0012],
("ref_size","TRUCK")=>[0.001,0.0008,0.0012],
("ref_size","PUMPED_HYDRO")=>[0.001,0.0008,0.0012],
("ref_size","POWER2GAS")=>[0.001,0.0008,0.0012],
("ref_size","POWER2GAS_1")=>[0.001,0.0008,0.0012],
("ref_size","POWER2GAS_2")=>[0.001,0.0008,0.0012],
("ref_size","POWER2GAS_3")=>[0.001,0.0008,0.0012],
("avail","WOOD")=>[12279,8340.45,16217.55],
("avail","WASTE") => [11142,7568.15,14715.85],
("1","share_mobility_public_min")=>[0.3,0.24,0.36],
("1","share_mobility_public_max")=>[0.5,0.4,0.6],
("1","share_freight_train_min")=>[0.4,0.32,0.479],
("1","share_freight_train_max")=>[0.6,0.48,0.72],
("1","share_heat_dhn_min")=>[0.1,0.08,0.12],
("1","share_heat_dhn_max")=>[0.3,0.24,0.36],
("loss_coeff","ELECTRICITY")=>[0.07,0.0686,0.0714],
("loss_coeff","HEAT_LOW_T_DHN")=>[0.05,0.049,0.051],
("lighting_month","1")=>[0.124,0.0992,0.1488],
("lighting_month","2")=>[0.102,0.0816,0.1224],
("lighting_month","3")=>[0.0918,0.0734,0.1102],
("lighting_month","4")=>[0.0708,0.0566,0.085],
("lighting_month","5")=>[0.0588,0.047,0.0706],
("lighting_month","6")=>[0.0438,0.035,0.0526],
("lighting_month","7")=>[0.0408,0.0326,0.049],
("lighting_month","8")=>[0.0398,0.0318,0.0478],
("lighting_month","9")=>[0.0678,0.0542,0.0814],
("lighting_month","10")=>[0.1048,0.0838,0.1258],
("lighting_month","11")=>[0.1238,0.099,0.1486],
("lighting_month","12")=>[0.1318,0.1054,0.1582],
("heating_month","1")=>[0.1984,0.1587,0.2381],
("heating_month","2")=>[0.1654,0.1323,0.1985],
("heating_month","3")=>[0.1421,0.1137,0.1705],
("heating_month","4")=>[0.0319,0.0255,0.0383],
("heating_month","5")=>[0,0,0],
("heating_month","6")=>[0,0,0],
("heating_month","7")=>[0,0,0],
("heating_month","8")=>[0,0,0],
("heating_month","9")=>[0.0147,0.0118,0.0176],
("heating_month","10")=>[0.0898,0.0718,0.1078],
("heating_month","11")=>[0.1383,0.1106,0.166],
("heating_month","12")=>[0.2194,0.1755,0.2633]
)


Parametros_GSA_big = Dict{Tuple{String,String,String},Array{Float64,1}}(
("end_uses_demand_year","ELECTRICITY","HOUSEHOLDS")=>[10848.1,10095.5694575898,11317.4668168738],
("end_uses_demand_year","ELECTRICITY","SERVICES")=>[15026.5,13916.6181665504,15637.1646491666],
("end_uses_demand_year","ELECTRICITY","INDUSTRY")=>[10443.5,9343.7705390651,11064.8764607173],
("end_uses_demand_year","LIGHTING","HOUSEHOLDS")=>[425.1,395.610897431,443.4928829798],
("end_uses_demand_year","LIGHTING","SERVICES")=>[3805.2,3524.1417127979,3959.8402104954],
("end_uses_demand_year","LIGHTING","INDUSTRY")=>[1263.8,1130.718361399,1338.9946733427],
("end_uses_demand_year","HEAT_HIGH_T","INDUSTRY")=>[19021.5,17018.4833924285,20153.2577773288],
("end_uses_demand_year","HEAT_LOW_T_SH","HOUSEHOLDS")=>[29489.2,27443.5400529823,30765.1148547813],
("end_uses_demand_year","HEAT_LOW_T_SH","SERVICES")=>[14524.8,13451.9745479992,15115.0759721968],
("end_uses_demand_year","HEAT_LOW_T_SH","INDUSTRY")=>[4947.5,4426.5145537439,5241.8706649494],
("end_uses_demand_year","HEAT_LOW_T_HW","HOUSEHOLDS")=>[7537.8,7014.9043111163,7863.9394338392],
("end_uses_demand_year","HEAT_LOW_T_HW","SERVICES")=>[3256,3015.5065218306,3388.3211724411],
("end_uses_demand_year","HEAT_LOW_T_HW","INDUSTRY")=>[1281.8,1146.8229115693,1358.0656530232],
("end_uses_demand_year","MOBILITY_PASSENGER","TRANSPORTATION")=>[146049.3,141126.251037864,150972.348962136],
("end_uses_demand_year","MOBILITY_FREIGHT","TRANSPORTATION")=>[39966.7,38619.4972338449,41313.9027661551],
("layers_in_out","NUCLEAR","URANIUM")=>[-2.7027,-2.8649,-2.5578],
("layers_in_out","CCGT","NG")=>[-1.5873,-1.6825,-1.5022],
("layers_in_out","CCGT_CCS","NG_CCS")=>[-1.7544,-1.8597,-1.6603],
("layers_in_out","COAL_US","COAL")=>[-2.0408,-2.1633,-1.9314],
("layers_in_out","COAL_IGCC","COAL")=>[-1.8519,-1.9630,-1.7526],
("layers_in_out","COAL_US_CCS","COAL_CCS")=>[-2.381,-2.5239,-2.2533],
("layers_in_out","COAL_IGCC_CCS","COAL_CCS")=>[-2.0833,-2.2083,-1.9716],
("layers_in_out","IND_COGEN_GAS","ELECTRICITY")=>[0.9565,0.9023,1.0107],
#("layers_in_out","IND_COGEN_GAS","NG")=>[-2.1739,-2.3044,-2.0574],
("layers_in_out","IND_COGEN_WOOD","ELECTRICITY")=>[0.3396,0.3203,0.3588],
#("layers_in_out","IND_COGEN_WOOD","WOOD")=>[-1.8868,-2.0001,-1.7857],
("layers_in_out","IND_COGEN_WASTE","ELECTRICITY")=>[0.4444,0.4192,0.4695],
#("layers_in_out","IND_COGEN_WASTE","WASTE")=>[-2.2222,-2.3556,-2.1031],
("layers_in_out","IND_BOILER_GAS","NG")=>[-1.0785,-1.1432,-1.0206],
("layers_in_out","IND_BOILER_WOOD","WOOD")=>[-1.1568,-1.2262,-1.0947],
("layers_in_out","IND_BOILER_OIL","LFO")=>[-1.1461,-1.2149,-1.0846],
("layers_in_out","IND_BOILER_COAL","COAL")=>[-1.2195,-1.2927,-1.1541],
("layers_in_out","IND_BOILER_WASTE","WASTE")=>[-1.2195,-1.2927,-1.1541],
("layers_in_out","IND_DIRECT_ELEC","ELECTRICITY")=>[-1,-1.06,-1], # direct electric heating pag 114
("layers_in_out","DHN_HP_ELEC","ELECTRICITY")=>[-0.25,-0.265,-0.2366],
("layers_in_out","DHN_COGEN_GAS","ELECTRICITY")=>[1.25,1.1792,1.3207],
#("layers_in_out","DHN_COGEN_GAS","NG")=>[-2.5,-2.6501,-2.366],
("layers_in_out","DHN_COGEN_WOOD","ELECTRICITY")=>[0.3396,0.3203,0.3588],
#("layers_in_out","DHN_COGEN_WOOD","WOOD")=>[-1.8868,-2.0001,-1.7857],
("layers_in_out","DHN_COGEN_WASTE","ELECTRICITY")=>[0.4444,0.3529,0.5358],
#("layers_in_out","DHN_COGEN_WASTE","WASTE")=>[-2.2222,-2.7979,-1.8429],
("layers_in_out","DHN_BOILER_GAS","NG")=>[-1.0785,-1.1432,-1.0206],
("layers_in_out","DHN_BOILER_WOOD","WOOD")=>[-1.1568,-1.2262,-1.0947],
("layers_in_out","DHN_BOILER_OIL","LFO")=>[-1.1461,-1.2149,-1.0846],
("layers_in_out","DEC_HP_ELEC","ELECTRICITY")=>[-0.3333,-0.3533,-0.3154],
("layers_in_out","DEC_THHP_GAS","NG")=>[-0.6667,-0.7067,-0.6309],
("layers_in_out","DEC_COGEN_GAS","ELECTRICITY")=>[0.9565,0.9023,1.0106],
#("layers_in_out","DEC_COGEN_GAS","NG")=>[-2.1739,-2.3044,-2.0573],
("layers_in_out","DEC_COGEN_OIL","ELECTRICITY")=>[0.907,0.8556,0.9583],
#("layers_in_out","DEC_COGEN_OIL","LFO")=>[-2.3256,-2.4652,-2.2009],
("layers_in_out","DEC_ADVCOGEN_GAS","ELECTRICITY")=>[2.6364,2.0883,3.1844],
#("layers_in_out","DEC_ADVCOGEN_GAS","NG")=>[-4.5455,-5.7382,-3.7632],
("layers_in_out","DEC_ADVCOGEN_H2","ELECTRICITY")=>[2.6364,2.0883,3.1844],
#("layers_in_out","DEC_ADVCOGEN_H2","H2")=>[-4.5455,-5.7382,-3.7632],
("layers_in_out","DEC_BOILER_GAS","NG")=>[-1.1111,-1.1778,-1.0515],
("layers_in_out","DEC_BOILER_WOOD","WOOD")=>[-1.1765,-1.2471,-1.1134],
("layers_in_out","DEC_BOILER_OIL","LFO")=>[-1.1765,-1.2471,-1.1134],
("layers_in_out","DEC_DIRECT_ELEC","ELECTRICITY")=>[-1,-1.06,-1],
("layers_in_out","TRAMWAY_TROLLEY","ELECTRICITY")=>[-0.1653,-0.1752,-0.1564],
("layers_in_out","BUS_COACH_DIESEL","DIESEL")=>[-0.2655,-0.2814,-0.2512],
("layers_in_out","BUS_COACH_HYDIESEL","DIESEL")=>[-0.1828,-0.2307,-0.1513],
("layers_in_out","BUS_COACH_CNG_STOICH","NG")=>[-0.3062,-0.3245,-0.2897],
("layers_in_out","BUS_COACH_FC_HYBRIDH2","H2")=>[-0.2255,-0.2846,-0.1866],
("layers_in_out","TRAIN_PUB","ELECTRICITY")=>[-0.0917,-0.0972,-0.0867],
("layers_in_out","CAR_GASOLINE","GASOLINE")=>[-0.4297,-0.541,-0.3563],
("layers_in_out","CAR_DIESEL","DIESEL")=>[-0.3868,-0.487,-0.3207],
("layers_in_out","CAR_NG","NG")=>[-0.4826,-0.6076,-0.4002],
("layers_in_out","CAR_HEV","GASOLINE")=>[-0.2471,-0.3465,-0.192],
#("layers_in_out","CAR_PHEV","ELECTRICITY")=>[-0.0451,-0.0632,-0.035],
("layers_in_out","CAR_PHEV","GASOLINE")=>[-0.176,-0.2468,-0.1367],
("layers_in_out","CAR_BEV","ELECTRICITY")=>[-0.1066,-0.1495,-0.0828],
("layers_in_out","CAR_FUEL_CELL","H2")=>[-0.1794,-0.2515,-0.1393],
("layers_in_out","TRAIN_FREIGHT","ELECTRICITY")=>[-0.0683,-0.0724,-0.0646],
("layers_in_out","TRUCK","DIESEL")=>[-0.5126,-0.5433,-0.4851],
("layers_in_out","POWER2GAS_1","LNG")=>[0.792,0.7471,0.8368],
("layers_in_out","H2_ELECTROLYSIS","ELECTRICITY")=>[-1.1765,-1.2471,-1.1134],
("layers_in_out","H2_NG","NG")=>[-1.3646,-1.4465,-1.2914],
("layers_in_out","H2_BIOMASS","WOOD")=>[-2.3121,-2.4509,-2.1881],
#("layers_in_out","GASIFICATION_SNG","WOOD")=>[-1.3514,-1.4325,-1.279],
#("layers_in_out","PYROLYSIS","WOOD")=>[-1.5017,-1.5919,-1.4212]
)


# Estapa I
        
##################### i_rate

LISTA_irate= ["i_rate"]
LISTA_irate = []

println("Parameter i_rate $(length(LISTA_irate))")

M_i_rate = Dict{String,Array{Float64}}();

for (k,v) in Dict(Parametros_GSA)
    if in(k[2],LISTA_irate) 
        d1 = Uniform(Parametros_GSA[k][2],Parametros_GSA[k][3])
        M_i_rate[k[2]] = [rand(d1,1)[1] for i=1:NS]
    end
end

try
   M_i_rate["i_rate"]
catch error
   if isa(error, KeyError)
       M_i_rate["i_rate"] = [i_rate for j=1:NS]; # default value
   end
end


##################### lifetime

LISTA_lifetime = [k[2] for (k,v) in Parametros_GSA if k[1] == "lifetime"]
LISTA_lifetime = []

println("Parameter lifetime $(length(LISTA_lifetime))")

M_lifetime = Dict{String,Array{Float64}}();

for (k,v) in Dict(Parametros_GSA)
    if k[1] == "lifetime" && in(k[2],LISTA_lifetime)
        d1 = Uniform(Parametros_GSA[k][2],Parametros_GSA[k][3])
        M_lifetime[k[2]] = [rand(d1,1)[1] for i=1:NS]
    end
end

for i in TECHNOLOGIES
    try
    M_lifetime[i]
    catch error
       if isa(error, KeyError)
           M_lifetime[i] = [lifetime[i] for j=1:NS]; # default value
       end
    end
end

##################### c_inv 

LISTA_c_inv = [k[2] for (k,v) in Parametros_GSA if k[1] == "c_inv"]
LISTA_c_inv = []

println("Parameter c_inv $(length(LISTA_c_inv))")

M_c_inv = Dict{String,Array{Float64}}();

for (k,v) in Dict(Parametros_GSA)
    if k[1] == "c_inv" && in(k[2],LISTA_c_inv)
        d1 = Uniform(Parametros_GSA[k][2],Parametros_GSA[k][3])
        M_c_inv[k[2]] = [rand(d1,1)[1] for i=1:NS]
    end
end

for i in TECHNOLOGIES
    try
    M_c_inv[i]
    catch error
       if isa(error, KeyError)
           M_c_inv[i] = [c_inv[i] for j=1:NS]; # default value
       end
    end
end


##################### c_maint

LISTA_c_maint = [k[2] for (k,v) in Parametros_GSA if k[1] == "c_maint"]
LISTA_c_maint = []

println("Parameter c_maint $(length(LISTA_c_maint))")

M_c_maint = Dict{String,Array{Float64}}();

for (k,v) in Dict(Parametros_GSA)
    if k[1] == "c_maint" && in(k[2],LISTA_c_maint)
        d1 = Uniform(Parametros_GSA[k][2],Parametros_GSA[k][3])
        M_c_maint[k[2]] = [rand(d1,1)[1] for i=1:NS]
    end
end

for i in TECHNOLOGIES
    try
    M_c_maint[i]
    catch error
       if isa(error, KeyError)
           M_c_maint[i] = [c_maint[i] for j=1:NS]; # default value
       end
    end
end

##################### tau 

M_tau = Dict{String,Array{Float64}}();

for i in TECHNOLOGIES
    try
    M_tau[i]
    catch error
       if isa(error, KeyError)
           M_tau[i] = [M_i_rate["i_rate"][j]*(1 + M_i_rate["i_rate"][j])^M_lifetime[i][j]/(((1 + M_i_rate["i_rate"][j])^M_lifetime[i][j])-1) for j=1:NS]; # default value
       end
    end
end

#f_max = 

#f_min =

#f_ref = 

############################################ II stage ##############################################

##################### end_uses_demand_year

LISTA_demand = [(k[2],k[3]) for (k,v) in Parametros_GSA_big if k[1]=="end_uses_demand_year"]
#LISTA_demand = []

println("Parameter end_uses_demand_year $(length(LISTA_demand))")

M_end_uses_demand_year = Dict{Tuple{String,String},Array{Float64}}()

for (k,v) in Dict(Parametros_GSA_big)
    if k[1] == "end_uses_demand_year" && in((k[2],k[3]),LISTA_demand)
        d1 = Uniform(Parametros_GSA_big[k][2],Parametros_GSA_big[k][3])
        M_end_uses_demand_year[k[2],k[3]] = [rand(d1,1)[1] for i=1:NS]
    end
end

for i in END_USES_INPUT, s in SECTORS
    try
    M_end_uses_demand_year[i,s]
    catch error
       if isa(error, KeyError)
           M_end_uses_demand_year[i,s] = [end_uses_demand_year[i,s] for j=1:NS]; # default value
       end
    end
end

M_end_uses_input = Dict(i => [sum([M_end_uses_demand_year[i,s][k] for s in SECTORS if end_uses_demand_year[i,s] > 0.0]) for k = 1:NS] for i in END_USES_INPUT);

##################### layers_in_out

M_layers_in_out = Dict{Tuple{String,String},Array{Float64}}()

LISTA_layers_in_out = [(k[2],k[3]) for (k,v) in Parametros_GSA_big if k[1]=="layers_in_out"]
LISTA_layers_in_out = []

println("Parameter Layers_in_out $(length(LISTA_layers_in_out))")

##################### c_p_t 

M_c_p_t = Dict{Tuple{String,Int64},Array{Float64}}();

LISTA_c_p_t = ["PV","WIND","HYDRO_DAM","HYDRO_RIVER","DEC_SOLAR"]
#LISTA_c_p_t =[]

println("Parametro Cpt $(length(LISTA_c_p_t))")

for (k,v) in Dict(Parametros_GSA)
    if k[1] == "c_p_t" && in(k[2],LISTA_c_p_t)
            d1 = Uniform(Parametros_GSA[k][2],Parametros_GSA[k][3])
            aux = [rand(d1,1)[1] for i=1:NS]
        for i = PERIODS
            M_c_p_t[k[2],i] = aux*c_p_t[k[2],i]
	    if k[2] == "HYDRO_DAM"
		M_c_p_t["NEW_HYDRO_DAM",i] = aux*c_p_t[k[2],i]
            end
	    if k[2] == "HYDRO_RIVER"
		M_c_p_t["NEW_HYDRO_RIVER",i] = aux*c_p_t[k[2],i]
            end
        end
    end
end


for i in TECHNOLOGIES, t in PERIODS
    try
    M_c_p_t[i,t]
    catch error
       if isa(error, KeyError)
           M_c_p_t[i,t] = [c_p_t[i,t] for j=1:NS]; # default value
       end
    end
end


##################### c_p

M_c_p = Dict{String,Array{Float64}}();

LISTA_c_p = ["NUCLEAR","CCGT","CCGT_CCS","COAL_US","COAL_IGCC",
             "COAL_US_CCS","COAL_IGCC_CCS","GEOTHERMAL","IND_COGEN_GAS",
             "IND_COGEN_WOOD","IND_COGEN_WASTE","IND_BOILER_GAS",
             "IND_BOILER_WOOD","IND_BOILER_OIL","IND_BOILER_COAL",
             "IND_BOILER_WASTE","IND_DIRECT_ELEC","DHN_HP_ELEC",
             "DHN_COGEN_GAS","DHN_COGEN_WOOD","DHN_COGEN_WASTE",
             "DHN_BOILER_GAS","DHN_BOILER_WOOD","DHN_BOILER_OIL",
             "DHN_DEEP_GEO","POWER2GAS_1","POWER2GAS_2","H2_ELECTROLYSIS",
             "H2_NG","H2_BIOMASS","GASIFICATION_SNG","PYROLYSIS"]
#LISTA_c_p =[]

println("Parameter Cp $(length(LISTA_c_p))")

for (k,v) in Dict(Parametros_GSA)
    if k[1] == "c_p" && in(k[2],LISTA_c_p)
        d1 = Uniform(Parametros_GSA[k][2],Parametros_GSA[k][3])
        M_c_p[k[2]] = [rand(d1,1)[1] for i=1:NS]
    end
end

for i in TECHNOLOGIES
    try
    M_c_p[i]
    catch error
       if isa(error, KeyError)
           M_c_p[i] = [c_p[i] for j=1:NS]; # default value
       end
    end
end

##################### c_op

M_c_op =  Dict{Tuple{String,Int64},Array{Float64}}();
                
LISTA_c_op = ["NG","COAL","LFO","DIESEL","GASOLINE","ELECTRICITY","WOOD","URANIUM"];
#LISTA_c_op = []
                				
println("Parameter Cop $(length(LISTA_c_op))")

cop_error = Dict{String,Array{Float64}}(
"ELECTRICITY" => [0.0 for i=1:NS], 
"GASOLINE"  => [0.0 for i=1:NS],
#"BIOETHANOL"  => [],
"DIESEL"  => [0.0 for i=1:NS],
#"BIODIESEL"  => [],
"LFO"  => [0.0 for i=1:NS],
"NG"  => [0.0 for i=1:NS],
"NG_CCS"  => [0.0 for i=1:NS],
#"SNG"  => [],
"COAL"  => [0.0 for i=1:NS],
"COAL_CCS"  => [0.0 for i=1:NS],
"URANIUM"  => [0.0 for i=1:NS]
)

for (k,v) in Dict(Parametros_GSA)
    if k[1] == "c_op" && in(k[2],LISTA_c_op)
        if in(k[2],["WOOD"])
            d1 = Uniform(Parametros_GSA[k][2],Parametros_GSA[k][3])
            aux = [rand(d1,1)[1] for i=1:NS]
            for i = PERIODS
	        M_c_op[k[2],i] = aux
	    end            
        else
	    d1 = Uniform(Parametros_GSA[k][2],Parametros_GSA[k][3])
		cop_error[k[2]] = [rand(d1,1)[1] for i=1:NS]
            for t = PERIODS
	            M_c_op[k[2],t] = cop_error[k[2]]
	    end
            if in(k[2],["NG"])
                for t=PERIODS
                    M_c_op["NG_CCS",t] = M_c_op["NG",t]
					M_c_op["SNG",t] = M_c_op["NG",t]
                end
            end
            if in(k[2],["COAL"])
                for t=PERIODS
                    M_c_op["COAL_CCS",t] = M_c_op["COAL",t]
                end
            end
			if in(k[2],["GASOLINE"])
                for t=PERIODS
                    M_c_op["BIOETHANOL",t] = M_c_op["GASOLINE",t]
                end
            end
			if in(k[2],["DIESEL"])
                for t=PERIODS
                    M_c_op["BIODIESEL",t] = M_c_op["DIESEL",t]
                end
            end
        end
    end
end

for i in RESOURCES, t in PERIODS
    try
    M_c_op[i,t]
    catch error
       if isa(error, KeyError)
           M_c_op[i,t] = [c_op[i,t] for j=1:NS]; # default value
       end
    end
end

##################### avail

LISTA_avail= ["WOOD","WASTE"]
#LISTA_avail = []

println("Parameter avail $(length(LISTA_avail))")

M_avail = Dict{String,Array{Float64}}();

for (k,v) in Dict(Parametros_GSA)
    if k[1] == "avail" && in(k[2],LISTA_avail)
        d1 = Uniform(Parametros_GSA[k][2],Parametros_GSA[k][3])
        M_avail[k[2]] = [rand(d1,1)[1] for i=1:NS]
    end
end

for i in RESOURCES
    try
    M_avail[i]
    catch error
       if isa(error, KeyError)
           M_avail[i] = [avail[i] for j=1:NS]; # default value
       end
    end
end

##################### loss 

M_loss_coeff = Dict{String,Array{Float64}}();

LISTA_losses = ["ELECTRICITY","HEAT_LOW_T_DHN"]
#LISTA_losses = []

println("Parameter %loss $(length(LISTA_losses))")

for (k,v) in Dict(Parametros_GSA)
    if k[1] == "loss_coeff" && in(k[2],LISTA_losses)
        d1 = Uniform(Parametros_GSA[k][2],Parametros_GSA[k][3])
        M_loss_coeff[k[2]] = [rand(d1,1)[1] for i=1:NS]
    end
end

for i in END_USES_TYPES
    try
    M_loss_coeff[i]
    catch error
       if isa(error, KeyError)
           M_loss_coeff[i] = [loss_coeff[i] for j=1:NS]; # default value
       end
    end
end

##################### %peak_DHN (un solo parametro)

M_peak_dhn_factor = Dict{String,Array{Float64}}();

LISTA_peak = []

println("Parameter %PeakDHN $(length(LISTA_peak))")

##################### fmin_%

M_fmin_perc = Dict{String,Array{Float64}}();

LISTA_fmin_perc = []

println("Parameter %fmin $(length(LISTA_fmin_perc))")

##################### fmax_%

M_fmax_perc = Dict{String,Array{Float64}}();

LISTA_fmax_perc = []

println("Parameter %fmax $(length(LISTA_fmax_perc))")

##################### lighting_month

M_lighting_month = Dict{Int64,Array{Float64}}();

LISTA_lighting_month = []

println("Parameter %lighting $(length(LISTA_lighting_month))")

##################### heating_month

M_heating_month = Dict{Int64,Array{Float64}}();

LISTA_heating_month = []

println("Parametro %sh $(length(LISTA_heating_month))")

count_param_uncer = length(LISTA_irate)+length(LISTA_lifetime)+length(LISTA_c_inv)+length(LISTA_c_maint)+length(LISTA_demand) + length(LISTA_layers_in_out) + length(LISTA_c_p_t) + length(LISTA_c_p) + length(LISTA_c_op) + length(LISTA_avail) + length(LISTA_losses) + length(LISTA_peak) + length(LISTA_fmin_perc) + length(LISTA_fmax_perc) + length(LISTA_lighting_month) + length(LISTA_heating_month);

println("Number of parameters uncertainty = $(count_param_uncer)")

