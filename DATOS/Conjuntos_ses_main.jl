
PERIODS = [i for i=1:12]; # time periods
SECTORS = ["HOUSEHOLDS","SERVICES","INDUSTRY","TRANSPORTATION"];
END_USES_INPUT = ["ELECTRICITY","LIGHTING","HEAT_HIGH_T","HEAT_LOW_T_SH","HEAT_LOW_T_HW","MOBILITY_PASSENGER","MOBILITY_FREIGHT"];
END_USES_CATEGORIES  = ["ELECTRICITY","HEAT_HIGH_T","HEAT_LOW_T","MOBILITY_PASSENGER","MOBILITY_FREIGHT"];
RESOURCES = ["ELECTRICITY","GASOLINE","DIESEL","BIOETHANOL","BIODIESEL","LFO","LNG","NG","NG_CCS","SNG","WOOD","COAL","COAL_CCS","URANIUM","WASTE","H2","ELEC_EXPORT"];
BIOFUELS = ["BIOETHANOL","BIODIESEL","SNG"];
EXPORT = ["ELEC_EXPORT"];

AUX = ["TRAMWAY_TROLLEY","BUS_COACH_DIESEL","BUS_COACH_HYDIESEL","BUS_COACH_CNG_STOICH","BUS_COACH_FC_HYBRIDH2","TRAIN_PUB","CAR_GASOLINE","CAR_DIESEL","CAR_NG","CAR_HEV","CAR_PHEV","CAR_BEV","CAR_FUEL_CELL","TRAIN_FREIGHT","TRUCK"]

END_USES_TYPES_OF_CATEGORY = Dict{String,Array{String}}(
    "ELECTRICITY" => ["ELECTRICITY"],
    "HEAT_HIGH_T" => ["HEAT_HIGH_T"],
    "HEAT_LOW_T"  => ["HEAT_LOW_T_DHN","HEAT_LOW_T_DECEN"],
    "MOBILITY_PASSENGER" => ["MOB_PUBLIC","MOB_PRIVATE"],
    "MOBILITY_FREIGHT" => ["MOB_FREIGHT_RAIL","MOB_FREIGHT_ROAD"]
);

TECHNOLOGIES_OF_END_USES_TYPE = Dict{String,Array{String}}(
    "ELECTRICITY" => ["NUCLEAR","CCGT","CCGT_CCS","COAL_US","COAL_IGCC","COAL_US_CCS","COAL_IGCC_CCS","PV","WIND","HYDRO_DAM","NEW_HYDRO_DAM","HYDRO_RIVER","NEW_HYDRO_RIVER","GEOTHERMAL"],
    "HEAT_HIGH_T" => ["IND_COGEN_GAS","IND_COGEN_WOOD","IND_COGEN_WASTE","IND_BOILER_GAS","IND_BOILER_WOOD","IND_BOILER_OIL","IND_BOILER_COAL","IND_BOILER_WASTE","IND_DIRECT_ELEC"],
    "HEAT_LOW_T_DHN" => ["DHN_HP_ELEC","DHN_COGEN_GAS","DHN_COGEN_WOOD","DHN_COGEN_WASTE","DHN_BOILER_GAS","DHN_BOILER_WOOD","DHN_BOILER_OIL","DHN_DEEP_GEO"],
    "HEAT_LOW_T_DECEN" => ["DEC_HP_ELEC","DEC_THHP_GAS","DEC_COGEN_GAS","DEC_COGEN_OIL","DEC_ADVCOGEN_GAS","DEC_ADVCOGEN_H2","DEC_BOILER_GAS","DEC_BOILER_WOOD","DEC_BOILER_OIL","DEC_SOLAR","DEC_DIRECT_ELEC"],
    "MOB_PUBLIC" => ["TRAMWAY_TROLLEY","BUS_COACH_DIESEL","BUS_COACH_HYDIESEL","BUS_COACH_CNG_STOICH","BUS_COACH_FC_HYBRIDH2","TRAIN_PUB"],
    "MOB_PRIVATE" => ["CAR_GASOLINE","CAR_DIESEL","CAR_NG","CAR_HEV","CAR_PHEV","CAR_BEV","CAR_FUEL_CELL"],
    "MOB_FREIGHT_RAIL" => ["TRAIN_FREIGHT"],
    "MOB_FREIGHT_ROAD" => ["TRUCK"]
);

STORAGE_TECH = ["PUMPED_HYDRO","POWER2GAS"];

INFRASTRUCTURE = ["EFFICIENCY","DHN","GRID","POWER2GAS_1","POWER2GAS_2","POWER2GAS_3","H2_ELECTROLYSIS","H2_NG","H2_BIOMASS","GASIFICATION_SNG","PYROLYSIS"];

COGEN = ["IND_COGEN_GAS","IND_COGEN_WOOD","IND_COGEN_WASTE","DHN_COGEN_GAS","DHN_COGEN_WOOD","DHN_COGEN_WASTE","DEC_COGEN_GAS","DEC_COGEN_OIL","DEC_ADVCOGEN_GAS","DEC_ADVCOGEN_H2"];

BOILERS = ["IND_BOILER_GAS","IND_BOILER_WOOD","IND_BOILER_OIL","IND_BOILER_COAL","IND_BOILER_WASTE","DHN_BOILER_GAS","DHN_BOILER_WOOD","DHN_BOILER_OIL","DEC_BOILER_GAS","DEC_BOILER_WOOD","DEC_BOILER_OIL"];

END_USES_TYPES = [];

for i in END_USES_CATEGORIES, j in END_USES_TYPES_OF_CATEGORY[i]
    push!(END_USES_TYPES,j)
end

LAYERS = union(setdiff(setdiff(RESOURCES,BIOFUELS),EXPORT),END_USES_TYPES);

Aux_set = [];
for i in END_USES_TYPES, j in TECHNOLOGIES_OF_END_USES_TYPE[i]
    push!(Aux_set,j)
end

TECHNOLOGIES = union(Aux_set,STORAGE_TECH,INFRASTRUCTURE);

TECHNOLOGIES_OF_END_USES_CATEGORY = Dict(i => [k for j in END_USES_TYPES_OF_CATEGORY[i] for k in TECHNOLOGIES_OF_END_USES_TYPE[j]] for i in END_USES_CATEGORIES);

function crearDict(conjunto)
    cont = 1;
    dict = Dict();
    for i in conjunto
        dict[i] = cont
        cont +=1;
    end
    return dict
end

RESO1 = ["GASOLINE","DIESEL","LFO","LNG","NG","NG_CCS","WOOD","COAL","COAL_CCS","URANIUM","WASTE","H2"]
RESO2 = setdiff(LAYERS,RESO1)

LAYERSD = crearDict(LAYERS);
TECHNOLOGIESD = crearDict(TECHNOLOGIES);
RESOURCESD = crearDict(RESOURCES);
TECHNOLOGIESdiffINFRASTRUCTURED = crearDict(setdiff(TECHNOLOGIES,INFRASTRUCTURE));
TECHNOLOGIES_OF_END_USES_TYPEHDD = crearDict(TECHNOLOGIES_OF_END_USES_TYPE["HEAT_LOW_T_DECEN"]);
STORAGE_TECHD = crearDict(STORAGE_TECH);
END_USES_TYPESD = crearDict(END_USES_TYPES);
op_strategy_mob_privateD = crearDict(union(TECHNOLOGIES_OF_END_USES_CATEGORY["MOBILITY_PASSENGER"],TECHNOLOGIES_OF_END_USES_CATEGORY["MOBILITY_FREIGHT"]));


