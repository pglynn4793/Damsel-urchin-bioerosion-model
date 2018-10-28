
#Cellular Automata Program
#
#
#This program uses a spatial model represented by cells to project
#the future state of a pocilloporid reef subjected to urchin bioerosion pressures
#
#Open files and load data structures
#
cellular_automata_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/CellularAutomataDataModel.csv"
reef_zone_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/Reef_Zone_Data.csv"
cellular_automata_params_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/Cellular_Automata_Params.csv"
homing_freq_vec_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/HomingDistanceFreqVec.csv"
aggregation_freq_vec_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/AggregationDistanceFreqVec.csv"
damselfish_resp_by_length_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/DamselfishResponsebyLength.csv"
damselfish_size_dist_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/DamselfishSizeDistData.csv"
circle_packing_sum_squares_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/CirclePackingClusterSumSquares.csv"
join_aggregation_freq_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/JoinAggregationFrequencyData.csv"
high_density_urchin_dist_sampling_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/HighDensityDistSamplingVec.csv"
mod_density_urchin_dist_sampling_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/ModDensityDistSamplingVec.csv"
agonistic_response_sampling_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/ExtendedAgonisticResponseSamplingVec.csv"
next_cell_transition_map_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/NextCellTransitionMap.csv"
symbol_to_coordinate_map_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/SymbolCoordinateMap.csv"
border_cell_id_zone1_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/BorderCellID_Zone1.csv"
border_cell_id_zone2_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/BorderCellID_Zone2.csv"
day_or_night_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/DayorNightSampleData.csv"
cellular_automata_output_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/CellularAutomataOutputData"
next_cell_route_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/NextCellRouteMap.csv"
pick_route_id_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/PickRouteIDFrequencyData.csv"
cellular_automata_summary_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/CellularAutomataSummaryData"
cellular_automata_simulation_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/CellularAutomataSimulationData"
proportional_disaggregation_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/ProportionalDisaggregationData.csv"
disaggregation_frequency_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/DisaggregationFrequencyData.csv"
damsel_size_by_dist_selection_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/DamselSizebyDistSelectionData.csv"
damsel_density_sampling_file_paths_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/DamselDensitySamplingFilePaths.csv"
damsel_responses_sampling_file_paths_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/DamselResponseSamplingFilePaths.csv"
num_aggregation_subgroups_sample_data_file_path<-"~/Desktop/Peter Docs/AgentBasedCarbonateModel/NumAggregationSubgroupsSampleData.csv"
#
cellular_automata_data<-read.csv(cellular_automata_file_path, stringsAsFactors = FALSE, header = T)
reef_zone_data<-read.csv(reef_zone_data_file_path, stringsAsFactors = FALSE, header = T)
cellular_automata_params<-read.csv(cellular_automata_params_file_path, stringsAsFactors = FALSE, header = T)
homing_dist_freq_data<-read.csv(homing_freq_vec_file_path, stringsAsFactors = FALSE, header = T)
aggregation_dist_freq_data<-read.csv(aggregation_freq_vec_file_path, stringsAsFactors = FALSE, header = T)
damselfish_resp_by_length_data<-read.csv(damselfish_resp_by_length_file_path, stringsAsFactors = FALSE, header = T)
damselfish_size_dist_data<-read.csv(damselfish_size_dist_data_file_path, stringsAsFactors = FALSE, header = T)
circle_packing_sum_squares_data<-read.csv(circle_packing_sum_squares_data_file_path, stringsAsFactors = FALSE, header = T)
join_aggregation_freq_data<-read.csv(join_aggregation_freq_data_file_path, stringsAsFactors = FALSE, header = T)
high_density_urchin_dist_sampling_data<-read.csv(high_density_urchin_dist_sampling_data_file_path, stringsAsFactors = FALSE, header = T)
mod_density_urchin_dist_sampling_data<-read.csv(mod_density_urchin_dist_sampling_data_file_path, stringsAsFactors = FALSE, header = T)
agonistic_response_sampling_data<-read.csv(agonistic_response_sampling_data_file_path, stringsAsFactors = FALSE, header = T)
next_cell_transition_map_data<-read.csv(next_cell_transition_map_data_file_path, stringsAsFactors = FALSE, header = T)
symbol_to_coordinate_map_data<-read.csv(symbol_to_coordinate_map_data_file_path, stringsAsFactors = FALSE, header = T)
border_cell_id_zone1_data<-read.csv(border_cell_id_zone1_data_file_path, stringsAsFactors = FALSE, header = T)
border_cell_id_zone2_data<-read.csv(border_cell_id_zone2_data_file_path, stringsAsFactors = FALSE, header = T)
day_or_night_data<-read.csv(day_or_night_data_file_path, stringsAsFactors = FALSE, header = T)
next_cell_route_data<-read.csv(next_cell_route_data_file_path, stringsAsFactors = FALSE, header = T)
pick_route_id_data<-read.csv(pick_route_id_data_file_path, stringsAsFactors = FALSE, header = T)
proportional_disaggregation_data<-read.csv(proportional_disaggregation_data_file_path, stringsAsFactors = FALSE, header = T)
disaggregation_frequency_data<-read.csv(disaggregation_frequency_data_file_path, stringsAsFactors = FALSE, header = T)
damsel_size_by_dist_selection_data<-read.csv(damsel_size_by_dist_selection_data_file_path, stringsAsFactors = FALSE, header = T)
damsel_density_sampling_file_paths_data<-read.csv(damsel_density_sampling_file_paths_data_file_path, stringsAsFactors = FALSE, header = T)
damsel_responses_sampling_file_paths_data<-read.csv(damsel_responses_sampling_file_paths_data_file_path, stringsAsFactors = FALSE, header = T)
num_aggregation_subgroups_sample_data<-read.csv(num_aggregation_subgroups_sample_data_file_path, stringsAsFactors = FALSE, header = T)
#
pick_route_id_data_vec<-pick_route_id_data[,1]
day_or_night_data_vec<-day_or_night_data[,1]
border_cell_id_zone1_vector<-border_cell_id_zone1_data[,1]
border_cell_id_zone2_vector<-border_cell_id_zone2_data[,1]
agonistic_response_sampling_vec<-agonistic_response_sampling_data[,1]
high_density_urchin_dist_sampling_vec<-high_density_urchin_dist_sampling_data[,1]
mod_density_urchin_dist_sampling_vec<-mod_density_urchin_dist_sampling_data[,1]
damselfish_size_vec<-damselfish_size_dist_data[,1]
proportional_disaggregation_data_vec<-proportional_disaggregation_data[,1]
disaggregation_frequency_data_vec<-disaggregation_frequency_data[,1]
disaggregation_frequency_data_vec
num_aggregation_subgroups_sample_data_vec<-num_aggregation_subgroups_sample_data[,1]
num_aggregation_subgroups_sample_data_vec
#
#
#initialize homing distance and aggregation distance frequency vectors 
#
join_aggregation_freq_data_vec<-join_aggregation_freq_data[,1]
homing_dist_freq_data_vec<-homing_dist_freq_data[,1]
aggregation_dist_freq_data_vec<-aggregation_dist_freq_data[,1]
#aggregation_dist_freq_data_vec
#homing_dist_freq_data_vec
#
#Randomly sample from urchin distribution specified for each urchin zone in configuration file
#to initialize cells with initial urchin densities
#
#Need to loop over all cells
#
#
total_sum_squares_data_dim_vec<-dim(circle_packing_sum_squares_data)
total_sum_squares_data_size<-total_sum_squares_data_dim_vec[1]
#total_sum_squares_data_size
cellular_automata_data_dim_vec<-dim(cellular_automata_data)
cellular_automata_num_cells<-cellular_automata_data_dim_vec[1]
#cellular_automata_num_cells
#cellular_automata_data
food_value_cell_vec<-numeric(cellular_automata_num_cells)
#food_value_cell_vec
#
reef_zone_data_dim_vec<-dim(reef_zone_data)
#reef_zone_data_dim_vec
reef_num_zones<-reef_zone_data_dim_vec[1]
#reef_num_zones
#
#
#Obtain parameters
#
#time periods
#
total_num_days<-cellular_automata_params$Num_Days[1]
time_period_intervals<-cellular_automata_params$Time_period_intervals[1]
#total_num_days
#time_period_intervals
#Urchin daily bioerosion parameter
urchin_daily_bioerosion_rate<-cellular_automata_params$Urchin_Mean_Ind_Daily_Bioerosion_g[1]
#
#coralline algae parameters
coralline_algae_annual_growth_rate<-cellular_automata_params$Corr_Algae_Carbonate_Dep_Rate_kg_yr[1]
#coralline_algae_annual_growth_rate
coralline_algae_daily_growth_rate_g_sq_meter<-(coralline_algae_annual_growth_rate/365.0)*1000
coralline_algae_mean_density<-cellular_automata_params$Corr_Algae_Mean_Density_g_cm2[1]
#coralline_algae_mean_density
urchin_daily_bioeorosion_rate_corr_algae<-cellular_automata_params$Urchin_Daily_Corr_Algae_Bioerosion_g[1]
#
#Live coral parameters
#
live_coral_annual_carbon_prod_kg_yr_sq_meter<-cellular_automata_params$Live_Coral_Carbon_Prod_kg_yr[1]
#live_coral_annual_carbon_prod_kg_yr_sq_meter
live_coral_daily_growth_rate_g_sq_meter<-(live_coral_annual_carbon_prod_kg_yr_sq_meter/365)*(1000.0)
#live_coral_daily_growth_rate_g_sq_meter
live_coral_mean_density<-cellular_automata_params$Live_Coral_Mean_Density_g_cm2[1]
#live_coral_mean_density
urchin_daily_bioerosion_rate_live_coral<-cellular_automata_params$Urchin_Daily_Live_Coral_Bioerosion_g[1]
#urchin_daily_bioerosion_rate_live_coral
live_coral_cover_growth_rate<-cellular_automata_params$Live_Coral_Cover_Growth_Rate[1]
#
#Algal lawn parameters
#
algal_lawn_daily_growth_rate_g_m2<-cellular_automata_params$Algal_Lawn_Annual_Growth_Rate_g_m2[1]/365
algal_lawn_biomass_g_per_cm2<-cellular_automata_params$Algal_Lawn_Biomass_g_per_cm2[1]
#
#Minimum urchin density at which aggregations form
min_urchin_aggreg_density<-cellular_automata_params$Min_Urchin_Aggreg_Density[1]
#
#Vector containing indeces of randomly selected points from a cluster with aggregating urchins
#
random_points_index_vec<-numeric(cellular_automata_params$Num_Random_Points[1])
#
#Summary and detail data file id offsets
#
file_id_offset<-cellular_automata_params$Simulation_Run_File_ID_Offset[1]
#
#Number of random points selected from cluster of aggregating urchins that will
#be used to calculate distance with urchins that are not aggregating in current cell
#
num_random_points<-cellular_automata_params$Num_Random_Points[1]
#
#Maximum number of urchins allowed in a cell
#
max_num_urchins_per_cell<-cellular_automata_params$Max_Number_Urchins_Cell[1]
#calculate food preference proportion factors
food_pref_live_coral_factor<-(cellular_automata_params$Urchin_Daily_Live_Coral_Bioerosion_g[1])/(cellular_automata_params$Urchin_Daily_Live_Coral_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Corr_Algae_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Protected_Algal_Lawn_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Unprotected_Algal_Lawn_Bioerosion_g[1])
food_pref_corr_algae_factor<-(cellular_automata_params$Urchin_Daily_Corr_Algae_Bioerosion_g[1])/(cellular_automata_params$Urchin_Daily_Live_Coral_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Corr_Algae_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Protected_Algal_Lawn_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Unprotected_Algal_Lawn_Bioerosion_g[1])
food_pref_protected_algal_lawn_factor<-(cellular_automata_params$Urchin_Daily_Protected_Algal_Lawn_Bioerosion_g[1])/(cellular_automata_params$Urchin_Daily_Live_Coral_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Corr_Algae_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Protected_Algal_Lawn_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Unprotected_Algal_Lawn_Bioerosion_g[1])
food_pref_unprotected_algal_lawn_factor<-(cellular_automata_params$Urchin_Daily_Unprotected_Algal_Lawn_Bioerosion_g[1])/(cellular_automata_params$Urchin_Daily_Live_Coral_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Corr_Algae_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Protected_Algal_Lawn_Bioerosion_g[1] + cellular_automata_params$Urchin_Daily_Unprotected_Algal_Lawn_Bioerosion_g[1])
#
#
#
urchin_density_data_zone_1<-read.csv(reef_zone_data$Init_Urchin_Distribution_Files[1], stringsAsFactors = FALSE, header = T)
urchin_density_data_zone_1_vec<-urchin_density_data_zone_1[,2]
#urchin_density_data_zone_1_vec
urchin_density_data_zone_2<-read.csv(reef_zone_data$Init_Urchin_Distribution_Files[2], stringsAsFactors = FALSE, header = T)
urchin_density_data_zone_2_vec<-urchin_density_data_zone_2[,2]
urchin_density_data_zone_3<-read.csv(reef_zone_data$Init_Urchin_Distribution_Files[3], stringsAsFactors = FALSE, header = T)
urchin_density_data_zone_3_vec<-urchin_density_data_zone_3[,2]
urchin_density_data_zone_4<-read.csv(reef_zone_data$Init_Urchin_Distribution_Files[4], stringsAsFactors = FALSE, header = T)
urchin_density_data_zone_4_vec<-urchin_density_data_zone_4[,2]
#urchin_density_data_zone_4_vec
#
#
#create data frame for transition of urchins to next cell
#
cell_id<-numeric(cellular_automata_num_cells)
current_urchin_density<-numeric(cellular_automata_num_cells)
number_input_urchins<-numeric(cellular_automata_num_cells)
number_output_urchins<-numeric(cellular_automata_num_cells)
next_cell_urchin_transition_data<-data.frame(cell_id,current_urchin_density,number_input_urchins,number_output_urchins)
next_cell_urchin_transition_data$cell_id<-c(1:cellular_automata_num_cells)
#
#create a data frame for outputing simulation data
#
num_rows<-total_num_days*cellular_automata_num_cells
cell_id<-numeric(num_rows)
time_id<-numeric(num_rows)
zone_id<-numeric(num_rows)
urchin_density_before<-numeric(num_rows)
net_carbonate<-numeric(num_rows)
damselfish_size<-numeric(num_rows)
aggregation_size<-numeric(num_rows)
num_urchin_ejections<-numeric(num_rows)
cellular_automata_output_data<-data.frame(cell_id,time_id,zone_id,urchin_density_before,net_carbonate,damselfish_size,aggregation_size,num_urchin_ejections)
#
#create data frame for outputing summary data on urchin density 
#
zone1_mean_urchin_density<-numeric(total_num_days)
zone2_mean_urchin_density<-numeric(total_num_days)
zone1_mean_damselfish_density<-numeric(total_num_days)
zone2_mean_damselfish_density<-numeric(total_num_days)
zone1_mean_damselfish_size<-numeric(total_num_days)
zone2_mean_damselfish_size<-numeric(total_num_days)
diff_mean_urchin_density<-numeric(total_num_days)
coral_patch_mean_urchin_density<-numeric(total_num_days)
zone1_exclude_patch_damselfish_urchin_density<-numeric(total_num_days)
zone2_exclude_patch_damselfish_urchin_density<-numeric(total_num_days)
cellular_automata_summary_data<-data.frame(zone1_mean_urchin_density,zone2_mean_urchin_density,diff_mean_urchin_density,zone1_mean_damselfish_density,zone2_mean_damselfish_density,zone1_mean_damselfish_size,zone2_mean_damselfish_size,coral_patch_mean_urchin_density,zone1_exclude_patch_damselfish_urchin_density,zone2_exclude_patch_damselfish_urchin_density)
#
#create data frame for outputing simulation-level summary data
#
num_simulation_runs<-cellular_automata_params$Num_Simulation_Runs[1]
#
simulation_id<-numeric(num_simulation_runs)
number_of_days<-numeric(num_simulation_runs)
zone1_mean_urchin_density<-numeric(num_simulation_runs)
zone2_mean_urchin_density<-numeric(num_simulation_runs)
aggregation_affinity_param<-numeric(num_simulation_runs)
max_size_urchin_aggregate<-numeric(num_simulation_runs)
mean_damselfish_density_zone1<-numeric(num_simulation_runs)
mean_damselfish_density_zone2<-numeric(num_simulation_runs)
mean_damselfish_size_zone1<-numeric(num_simulation_runs)
mean_damselfish_size_zone2<-numeric(num_simulation_runs)
coral_patch_init_mean_urchin_density<-numeric(num_simulation_runs)
coral_patch_final_mean_urchin_density<-numeric(num_simulation_runs)
exclude_patch_damselfish_init_zone1_mean_urchin_density<-numeric(num_simulation_runs)
exclude_patch_damselfish_final_zone1_mean_urchin_density<-numeric(num_simulation_runs)
exclude_patch_damselfish_init_zone2_mean_urchin_density<-numeric(num_simulation_runs)
exclude_patch_damselfish_final_zone2_mean_urchin_density<-numeric(num_simulation_runs)
#
cellular_automata_simulation_data<-data.frame(simulation_id,number_of_days,aggregation_affinity_param,max_size_urchin_aggregate,mean_damselfish_density_zone1,mean_damselfish_density_zone2,mean_damselfish_size_zone1,mean_damselfish_size_zone2,zone1_mean_urchin_density,zone2_mean_urchin_density,coral_patch_init_mean_urchin_density,coral_patch_final_mean_urchin_density,exclude_patch_damselfish_init_zone1_mean_urchin_density,exclude_patch_damselfish_final_zone1_mean_urchin_density,exclude_patch_damselfish_init_zone2_mean_urchin_density,exclude_patch_damselfish_final_zone2_mean_urchin_density)
#
#next_cell_urchin_transition_data
#
#Damselfish ejection checking function
#
DamselfishEjectionCheck<-function(x)
{
  #x - length of damselfish
  #this function will return T or F indication of whether damselfish
  #ejection is activated; will also check to ensure that an aggregation of 
  #urchins that occurs following a cell transition does not exceed maximum
  #specified parameter value
  #
  #
  #
  result<-"F"
  if (x <= 0.0)
  {
    result<-"F"
    return (result)
  }
  else
  {
   #damselfish is larger than 0 cm
   #perform day/night check
   #if night time, no ejection
    day_or_night_indicator<-sample(day_or_night_data_vec,1)
    if (day_or_night_indicator == 0)
    {
      #night time
      result<-"F"
      return(result)
    }
    else
    {
      #it is day time
      #now check if damselfish agonistic according to its size
      #
      agonism_freq<-damselfish_resp_by_length_data[damselfish_resp_by_length_data$Damselfish_length_cm==x,]$Percent_response
      if (agonism_freq <= 0.0)
      {
        #damselfish has no agonistic response
        result<-"F"
        return(result)
      }
      else
      {
        #damselfish may have agonistic response
        #produce a sampling vector to check for agonistic response
        #must have probability be a whole number
        agonism_freq<-round(agonism_freq)
        agonism_sample_vec<-numeric(100)
        #
        #repeat ones into agonism sample vec according to agonism frequency
        for (k in 1:agonism_freq)
        {
          agonism_sample_vec[k]<-1
        }
        #
        #randomly sample from agonism sampling vector
        #
        agonism_indicator<-sample(agonism_sample_vec,1)
        if (agonism_indicator == 0)
        {
          #damselfish is not displaying agonism
          result<-"F"
          return(result)
        }
        else
        {
          #damselfish is agonistic
          #now check if ejection occurs
          #
          urchin_ejection_indicator<-sample(agonistic_response_sampling_vec,1)
          if (urchin_ejection_indicator == "Y")
          {
            #urchin will be ejected
            result<-"T"
            return(result)
          }
          else
          {
            #urchin will not be ejected
            result<-"F"
            return(result)
          }#end if agonistic damselfish did not eject urchin
        }#end if urchin is agonistic    
      }#end if urchin may be agonistic
    }#end if it is day time
  }#end if damselfish not zero length
}#end function
#
#
#
NextCellTransition<-function(cur_cell,nxt_cell,urchin_density)
{
  #this function will return a list with the physical id of next cell, and a flag indicating if 
  #an ejection of an urchin (or aggregation) occurred
  #
  #inputs
  #
  #cur_cell - physical id of current cell
  #nxt_cell - symbol representing next cell proposed for transition
  #urchin_density - the amount of urchins proposed for transition to the next cell 
  #
  #
  #use symbol coordinate map to translate the symbol representing next-cell transition
  #to logical x and y coordinates
  nxt_cell_xcoord<-symbol_to_coordinate_map_data[symbol_to_coordinate_map_data$symbol==nxt_cell,]$x.coord
  nxt_cell_ycoord<-symbol_to_coordinate_map_data[symbol_to_coordinate_map_data$symbol==nxt_cell,]$y.coord
  #print("begin function call")
  #print("current cell")
  #print(cur_cell)
  #print("next cell symbol")
  #print(nxt_cell)
  #print("next cell x coord for input next cell")
  #print(nxt_cell_xcoord)
  #print("next cell y coord for input next cell")
  #print(nxt_cell_ycoord)
  #
  #Translate x and y coordinates of next-cell to a physical cell id
  #
  nxt_cell_physical_id<-next_cell_transition_map_data[next_cell_transition_map_data$cell_id==cur_cell&next_cell_transition_map_data$x_coord==nxt_cell_xcoord&next_cell_transition_map_data$y_coord==nxt_cell_ycoord,]$next_cell_id
  #print("next cell physical id")
  #print(nxt_cell_physical_id)
  #
  #Check if next cell physical id is zero; if so, we need to find a border cell
  #for transition, and check for damselfish ejection behavior
  #
  if (nxt_cell_physical_id == 0)
  {
    #
    #next cell is leaving lattice; allow urchins to leave lattice
    #
    #print("next cell physical id = 0")
    urchin_ejection_flag<-F
    return_list<-list("next_cell_id" = nxt_cell_physical_id, "damselfish_eject_flag" = urchin_ejection_flag)
    return(return_list)
  }#end if next-cell is zero
  else
  {
    #next cell is not zero
    #
    #use x & y coordinates for next cell to lookup
    #number of intermediate cells between current cell
    #and the next cell
    #
    number_int_cells_vec<-next_cell_route_data[next_cell_route_data$nxt_cell_xcoord==nxt_cell_xcoord&next_cell_route_data$nxt_cell_ycoord==nxt_cell_ycoord,]$num_int_cells
    number_int_cells<-number_int_cells_vec[1]
    if (number_int_cells == 0)
    {
      #no path; can check next cell directly for damselfish ejection
      urchin_ejection_flag<-DamselfishEjectionCheck(cellular_automata_data$DamselfishLength_cm[nxt_cell_physical_id])
      if (urchin_ejection_flag)
      {
        nxt_cell_physical_id<-cur_cell
      }  
      #create list for return of function
      if (nxt_cell_physical_id != cur_cell)
      {
        total_num_urchins_next_cell<-cellular_automata_data$Urchin_Density[cur_cell]+cellular_automata_data$Urchin_Density[nxt_cell_physical_id] + next_cell_urchin_transition_data$number_input_urchins[nxt_cell_physical_id]
        if (total_num_urchins_next_cell > max_num_urchins_per_cell)
        {
          #combined urchins exceed threshold; urchins will be ejected
          urchin_ejection_flag<-T
          nxt_cell_physical_id<-cur_cell
        }
      }
      return_list<-list("next_cell_id" = nxt_cell_physical_id, "damselfish_eject_flag" = urchin_ejection_flag)
      return(return_list)
    }# end if number of intermediate cells is zero, i.e., no path
    else
    {
      #a route exists
      route_id<-1
      #if more than one route exists, need to randomly select a route id
      #
      num_routes_vec<-next_cell_route_data[next_cell_route_data$nxt_cell_xcoord==nxt_cell_xcoord&next_cell_route_data$nxt_cell_ycoord==nxt_cell_ycoord,]$num_routes
      num_routes<-num_routes_vec[1]
      #print("number of routes")
      #print(num_routes)
      if (num_routes > 1)
      {
        #need to randomly select a route id
        #
        route_id<-sample(pick_route_id_data_vec,1)
      }#end if number of routes > 1
      #loop over number of cells in route
      #check for damselfish ejection at each cell
      #if rejected at a cell, set next cell for transition to current cell
      #if not rejected on route all the way to the destination, return the next cell
      #
      total_cells_in_route<-number_int_cells + 1
      #print("total number intermediate cells")
      #print(number_int_cells)
      for (i in 1:total_cells_in_route)
      {
        #loop over number of cells in route and check for damselfish ejection
        #
        route_cell_physical_id<-nxt_cell_physical_id
        if (i < total_cells_in_route)
        {
          #next cell is not end of route
          #find logical coordinates for next cell from route map
          #and then find physical cell id
          #
          #print("current cell id")
          #print(cur_cell)
          #print("i")
          #print(i)
          #print("next cell xcoord")
          #print(nxt_cell_xcoord)
          #print("next cell ycoord")
          #print(nxt_cell_ycoord)
          #print("route id")
          #print(route_id)
          logical_xcoord<-next_cell_route_data[next_cell_route_data$nxt_cell_xcoord==nxt_cell_xcoord&next_cell_route_data$nxt_cell_ycoord==nxt_cell_ycoord&next_cell_route_data$route_id==route_id&next_cell_route_data$route_cell_id==i,]$int_x_coord
          logical_ycoord<-next_cell_route_data[next_cell_route_data$nxt_cell_xcoord==nxt_cell_xcoord&next_cell_route_data$nxt_cell_ycoord==nxt_cell_ycoord&next_cell_route_data$route_id==route_id&next_cell_route_data$route_cell_id==i,]$int_y_coord
          #
          #find physical cell id
          #
          route_cell_physical_id<-next_cell_transition_map_data[next_cell_transition_map_data$cell_id==cur_cell&next_cell_transition_map_data$x_coord==logical_xcoord&next_cell_transition_map_data$y_coord==logical_ycoord,]$next_cell_id
          #print("next cell is not end of route")
          #print(route_cell_physical_id)
        }#end if current route cell is not destination
        #
        #check to see if damselfish ejection at current route cell
        #
        #print("route cell physical id")
        #print(route_cell_physical_id)
        if (cellular_automata_data$DamselfishLength_cm[route_cell_physical_id] > 0.0)
        {
          #
          #a damselfish is located in next-cell; check for damselfish ejection
          urchin_ejection_flag<-DamselfishEjectionCheck(cellular_automata_data$DamselfishLength_cm[route_cell_physical_id])
          if (urchin_ejection_flag)
          {
            #urchin has been ejected; set next-cell to current cell
            #create list for return of function
            return_list<-list("next_cell_id" = cur_cell, "damselfish_eject_flag" = urchin_ejection_flag)
            return(return_list)
          }#end if urchin is ejected
          #
        }#end if damselfish is non-zero length
        #
      }#end for loop on looping over number of cells in route
      #if reached end of for loop, this means no ejections; can transition to next-cell
      #Need to check to see if combination of urchins exceeds threshold
      #
      urchin_ejection_flag<-F
      total_num_urchins_next_cell<-cellular_automata_data$Urchin_Density[cur_cell]+cellular_automata_data$Urchin_Density[route_cell_physical_id]
      if (total_num_urchins_next_cell >max_num_urchins_per_cell)
      {
        #combined urchins exceed threshold; urchins will be ejected
        urchin_ejection_flag<-T
        route_cell_physical_id<-cur_cell
      }
      return_list<-list("next_cell_id" = route_cell_physical_id, "damselfish_eject_flag" = urchin_ejection_flag)
      return(return_list)
      
    }#end if number of intermediate cells > 0   
  }#end if next cell id is zero
  #
}#end function
#Loop over total number of days specified in input parameter
#
#
#cellular_automata_num_cells<-10
for (curr_simulation_run in 1:num_simulation_runs)
{
  damselfish_length_zone_1<-0
  damselfish_length_zone_2<-0
  damselfish_count_zone_1<-0
  damselfish_count_zone_2<-0
  
  for (i in 1:cellular_automata_num_cells)
  {
    #
    #Read record for current cell and determine zone
    #
    zone_id<-cellular_automata_data$Reef_Zone[i]
    #print("zone id")
    #print(zone_id)
    #
    #sample from urchin distribution assigned to reef zone 
    #
    if (zone_id == 1)
    {
      urchin_density_vec<-urchin_density_data_zone_1_vec
    }
    if (zone_id == 2)
    {
      urchin_density_vec<-urchin_density_data_zone_2_vec
    }
    if (zone_id == 3)
    {
      urchin_density_vec<-urchin_density_data_zone_3_vec
    }
    if (zone_id == 4)
    {
      urchin_density_vec<-urchin_density_data_zone_4_vec
    }
    urchin_density_mean<-reef_zone_data$Urchin_Density_Mean[zone_id]
    #randomly sample urchin density from distribution and initialize urchin density
    #value for current cell
    #
    urchin_density_vec_size<-length(urchin_density_vec)
    urchin_density_vec_indx<-sample(urchin_density_vec_size,1)
    urchin_density_init_value<-urchin_density_vec[urchin_density_vec_indx]
    #
    if (urchin_density_init_value < 0.0)
    {
      #adjust urchin density using mean urchin density for current distribution
      urchin_density_init_value<-(urchin_density_mean - urchin_density_init_value) + urchin_density_mean
    }
    cellular_automata_data$Urchin_Density[i]<-urchin_density_init_value
    #
    #Populate cells that are eligible for damselfish 
    #
    #
    #initialize damselfish length from previous simulation run before determining if cell eligible for damselfish 
    cellular_automata_data$DamselfishLength_cm[i]<-0.0
    #
    if ((cellular_automata_data$Dist_ID_To_Coral_Patch[i] >= 1) & (cellular_automata_data$Dist_ID_To_Coral_Patch[i] <= 6))
    {
      #current cell is eligible for damselfish
      #
      #
      #randomly sample from damselfish density sampling distribution that corresponds
      #to the distance id from coral patch
      #
      #obtain path name of damsel density sampling data file that corresponds to distance to coral patch id
      #
      damsel_density_file_path<-damsel_density_sampling_file_paths_data$damsel_density_file_path_name[cellular_automata_data$Dist_ID_To_Coral_Patch[i]]
      damsel_density_sample_data<-read.csv(damsel_density_file_path, stringsAsFactors = FALSE, header = T)
      #
      #read damsel density sample data into a vector and sample to determine if a damselfish will occupy
      #the current cell
      #
      damsel_density_data_vec<-damsel_density_sample_data[,2]
      #
      damselfish_found<-sample(damsel_density_data_vec,1)
      #
      #if a damselfish found for current cell, now need to determine the size of the damselfish according
      #to agonism exhibited by the damselfish
      if (damselfish_found == 1)
      {
        #
        #Zero out urchin density for a damselfish protected lawn
        cellular_automata_data$Urchin_Density[i]<-0.0
        #a damselfish has been found for current cell
        #find agonism severity of damselfish according to appropriate damselfish responsiveness sample data
        #and set the damselfish size for the current cell
        #
        damsel_response_file_path<-damsel_responses_sampling_file_paths_data$damsel_resp_file_path_name[cellular_automata_data$Dist_ID_To_Coral_Patch[i]]
        damsel_response_sample_data<-read.csv(damsel_response_file_path, stringsAsFactors = FALSE, header = T)
        #
        #read damsel responsiveness sample data into a vector and sample to determine the degree of agonism
        #if higher degree of agonism will select largest fish size available for distance id
        #
        damsel_response_data_vec<-damsel_response_sample_data[,2]
        #
        damselfish_degree_agonism<-sample(damsel_response_data_vec,1)
        #
        if (damselfish_degree_agonism == 1)
        {
          #most agonistic response; select largest damselfish for current distance id
          damselfish_size<-damsel_size_by_dist_selection_data$High_Resp_Length[cellular_automata_data$Dist_ID_To_Coral_Patch[i]]
        }#end damselfish degree of agonism is max
        else{
          #least agonistic response;select shortest damselfish for current distance id
          damselfish_size<-damsel_size_by_dist_selection_data$Low_Resp_Length[cellular_automata_data$Dist_ID_To_Coral_Patch[i]]
        }
        #set damselfish length for current cell
        cellular_automata_data$DamselfishLength_cm[i]<-damselfish_size
      }#end if damselfish found for current cell
      #
    }#end if distance id between 1 and 6, inclusive
    #
    #Calculate initial carbonate per cell and per reef zone based on substrate cover % and mean annual 
    #carbonate growth
    #
    #Carbonate production contribution only from live coral cover and coralline algae substrate
    #
    carbonate_production_curr_cell<-0.0
    biomass_coralline_algae<-0.0
    biomass_live_coral<-0.0
    if (cellular_automata_data$DeadCoralPercent[i] > 0.0)
    {
      #dead coral encrusted with coralline algae is present in cell
      #
      #calculate square cm of corraline algal cover
      #
      #print("coralline algae found")
      #print("current cell id")
      #print(i)
      #
      #compute grams of coralline algae produced in 1 yr as product of annual 
      #carbonate production adjusted by % cover and by areal adjustment factor
      #
      areal_adj_factor<-4/((1.0/cellular_automata_data$ChainLengthX[i] + 1.0/cellular_automata_data$ChainLengthY[i])**2)
      #print("areal adj factor")
      #print(areal_adj_factor)
      #print("percent cover")
      #print(cellular_automata_data$DeadCoralPercent[i])
      biomass_coralline_algae<-(coralline_algae_annual_growth_rate*1000.0)*(cellular_automata_data$DeadCoralPercent[i]/100.0)*areal_adj_factor
      #
      #compute biomass (g) of corraline algae by multiplying planar area adjusted by areal
      #factor by density
      #
      #print("initial corraline algae biomass")
      #print(biomass_coralline_algae)
      carbonate_production_curr_cell<-carbonate_production_curr_cell + biomass_coralline_algae
    }#end if current cell contains corraline algae
    if (cellular_automata_data$LiveCoralPercent[i] > 0.0)
    {
      #live coral is present in cell
      #
      #print("live coral found")
      #print("current cell id")
      #print(i)
      #
      #compute grams of live coral produced in 1 yr as product of annual 
      #carbonate production adjusted by % cover
      #
      #print("percent cover")
      #print(cellular_automata_data$LiveCoralPercent[i])
      biomass_live_coral<-(live_coral_annual_carbon_prod_kg_yr_sq_meter*1000.0)*(cellular_automata_data$LiveCoralPercent[i]/100.0)
      #
      #
      #print("initial live coral biomass")
      #print(biomass_live_coral)
      carbonate_production_curr_cell<-carbonate_production_curr_cell + biomass_live_coral
    }#end if current cell contains live coral
    #
    #Calculate net carbonate for current cell
    #
    cellular_automata_data$Initial_Carbonate_g[i]<-carbonate_production_curr_cell
    #
    #update count of damselfish to produce damselfish density
    if (cellular_automata_data$Reef_Zone[i] == 1)
    {
      if (cellular_automata_data$DamselfishLength_cm[i] > 0.0)
      {
        damselfish_count_zone_1<-damselfish_count_zone_1 + 1
        damselfish_length_zone_1<-damselfish_length_zone_1 + cellular_automata_data$DamselfishLength_cm[i]
      }
    }
    if (cellular_automata_data$Reef_Zone[i] == 2)
    {
      if (cellular_automata_data$DamselfishLength_cm[i] > 0.0)
      {
        damselfish_count_zone_2<-damselfish_count_zone_2 + 1
        damselfish_length_zone_2<-damselfish_length_zone_2 + cellular_automata_data$DamselfishLength_cm[i]
      }
    }
    #
  }#end loop on cells to initialize urchin distributions per reef zone type
  #
  #compute mean damselfish length and density per zone
  #
  cellular_automata_simulation_data$mean_damselfish_density_zone1[curr_simulation_run]<-(damselfish_count_zone_1/cellular_automata_num_cells)*2.0
  cellular_automata_simulation_data$mean_damselfish_density_zone2[curr_simulation_run]<-(damselfish_count_zone_2/cellular_automata_num_cells)*2.0
  cellular_automata_simulation_data$mean_damselfish_size_zone1[curr_simulation_run]<-(damselfish_length_zone_1/damselfish_count_zone_1)
  cellular_automata_simulation_data$mean_damselfish_size_zone2[curr_simulation_run]<-(damselfish_length_zone_2/damselfish_count_zone_2)
  #
output_data_ctr<-0
for (curr_day in 1:total_num_days)
{
  #
  #loop over total number of cells
  #
  #
  zone1_total_urchins<-0.0
  zone2_total_urchins<-0.0
  zone1_total_damselfish<-0
  zone2_total_damselfish<-0
  zone1_total_damselfish_size<-0.0
  zone2_total_damselfish_size<-0.0
  coral_cover_total_urchins<-0.0
  zone1_total_urchins_exclude_patches_damselfish<-0.0
  zone2_total_urchins_exclude_patches_damselfish<-0.0
  zone1_total_quadrats_with_damselfish<-0
  zone2_total_quadrats_with_damselfish<-0
  zone1_total_quadrats_coral_patch<-0
  zone2_total_quadrats_coral_patch<-0
  zone1_total_quadrats<-0
  zone2_total_quadrats<-0
  #
  for (curr_cell in 1:cellular_automata_num_cells)
  {
   #
    output_data_ctr<-output_data_ctr+1
    urchin_ejection_ctr<-0
    #Calculate the total biomass in grams that can be consumed by urchins in current cell
   #
    #print("current day")
    #print(curr_day)
    #print("current cell")
    #print(curr_cell)
   total_biomass_consumption_urchins<-cellular_automata_data$Urchin_Density[curr_cell]*urchin_daily_bioerosion_rate
   #
   #calculate total biomass available for consumption per substrate type
   #
   #calculate biomass for coralline algae
   #
   #
    live_coral_cover_sq_cm<-0.0
    corraline_cover_sq_cm<-0.0
    biomass_live_coral<-0.0
    biomass_coralline_algae<-0.0
    final_biomass_corraline_algae<-0.0
    biomass_unprotected_algal_lawn<-0.0
    biomass_protected_algal_lawn<-0.0
    live_coral_found_factor<-0
    coralline_algae_found_factor<-0
    protected_algal_lawn_found_factor<-0
    unprotected_algal_lawn_found_factor<-0
    total_net_biomass_current_cell<-0.0
    #
    #
    if (cellular_automata_data$DeadCoralPercent[curr_cell] > 0.0)
    {
      #dead coral encrusted with coralline algae is present in cell
      #
      coralline_algae_found_factor<-1
      #
      #print("dead coral in current cell")
      #print("current cell id")
      #print(curr_cell)
      #
      #Calculate increased biomass in grams for coralline algae
      #
      areal_adj_factor<-4/((1.0/cellular_automata_data$ChainLengthX[curr_cell] + 1.0/cellular_automata_data$ChainLengthY[curr_cell])**2)
      biomass_coralline_algae<-((coralline_algae_annual_growth_rate*1000.0)/time_period_intervals)*(cellular_automata_data$DeadCoralPercent[curr_cell]/100.0)*areal_adj_factor
      #print("coralline algae biomass")
      #print(biomass_coralline_algae)
      #
      #
    }#end if corraline algae present in current cell
    #
    #Determine if live coral is located with current cell
    if (cellular_automata_data$LiveCoralPercent[curr_cell] > 0.0)
    {
      #Live coral is present in cell
      #
      live_coral_found_factor<-1
      #calculate square cm of live coral cover
      #
      #print("live coral in current cell")
      #print("current cell id")
      #print(curr_cell)
      #
      #compute biomass of live coral cover in grams
      #
      biomass_live_coral<-((live_coral_annual_carbon_prod_kg_yr_sq_meter*1000.0)/time_period_intervals)*(cellular_automata_data$LiveCoralPercent[curr_cell]/100.0)
      #print("live coral biomass")
      #print(biomass_live_coral)
      #
      #Increase live coral cover by cover growth rate
      #
      live_coral_percent_current<-cellular_automata_data$LiveCoralPercent[curr_cell]
      cellular_automata_data$LiveCoralPercent[curr_cell]<-(1.0 + (live_coral_cover_growth_rate/time_period_intervals))*cellular_automata_data$LiveCoralPercent[curr_cell]
      live_coral_percent_next<-cellular_automata_data$LiveCoralPercent[curr_cell]
      live_coral_percent_difference<-live_coral_percent_next - live_coral_percent_current
      #
      #
      #Proportionally reduce cover of other substrate types
      #
      #compute total cover % excluding live coral
      #
      total_non_live_coral_percent<-cellular_automata_data$DeadCoralPercent[curr_cell] + cellular_automata_data$NoNutritionalValuePercent[curr_cell] + cellular_automata_data$UnprotectedAlgalLawnPercent[curr_cell] + cellular_automata_data$ProtectedAlgalLawnPercent[curr_cell]
      #
      #
      #Update % cover for each substrate type that is not live coral if any present
      #
      if (total_non_live_coral_percent > 0.0)
      {
        #print("non-live coral % is non-zero")
        #print("dead coral % before")
        #print(cellular_automata_data$DeadCoralPercent[curr_cell])
        #print("no nutritional % before")
        #print(cellular_automata_data$NoNutritionalValuePercent[curr_cell])
        #print("unprotected algae % before")
        #print(cellular_automata_data$UnprotectedAlgalPercent[curr_cell])
        #print("protected algae % before")
        #print(cellular_automata_data$ProtectedAlgalPercent[curr_cell])
        cellular_automata_data$DeadCoralPercent[curr_cell]<-cellular_automata_data$DeadCoralPercent[curr_cell] - (cellular_automata_data$DeadCoralPercent[curr_cell]/total_non_live_coral_percent)*live_coral_percent_difference
        cellular_automata_data$NoNutritionalValuePercent[curr_cell]<-cellular_automata_data$NoNutritionalValuePercent[curr_cell] - (cellular_automata_data$NoNutritionalValuePercent[curr_cell]/total_non_live_coral_percent)*live_coral_percent_difference
        cellular_automata_data$UnprotectedAlgalPercent[curr_cell]<-cellular_automata_data$UnprotectedAlgalLawnPercent[curr_cell] - (cellular_automata_data$UnprotectedAlgalLawnPercent[curr_cell]/total_non_live_coral_percent)*live_coral_percent_difference
        cellular_automata_data$ProtectedAlgalPercent[curr_cell]<-cellular_automata_data$ProtectedAlgalLawnPercent[curr_cell] - (cellular_automata_data$ProtectedAlgalLawnPercent[curr_cell]/total_non_live_coral_percent)*live_coral_percent_difference
        #print("dead coral % after")
        #print(cellular_automata_data$DeadCoralPercent[curr_cell])
        #print("no nutritional % after")
        #print(cellular_automata_data$NoNutritionalValuePercent[curr_cell])
        #print("unprotected algae % after")
        #print(cellular_automata_data$UnprotectedAlgalPercent[curr_cell])
        #print("protected algae % after")
        #print(cellular_automata_data$ProtectedAlgalPercent[curr_cell])
      } 
      #
    }#end if live coral present in current cell
    #
   #compute net biomass for cell 
   #print("cell id")
   #print(curr_cell)
   #print("total net biomass before for current cell")
   #print(cellular_automata_data$Initial_Carbonate_g[curr_cell])
   #print("biomass live coral")
   #print(biomass_live_coral)
   #print("biomass coralline algae")
   #print(biomass_coralline_algae)
   cellular_automata_data$Initial_Carbonate_g[curr_cell]<-cellular_automata_data$Initial_Carbonate_g[curr_cell] + biomass_live_coral + biomass_coralline_algae - total_biomass_consumption_urchins
   #print("urchin density")
   #print(cellular_automata_data$Urchin_Density[curr_cell])
   #print("total bioerosion by urchins")
   #print(total_biomass_consumption_urchins)
   #print("total net biomass after for current cell")
   #print(cellular_automata_data$Initial_Carbonate_g[curr_cell])
   #
   #
   #
   #
   #
   join_aggregation_count<-0
   total_aggregating_urchins_curr_cell<-0
   urchin_density_curr_cell<-round(cellular_automata_data$Urchin_Density[curr_cell])
   #print("urchin density current cell")
   #print(urchin_density_curr_cell)
   #determine if aggregations will form in current cell
   #
   if (urchin_density_curr_cell >= min_urchin_aggreg_density)
   {
     #
     #minimum urchin density for aggregation exists
     #
     #
     #randomly generate (x,y) coordinates within 100 cm x 100 cm quadrat for as many urchins
     #in current cell 
     #use a 1 cm resolution for coordinates
     #
     x<-1:100
     y<-1:100
     x_coord_vec<-rep(0,urchin_density_curr_cell)
     y_coord_vec<-rep(0,urchin_density_curr_cell)
     #
     for (j in 1:urchin_density_curr_cell)
     {
      #
      #randomly select both x and y coordinates 
      #
       x_coord_vec_indx<-sample(x,1)
       y_coord_vec_indx<-sample(y,1)
       x_coord_vec[j]<- x_coord_vec_indx
       y_coord_vec[j]<- y_coord_vec_indx
      # 
     }
     #print("x comp vector for random points")
     #print(x_coord_vec)
     #print("y comp vector for random points")
     #print(y_coord_vec)
     #
     #perform k-means clustering using k = 1,....,m for random assortment of urchins in 100 cm x 100 cm quadrat
     #will stop if aggregation is found or if largest cluster is smaller than minimum aggregation size 
     #
     urchins_aggregating<-FALSE
     cluster_too_small<-FALSE
     m<-1
     while (urchins_aggregating == FALSE | cluster_too_small == FALSE)
     {
      random_urchin_centroids_data<-data.frame(x_coord_vec,y_coord_vec)
      kmeans_results<-kmeans(random_urchin_centroids_data,m)
      #print("kmeans results for random points")
      #print(kmeans_results)
     #
     #
     #Determine if of m clusters, there is one that contains at least as many urchins
     #as the minimum number to form an aggregation
     #
     num_clusters_too_small<-0
     for (w in 1:m)
     {
       if (kmeans_results$size[w] < min_urchin_aggreg_density)
       {
         num_clusters_too_small<-num_clusters_too_small + 1
       }
     }
     if (num_clusters_too_small == m)
     {
       #
       #all clusters smaller than minimum size; break out of while loop
       #
       cluster_too_small<-TRUE
       break;
       #
     }
     #
     #here if at least one cluster had a sufficient number of urchins
     #
     #print("at least one cluster has sufficient number of urchins")
     #
     #now for largest cluster, determine if aggregation exists
     #
     max_cluster_size<-max(kmeans_results$size)
     #print("max cluster size")
     #print(max_cluster_size)
     #
     #loop over clusters to find index of cluster with max number of points
     for (p in 1:m)
     {
       if (kmeans_results$size[p] == max_cluster_size)
       {
         #
         #found cluster with largest number of points
         #
         max_cluster_index<-p;
         break;
         #
       }
     }
     #
     #obtain within cluster sum of squares for largest cluster
     #
     #print("cluster with max points index")
     #print(max_cluster_index)
      random_sum_of_squares<-kmeans_results$withinss[max_cluster_index]
     #
     #lookup sum of squares for rectangular packing using urchin density as index
     #
      if (urchin_density_curr_cell<=total_sum_squares_data_size)
      {
        circle_packing_sum_squares<-circle_packing_sum_squares_data$WCSS[max_cluster_size]
      }
      else
      {
        circle_packing_sum_squares<-circle_packing_sum_squares_data$WCSS[total_sum_squares_data_size]
      }
     #
     #compute ratio of total sum square of random assortment of urchin centroids to 
     #total sum square of optimal rectangular circle packing
     #print("random sum of squares")
     #print(random_sum_of_squares)
     #print("circle packing sum of squares")
     #print(circle_packing_sum_squares)
      urchin_packing_ratio<-circle_packing_sum_squares/random_sum_of_squares
     #print("urchin packing ratio")
     #print(urchin_packing_ratio)
     #
     #determine if urchin packimg ratio lies within aggregation affinity parameter
     #
      if (urchin_packing_ratio >= cellular_automata_params$Aggregation_Affinity_Param[1])
      {
        #urchins in one of the clusters are aggregating!
        #
        urchins_aggregating<-TRUE
        #print("urchin in one of clusters is aggregating")
        #print("max cluster index")
        #print(max_cluster_index)
        break;
        #
      }
     m<-m + 1 
     }#end while loop
     #
     #if urchins are aggregating in largest cluster, then randomly select k points in aggregating cluster
     #and loop through urchins in non-aggregating clusters, and each time compute distance between 
     #each non-aggregating urchin and k points randomly selected from aggregating cluster. If we find
     #a distance <= 10 cm, we can apply 6/9 probability of urchin joining the aggregation
     #
     #
     if (urchins_aggregating == TRUE)
     {
       #
       #produce index vector that will contain all indeces of cluster that contains aggregating urchins
       #
       #print("URCHINS ARE AGGREGATING IN CURRENT CELL")
       aggregate_indx_vec<-rep(0,kmeans_results$size[max_cluster_index])
       #print("init aggregation index vector")
       #print(aggregate_indx_vec)
       cntr<-1
       for (indx in 1:length(kmeans_results$cluster))
       {
         if (kmeans_results$cluster[indx] == max_cluster_index)
         {
           #
           #current entry in cluster data has index matching the aggregating cluster index
           aggregate_indx_vec[cntr]<-indx
           cntr<-cntr+1
           #
         }
       }#end for loop on iterating through cluster vector
       #
       #Randomly select indeces for points that lie within the aggregation
       #
       for (f in 1:num_random_points)
       {
         random_points_index_vec[f]<-sample(aggregate_indx_vec,1)
       }
       #print("random points")
       #print(num_random_points)
       #print("random points index vec")
       #print(random_points_index_vec)
       #
       #Find distance between random points and each point in quadrat that is not contained in cluster
       #with aggregating urchins
       #
       #Loop over all points over all clusters, and if a point is not contained in aggregating cluster,
       #compute distance from that point to randomly selected points from aggregating cluster
       #if distance is within maximum distance for urchins to join aggregations, apply probability
       #of joining aggregation 
       #
       #
       #
       for (g in 1:length(kmeans_results$cluster))
       {
         #
         #determine if current point is in aggregating cluster
         #
         if (kmeans_results$cluster[g] != max_cluster_index)
         {
           #
           #Found a point that does not belong to aggregating cluster
           #Compute distance between this point and random set of points
           #obtained from aggregating cluster. Will break out of loop
           #if we find distance between points is within threshold for 
           #joining aggregation
           #
           dist_to_aggregation_vec<-numeric(num_random_points)
           for (w in 1:num_random_points)
           {
             #
             #compute distance between current non-aggregate point and each randomly
             #selected point from aggregation
             #
             x1_coord<-random_urchin_centroids_data$x_coord_vec[g]
             x0_coord<-random_urchin_centroids_data$x_coord_vec[random_points_index_vec[w]]
             y1_coord<-random_urchin_centroids_data$y_coord_vec[g]
             y0_coord<-random_urchin_centroids_data$y_coord_vec[random_points_index_vec[w]]
             dist_to_aggregation_vec[w]<-sqrt((y1_coord-y0_coord)^2 + (x1_coord-x0_coord)^2)
             #
           }#end for loop 
           #
           #check distance to aggregation vector for current non-aggregation point
           #to determine if within maximum distance for joining aggregation
           #
           #print("distance vector to aggregation")
           #print(dist_to_aggregation_vec)
           if (min(dist_to_aggregation_vec) < cellular_automata_params$Max_Aggregation_Join_Dist_cm[1])
           {
             #found a point that could potentially join the aggregation
             #print("found point that could potentially join aggregation")
             #print("vector of distances for 10 randomly selected points in aggregation")
             #print(dist_to_aggregation_vec)
             #print("minimum distance found")
             #print(min(dist_to_aggregation_vec))
             #
             #randomly sample probability vector for joining aggregation
             #
             join_aggregation_indicator<-sample(join_aggregation_freq_data_vec,1)
             #
             if (join_aggregation_indicator == 1)
             {
               #
               #urchin will join aggregation
               #update count of urchins that will join aggregation
               #
               join_aggregation_count<-join_aggregation_count + 1 
               #
             }#end check if urchins join aggregation
           }#end if urchins within max aggregation join distance
           #
         }# end if found point that does not belong to aggregating cluster
         #
       }#end loop on all points across all clusters
       #
       #Update total count of urchins aggregating in current cell
       #
       total_aggregating_urchins_curr_cell<- join_aggregation_count + kmeans_results$size[max_cluster_index]
       #
       #
     }#end if urchins are aggregating
     #
   }# end if minimum threshold of urchins for aggregation exists
   #
   #Determine if non-aggregating urchins are homing
   #print("total aggregating urchins current cell")
   #print(total_aggregating_urchins_curr_cell)
   total_non_aggregating_urchins_curr_cell<-urchin_density_curr_cell - total_aggregating_urchins_curr_cell
   #
   total_homing_urchins_current_cell<-0
   if (total_non_aggregating_urchins_curr_cell > 0)
   {
     #
     #compute probability that urchins may be homing as a function of urchin density in cell
     #
     #will use linear function that estimates percentage of crevice fidelity as a function
     #of urchin density (Carpenter 1984)
     #     
     homing_behavior_frequency<-(-2.0)*total_non_aggregating_urchins_curr_cell + 64.0
     if ((homing_behavior_frequency < 0.0) | (cellular_automata_data$DamselfishLength_cm[curr_cell] > 0.0))
     {
      #if a cell contains a damselfish, will not permit homing behavior
      #
      homing_behavior_frequency<-0.0
     }# end if homing behavior frequency less than zero
     #
     if (homing_behavior_frequency > 0.0)
     {
       #There is a non-zero probability that urchins
       #in current cell may exhibit homing behavior
       #
       #
       #create a sampling distribution using bootstrapping
       #create a vector initialized with all zero's
       #
       homing_sampling_vec<-numeric(100)
       #
       #set number of entries in sampling vector to "1" to represent homing
       #according to the homing probability
       #
       for (j in 1:homing_behavior_frequency)
       {
        homing_sampling_vec[j]<-1 
       }
       #
       #loop over all non-aggregating urchins to determine how many are homing
       #
       for (w in 1:total_non_aggregating_urchins_curr_cell)
       {
         #
         #randomly sample from homing sampling vector to determine if urchin is homing
         #
         urchin_homing_indicator<-sample(homing_sampling_vec,1)
         if (urchin_homing_indicator == 1)
         {
           #current urchin is homing
           #
           total_homing_urchins_current_cell<-total_homing_urchins_current_cell + 1
           #
         }#end if current urchin is homing
         #   
       }# end loop to find how many non-aggregating urchins are homing
       #print("total homing urchins for current cell")
       #print(total_homing_urchins_current_cell)
       #
     }# end if homing behavior possible for current cell
     
   }# end if there are non-zero number of non-aggregating urchins
   #
   #Know number of urchins that are aggregating, homing, and not homing for current cell
   #
   #prepare for transition to next cell
   #
   #Homing urchins will remain in current cell; non-homing urchins may transition to other cells
   #
   #
   next_cell_urchin_transition_data$current_urchin_density[curr_cell]<-urchin_density_curr_cell
   #
   total_non_homing_urchins<-urchin_density_curr_cell - total_homing_urchins_current_cell - total_aggregating_urchins_curr_cell
   #
   #loop over total non-homing urchins if non-zero to determine their next-cell transition
   #
   if (total_non_homing_urchins > 0)
   {
     #non-homing, non-aggregating urchins found
     #
     #randomly sample next-cell frequency data; determine which frequency data
     #to sample according to the number of urchins in current cell
     #
     next_cell_trans_sample_vec<-numeric(total_non_homing_urchins)
     if (urchin_density_curr_cell > cellular_automata_params$Urchin_Density_Sampling_Vec_Threshold[1])
     {
       #
       #use high density sampling vector
       #
       next_cell_trans_sample_vec<-sample(high_density_urchin_dist_sampling_vec,total_non_homing_urchins)
      }
    else
    {
     #use moderate density sampling vector
     #
      next_cell_trans_sample_vec<-sample(mod_density_urchin_dist_sampling_vec,total_non_homing_urchins)
     }
     #
     #Create vector that will contain the physical cell id's that urchins will transition to
     #
     next_cell_id_transition_vec<-numeric(total_non_homing_urchins)  
     for (k in 1:total_non_homing_urchins)
     {
       #
       #loop over non-homing urchins and determine next-cell to transition too
       #will check for damselfish avoidance behavior for each next-cell 
       if (next_cell_trans_sample_vec[k] == "H")
       {
         #
         #urchin is to remain in current cell
         next_cell_id_transition_vec[k]<-curr_cell
         #
       }
       else
       { 
         #call next cell transition function to determine next cell for transition
         #print("current cell")
         #print(curr_cell)
         #print("next cell")
         #print(next_cell_trans_sample_vec[k])
         #print("call next cell transition function")
         return_list<-NextCellTransition(curr_cell,next_cell_trans_sample_vec[k],urchin_density_curr_cell)
         next_cell_id_transition_vec[k]<-return_list$next_cell_id
         if (return_list$damselfish_eject_flag)
         {
           urchin_ejection_ctr<-urchin_ejection_ctr + 1
         }
       }#end next cell is not home cell
       #
     }#end loop on non-homing urchins
      #
      #loop over next-cell transition vector for non-homing urchins and assign 
      #urchin counts
      #
      #if next-cell does not match current cell id this is an output urchin
      #if next-cell matches current cell id this is an input urchin 
      #
    #print("total non-homing urchins")
    #print(total_non_homing_urchins)
      for (w in 1:total_non_homing_urchins)
      {
        #
        if (next_cell_id_transition_vec[w] == 0)
        {
          #non-homing urchin will leave lattice
          #
          next_cell_urchin_transition_data$number_output_urchins[curr_cell]<-1 + next_cell_urchin_transition_data$number_output_urchins[curr_cell]
        }
        else
        {
         if (next_cell_id_transition_vec[w] != curr_cell)
         {
          #a non-homing urchin will leave home cell
          #print("non homing urchin leaving cell")
          #print("cell id")
          #print(curr_cell)
          #print("next cell transition data - number output urchins")
          #print(next_cell_urchin_transition_data$number_output_urchins[curr_cell])
          next_cell_urchin_transition_data$number_output_urchins[curr_cell]<-1 + next_cell_urchin_transition_data$number_output_urchins[curr_cell]
          next_cell_urchin_transition_data$number_input_urchins[next_cell_id_transition_vec[w]]<-1 + next_cell_urchin_transition_data$number_input_urchins[next_cell_id_transition_vec[w]]
          }#end if next cell is not home cell
        }#end else if next cell is not zero   
      }#end loop on non-homing urchins
    }#end if non-zero number of non-homing urchins
   #
   #Check if any aggregating urchins
   #
   if (total_aggregating_urchins_curr_cell > 0.0)
   {
     #there is an aggregation of urchins
     #
     #
     #sample from disaggregation subgroup vector
     #to determine how many subgroups to form
     #
     num_aggregation_subgroups<-sample(num_aggregation_subgroups_sample_data_vec,1)
     #
     #calculate proportion of disaggregation
     #
     proportion_disaggregated<-1.0/num_aggregation_subgroups
     #
     #determine whether aggregation should be disaggregated at this time step
     disaggregation_indicator<-sample(disaggregation_frequency_data_vec,1)
     if (disaggregation_indicator == 0)
     {
       #aggregation will remain as one complete group
       proportion_disaggregated<-1.0
       num_aggregation_subgroups<-1
       #print("NO DISAGGREGATION ALLOWED AT THIS TIME")
     }
     if(proportion_disaggregated > 0.0)
     {
      # 
      #
      #
      #
      #loop over aggregation subgroups and determine next-cell transition for each subgroup
      #
      total_urchin_subgroup_ctr<-0
      for (z in 1:num_aggregation_subgroups)
      {
       #
        if (z == num_aggregation_subgroups)
        {
         total_num_urchins_curr_subgroup<-round(total_aggregating_urchins_curr_cell - total_urchin_subgroup_ctr)
        }
        else
        {
          total_num_urchins_curr_subgroup<-round(proportion_disaggregated*total_aggregating_urchins_curr_cell)
          total_urchin_subgroup_ctr<-total_urchin_subgroup_ctr + total_num_urchins_curr_subgroup
        }
       #determine next-cell transition for aggregating urchin subgroup
       next_cell_trans_logical_cell<-sample(high_density_urchin_dist_sampling_vec,1)
       #
       #convert logical cell to physical cell id 
       #
       if (next_cell_trans_logical_cell == "H")
       {
         #aggregating urchins stay at current cell
         next_cell_trans_cell_id<-curr_cell
         #print("aggregating urchins staying at current cell")
       }
       else
       {
         #next-cell for transition is not current cell
         #
         #Call function to determine next cell for transition
         #print("current cell")
         #print(curr_cell)
         #print("next cell")
         #print(next_cell_trans_logical_cell)
         #print("call next cell transition function")
         #print("next cell for trans is not current cell")
         return_list<-NextCellTransition(curr_cell,next_cell_trans_logical_cell,total_num_urchins_curr_subgroup)
         next_cell_trans_cell_id<-return_list$next_cell_id
         if (return_list$damselfish_eject_flag)
         {
           urchin_ejection_ctr<-urchin_ejection_ctr + total_num_urchins_curr_subgroup
           #print("total num urchins ejected")
           #print(urchin_ejection_ctr)
         }
         if (next_cell_trans_cell_id == 0)
         {
           #aggregating urchins will leave lattice
           next_cell_urchin_transition_data$number_output_urchins[curr_cell]<- total_num_urchins_curr_subgroup + next_cell_urchin_transition_data$number_output_urchins[curr_cell]
         }
         else
         {
          if (next_cell_trans_cell_id != curr_cell)
          {
           #all aggregating urchins will leave home cell
           next_cell_urchin_transition_data$number_output_urchins[curr_cell]<- total_num_urchins_curr_subgroup + next_cell_urchin_transition_data$number_output_urchins[curr_cell]
           next_cell_urchin_transition_data$number_input_urchins[next_cell_trans_cell_id]<- total_num_urchins_curr_subgroup + next_cell_urchin_transition_data$number_input_urchins[next_cell_trans_cell_id]
          }#end if aggregating urchin will leave home cell
         }#end else statement that next cell is zero   
       }#end if next cell is not home cell
      }#end loop on aggregation subgroups
     }#end if proportion aggregation is non-zero
    }#end if there are aggregating urchins
   #
   #
   #
   #
   #
   #
   #set output data fields
   #
   cellular_automata_output_data$cell_id[output_data_ctr]<-curr_cell
   cellular_automata_output_data$time_id[output_data_ctr]<-curr_day
   cellular_automata_output_data$zone_id[output_data_ctr]<-cellular_automata_data$Reef_Zone[curr_cell]
   cellular_automata_output_data$urchin_density_before[output_data_ctr]<-urchin_density_curr_cell
   cellular_automata_output_data$net_carbonate[output_data_ctr]<-cellular_automata_data$Initial_Carbonate_g[curr_cell]
   cellular_automata_output_data$damselfish_size[output_data_ctr]<-cellular_automata_data$DamselfishLength_cm[curr_cell]
   cellular_automata_output_data$aggregation_size[output_data_ctr]<-total_aggregating_urchins_curr_cell
   cellular_automata_output_data$num_urchin_ejections[output_data_ctr]<-urchin_ejection_ctr
   #
   if (cellular_automata_data$Reef_Zone[curr_cell] == 1)
   {
     zone1_total_urchins<-zone1_total_urchins + urchin_density_curr_cell
     zone1_total_quadrats<-zone1_total_quadrats + 1
     if (cellular_automata_data$DamselfishLength_cm[curr_cell] > 0)
     {
       zone1_total_damselfish<-zone1_total_damselfish+1
       zone1_total_damselfish_size<- zone1_total_damselfish_size + cellular_automata_data$DamselfishLength_cm[curr_cell]
       zone1_total_quadrats_with_damselfish<-zone1_total_quadrats_with_damselfish + 1
     }
   }
   if (cellular_automata_data$Reef_Zone[curr_cell] == 2)
   {
     zone2_total_urchins<-zone2_total_urchins + urchin_density_curr_cell
     zone2_total_quadrats<-zone2_total_quadrats + 1
     if (cellular_automata_data$DamselfishLength_cm[curr_cell] > 0)
     {
       zone2_total_damselfish<-zone2_total_damselfish+1
       zone2_total_damselfish_size<- zone2_total_damselfish_size + cellular_automata_data$DamselfishLength_cm[curr_cell]
       zone2_total_quadrats_with_damselfish<-zone2_total_quadrats_with_damselfish + 1
     }
   }
   if (cellular_automata_data$Coral_Patch[curr_cell] == 1)
   {
     coral_cover_total_urchins<-coral_cover_total_urchins + urchin_density_curr_cell
     if (cellular_automata_data$Reef_Zone[curr_cell] == 1)
     {
       zone1_total_quadrats_coral_patch<- zone1_total_quadrats_coral_patch + 1
     }
     if (cellular_automata_data$Reef_Zone[curr_cell] == 2)
     {
       zone2_total_quadrats_coral_patch<- zone2_total_quadrats_coral_patch + 1
     }
   }
   else
   {
    if ((cellular_automata_data$Reef_Zone[curr_cell] == 1) & (cellular_automata_data$DamselfishLength_cm[curr_cell] <= 0.0))
     {
      zone1_total_urchins_exclude_patches_damselfish<-zone1_total_urchins_exclude_patches_damselfish + urchin_density_curr_cell
     }
    if ((cellular_automata_data$Reef_Zone[curr_cell] == 2) & (cellular_automata_data$DamselfishLength_cm[curr_cell] <= 0.0))
     {
       zone2_total_urchins_exclude_patches_damselfish<-zone2_total_urchins_exclude_patches_damselfish + urchin_density_curr_cell
     }
   }
   #
  }#end loop on cells
  #
  #initialize next cell transition data per cell for next time period
  #
  for (k in 1:cellular_automata_num_cells)
  {
    #update urchin counts in cellular automata data 
    cellular_automata_data$Urchin_Density[k]<-round(cellular_automata_data$Urchin_Density[k]) + next_cell_urchin_transition_data$number_input_urchins[k] - next_cell_urchin_transition_data$number_output_urchins[k]
    next_cell_urchin_transition_data$current_urchin_density[k]<-0.0
    next_cell_urchin_transition_data$number_input_urchins[k]<-0.0
    next_cell_urchin_transition_data$number_output_urchins[k]<-0.0
  }
  #next_cell_urchin_transition_data<-data.frame(cell_id,current_urchin_density,number_input_urchins,number_output_urchins)
  #next_cell_urchin_transition_data$cell_id<-c(1:cellular_automata_num_cells)
  #print("current day")
  #print(curr_day)
  #
  #find the mean difference in urchin density for each time period
  #will output to file, per time period: 
  #mean density zone 1
  #mean density zone 2
  #difference of mean density
  print("current day")
  print(curr_day)
  print("zone1 total quadrats")
  print(zone1_total_quadrats)
  print("zone2 total quadrats")
  print(zone2_total_quadrats)
  print("zone1 coral patch quadrats")
  print(zone1_total_quadrats_coral_patch)
  print("zone2 coral patch quadrats")
  print(zone2_total_quadrats_coral_patch)
  print("zone1 damselfish quadrats")
  print(zone1_total_quadrats_with_damselfish)
  print("zone2 damselfish quadrats")
  print(zone2_total_quadrats_with_damselfish)
  print("zone1 total urchins")
  print(zone1_total_urchins)
  print("zone2 total urchins")
  print(zone2_total_urchins)
  print("coral patch total urchins")
  print(coral_cover_total_urchins)
  print("zone 1 total urchins exclude coral & damselfish")
  print(zone1_total_urchins_exclude_patches_damselfish)
  print("zone 2 total urchins exclude coral & damselfish")
  print(zone2_total_urchins_exclude_patches_damselfish)
#
  cellular_automata_summary_data$zone1_mean_urchin_density[curr_day]<-2*zone1_total_urchins/(cellular_automata_num_cells)
  cellular_automata_summary_data$zone2_mean_urchin_density[curr_day]<-2*zone2_total_urchins/(cellular_automata_num_cells)
  cellular_automata_summary_data$diff_mean_urchin_density[curr_day]<-cellular_automata_summary_data$zone1_mean_urchin_density[curr_day] - cellular_automata_summary_data$zone2_mean_urchin_density[curr_day]
  cellular_automata_summary_data$zone1_mean_damselfish_density[curr_day]<-2*zone1_total_damselfish/(cellular_automata_num_cells)
  cellular_automata_summary_data$zone2_mean_damselfish_density[curr_day]<-2*zone2_total_damselfish/(cellular_automata_num_cells)
  cellular_automata_summary_data$zone1_mean_damselfish_size[curr_day]<-zone1_total_damselfish_size/zone1_total_damselfish
  cellular_automata_summary_data$zone2_mean_damselfish_size[curr_day]<-zone2_total_damselfish_size/zone2_total_damselfish
  cellular_automata_summary_data$coral_patch_mean_urchin_density[curr_day]<-coral_cover_total_urchins/cellular_automata_params$Coral_Patch_Size_sq_meter[1]
  cellular_automata_summary_data$zone1_exclude_patch_damselfish_urchin_density[curr_day]<-zone1_total_urchins_exclude_patches_damselfish/(zone1_total_quadrats - zone1_total_quadrats_coral_patch - zone1_total_quadrats_with_damselfish)
  cellular_automata_summary_data$zone2_exclude_patch_damselfish_urchin_density[curr_day]<-zone2_total_urchins_exclude_patches_damselfish/(zone2_total_quadrats - zone2_total_quadrats_coral_patch - zone2_total_quadrats_with_damselfish)
  #
  #
  #Update urchin counts per zone to match mean urchin density per zone; compensation for
  #urchins that migrate out of the square lattice
  #
  #obtain mean urchin density for zone 1 and zone 2 determine if each one is within its prescribed tolerance
  #
  #
  #zone 1 processing
  diff_urchin_means_zone1<-reef_zone_data$Urchin_Density_Mean[1] - cellular_automata_summary_data$zone1_mean_urchin_density[curr_day]
  if (diff_urchin_means_zone1 > cellular_automata_params$Urchin_Density_Difference_Tolerance[1])
  {
    #difference in urchin density between optimal level and current level exceeds tolerance
    #calculate number of urchins to supplement
    #
    num_urchins_supplement<-round(diff_urchin_means_zone1*cellular_automata_params$Zone1_Study_Area_Size_sq_meter[1])
    #
    num_urchins_ctr<-num_urchins_supplement
    #randomly sample from urchin density distribution for zone 1 and from border cells for zone 1 and populate
    #cells with urchins until urchin deficit has been fulfilled 
    while (num_urchins_ctr > 0)
    {
      #randomly sample from urchin distribution for zone 1
      num_urchins_curr_cell<-sample(urchin_density_data_zone_1_vec,1)
      #
      #randomly sample cell id from border and near border cells for zone 1
      #
      border_cell_id_for_supplement<-sample(border_cell_id_zone1_vector,1)
      #
      #Increment urchin count for border cell in zone 1 according to sampled density
      cellular_automata_data$Urchin_Density[border_cell_id_for_supplement]<-cellular_automata_data$Urchin_Density[border_cell_id_for_supplement] + num_urchins_curr_cell
      num_urchins_ctr<-num_urchins_ctr-num_urchins_curr_cell
      #check to see if number of urchins for current cell will exceed maximum with addition of urchins
      #
      if (cellular_automata_data$Urchin_Density[border_cell_id_for_supplement] > cellular_automata_params$Max_Number_Urchins_Cell[1])
      {
        #addition of urchins will exceed maximum;back out urchins
        cellular_automata_data$Urchin_Density[border_cell_id_for_supplement]<-cellular_automata_data$Urchin_Density[border_cell_id_for_supplement] - num_urchins_curr_cell
        num_urchins_ctr<-num_urchins_ctr+num_urchins_curr_cell
      }
    }#end while loop on number of urchins to supplement
  }#end if difference in mean urchin density exceeds tolerance  
  #
  #zone 2 processing
  diff_urchin_means_zone2<-reef_zone_data$Urchin_Density_Mean[2] - cellular_automata_summary_data$zone2_mean_urchin_density[curr_day]
  if (diff_urchin_means_zone2 > cellular_automata_params$Urchin_Density_Difference_Tolerance[1])
  {
    #difference in urchin density between optimal level and current level exceeds tolerance
    #calculate number of urchins to supplement
    #
    num_urchins_supplement<-round(diff_urchin_means_zone2*cellular_automata_params$Zone2_Study_Area_Size_sq_meter[1])
    #
    num_urchins_ctr<-num_urchins_supplement
    #randomly sample from urchin density distribution for zone 2 and from border cells for zone 2 and populate
    #cells with urchins until urchin deficit has been fulfilled 
    while (num_urchins_ctr > 0)
    {
      #randomly sample from urchin distribution for zone 2
      num_urchins_curr_cell<-sample(urchin_density_data_zone_2_vec,1)
      #
      #randomly sample cell id from border and near border cells for zone 2
      #
      border_cell_id_for_supplement<-sample(border_cell_id_zone2_vector,1)
      #
      #Increment urchin count for border cell in zone 1 according to sampled density
      cellular_automata_data$Urchin_Density[border_cell_id_for_supplement]<-cellular_automata_data$Urchin_Density[border_cell_id_for_supplement] + num_urchins_curr_cell
      num_urchins_ctr<-num_urchins_ctr-num_urchins_curr_cell
      #check to see if number of urchins for current cell will exceed maximum with addition of urchins
      #
      if (cellular_automata_data$Urchin_Density[border_cell_id_for_supplement] > cellular_automata_params$Max_Number_Urchins_Cell[1])
      {
        #addition of urchins will exceed maximum;back out urchins
        cellular_automata_data$Urchin_Density[border_cell_id_for_supplement]<-cellular_automata_data$Urchin_Density[border_cell_id_for_supplement] - num_urchins_curr_cell
        num_urchins_ctr<-num_urchins_ctr+num_urchins_curr_cell
      }
    }#end while loop on number of urchins to supplement
  }#end if   
  #
}#end time loop
#
output_data_file_path<-""
summary_data_file_path<-""
output_data_file_path<-paste0(cellular_automata_output_data_file_path,"_",(curr_simulation_run+file_id_offset),".csv")
summary_data_file_path<-paste0(cellular_automata_summary_data_file_path,"_",(curr_simulation_run+file_id_offset),".csv")
write.csv(cellular_automata_output_data,output_data_file_path,row.names = FALSE)
write.csv(cellular_automata_summary_data,summary_data_file_path,row.names = FALSE)
print("completed time loop")
#
cellular_automata_simulation_data$simulation_id[curr_simulation_run]<-curr_simulation_run
cellular_automata_simulation_data$number_of_days[curr_simulation_run]<-total_num_days
cellular_automata_simulation_data$aggregation_affinity_param[curr_simulation_run]<-cellular_automata_params$Aggregation_Affinity_Param[1]
cellular_automata_simulation_data$max_size_urchin_aggregate[curr_simulation_run]<-cellular_automata_params$Max_Number_Urchins_Cell[1]
cellular_automata_simulation_data$zone1_mean_urchin_density<-cellular_automata_summary_data$zone1_mean_urchin_density[total_num_days]
cellular_automata_simulation_data$zone2_mean_urchin_density<-cellular_automata_summary_data$zone2_mean_urchin_density[total_num_days]
cellular_automata_simulation_data$coral_patch_init_mean_urchin_density[curr_simulation_run]<-cellular_automata_summary_data$coral_patch_mean_urchin_density[1]
cellular_automata_simulation_data$coral_patch_final_mean_urchin_density[curr_simulation_run]<-cellular_automata_summary_data$coral_patch_mean_urchin_density[total_num_days]
cellular_automata_simulation_data$exclude_patch_damselfish_init_zone1_mean_urchin_density[curr_simulation_run]<-cellular_automata_summary_data$zone1_exclude_patch_damselfish_urchin_density[1]
cellular_automata_simulation_data$exclude_patch_damselfish_final_zone1_mean_urchin_density[curr_simulation_run]<-cellular_automata_summary_data$zone1_exclude_patch_damselfish_urchin_density[total_num_days]
cellular_automata_simulation_data$exclude_patch_damselfish_init_zone2_mean_urchin_density[curr_simulation_run]<-cellular_automata_summary_data$zone2_exclude_patch_damselfish_urchin_density[1]
cellular_automata_simulation_data$exclude_patch_damselfish_final_zone2_mean_urchin_density[curr_simulation_run]<-cellular_automata_summary_data$zone2_exclude_patch_damselfish_urchin_density[total_num_days]
#
}#end simulation run loop
simulation_run_stamp<-cellular_automata_params$Simulation_Run_Stamp[1]
simulation_data_file_path<-paste0(cellular_automata_simulation_data_file_path,"_",simulation_run_stamp,".csv")
write.csv(cellular_automata_simulation_data,simulation_data_file_path,row.names = FALSE)
#print("current cell id")
#print(curr_cell)
#print("total urchins current cell")
#print(urchin_density_curr_cell)
#print("total aggregating urchins current cell")
#print(total_aggregating_urchins_curr_cell)
#print("total homing urchins current cell")
#print(total_homing_urchins_current_cell)
