## Default Settings used in AnalysisParam for user adjustments. 


|Setting                                |Value                                      |Function                                                   |
| ------------------------------------- | :---------------------------------------: | --------------------------------------------------------- |
|neutral_loss_calc_threshold            |350                                        |Neutral loss calculation max Δ m/z 				|
|neutral_loss_FA_lib_filename           |BBLib_FA_Neutral_Loss.csv                  |Building block library (BBLib) for fatty acid neutral loss |
|neutral_loss_HG_lib_filename           |BBLib_Common_Neutral_Loss.csv              |Building block library (BBLib) for headgroup neutral loss  |
|anion_HG_lib_filename	                |BBLib_Common_Anions.csv	                	|Building block library (BBLib) for headgroup anions		|
|level_2_lib_filename	                |Lib_Level_2_ID.csv	                        |Level 2 identification format					|
|fragment_rules_filename	          |LipidAIDER_FragRule.csv                    |Fragmentation rules
|fragments_lib	                      |TFLib_ACZX_Fragment_Sugars.csv             |Theoretical fragment library (TFLib)
|HG_chemical_composition_filename       |HG.d_Formula.csv	                        |Headgroup formula for chemical composition calculation
|score_cut_off_filename	                |filter_lipid_hit_min_score_exceptions.csv	|Optional setitings for lipid A-type dependent score cut-off
|filter_FA_NL_Da_threshold	          |250	                                    |m/z threshold to filter data for neutral loss determination. Ions below this cutoff will not be used
|neutral_loss_mode	                |mDa	                                    |Threshold type for neutral loss matching. Accepts “mDa” or “ppm”
|neutral_loss_threshold	                |15	                                    |Numerical threshold for neutral loss matching
|fragments_type                         |04A2,B,C,02A2		                  |Fragments used for TFLib
|fragments_mode	                      |ppm	                                    |Threshold type for fragments matching. Accepts “mDa” or “ppm”
|fragments_threshold	                |30	                                    |Numerical threshold for fragments matching
|filter_ion_data_threshold	          |55	                                    |Number of ions to consider
|filter_ion_data_mode	                |Top Ions	                              |Accepts “Top Ions” only
|filter_NL_threshold	                |10	                                    |Threshold frequency (%) for neutral losses detected
|filter_lipid_components_threshold	    |5                                          |threshold frequency (%) for lipid A components detected
|filter_top_lipid_hits	                |5	                                    |No. of top hits to save
|filter_lipid_hit_min_score	          |35	                                    |Minimum score of identified lipid A
|filter_lipid_hit_min_peaks	          |4	                                    |Mininum number of peaks matched, where proposed structures with less than n+y peaks are removed. n is the number of labile components in the hit and y is a user-defined parameter (default: 4), which can be from either the backbone fragments or the headgroups. 
|LpxA_threshold	                      |2	                                    |Maximum difference in carbon count between C3’ FA and C3 FA
|LpxD_threshold	                      |0	                                    |Maximum difference in carbon count between C2’ FA and C2 FA
|carbon_count_FA_lower_threshold	    |10	                                    |Minimum carbon count for all acyls
|carbon_count_FA_upper_threshold	    |16	                                    |Maximum carbon count for all acyls
|non_labile_secondary_FA	          |C'2 Sec FA	                              |Name of non-labile secondary fatty acid(s). All other secondary acids and O-linked primary acids are assumed to be labile
|missing_fragment_score_penalty	    |5	                                    |Base score penalty
|backbone_fragment_penalty_multiplier   |20	                                    |Penalty multiplier for backbone fragments
|frag_types_with_penalty	 	    |                                           |Fragments that always have a score penalty regardless of lipid A subtype
|RT_bin_window	                      |18                                         |Retention time binning. Applicable to DDA data with no restrictions. "1" - Ion with max height/ area within the bin window. "0" - Ions with intensity < max within the bin window. 
|filter_penalty	                      |30	                                    |Penalty for absence of first neutral loss, a spectral quality parameter
|total_intensity_minimum_cutoff	    |500	                                    |Total peak area/ intensity cut-off, a spectral quality parameter.
|peaks_considered_top_n_ratio	          |0.7	                                    |Proportion of peaks analysed, a spectral quality parameter


