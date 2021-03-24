# Proj_PET_fMRI_Munich_Naples

Code: 

	main_code_paper_Palombit_et_al_Neuroimage.m  
		main code to run the analyses and generate the figures from the paper "Variability of regional glucose metabolism and the topology of functional networks in the human brain" by Palombit et al.

Data:

	dataset_A.mat (Munich)
		'FC_group_dataset_A': average FC (347x347)
		'SUVR_group_dataset_A': average SUVR (347x1)
		
	dataset_B.mat (Naples)
		'FC_group_dataset_B': average FC (347x347)
		'SUVR_group_dataset_B': average SUVR (347x1)

	dataset.mat  (Munich + Naples)
		'FC_group': average FC (347x347)
		'SUVR_group': average SUVR (347x1)
		'Euclidean_distance_group': average distance (347x347)

	networks_info_GL_subcort.mat
		'net_assignment_Gordon_table': table
		'net_assignment': name of RSN to which ROI is assigned (347x1)
		'net_separation': RSN label assigned to each ROI (347x1)
		'labels_ord': ordering of GL parcels according to RSNs
		'net_names_Gordon': names of RSNs (GL atlas + SUB)

Results:

	Fig_2A_boxplot_SUVR_RSNs.png                                        
	Fig_3B_SUVR_vs_strength_scatter.png                                 
	Fig_4A_FC_short_long_range_within_between_RSNs.png                  
	Fig_4B_scatter_SUVR_vs_short_long_range_within_between_RSNs.png     
	Fig_5_STR_SUVR_association_at_increasing_DEG_percentiles.png        
	Fig_6A_connector_provincial_HUBs.png                                
	Fig_7_SUVR_graph_metrics_association_connector_provincial_HUBs.png  
	SF_1A_SUVR_dataset_A_B_merged.png                                   
	SF_1B_SUVR_dataset_A_B_merged_distributions.png                     
	SF_2_FC_dataset_A_B_merged.png                                      
	SF_3B_boxplot_SUVR_intrinsic_vs_extrinsic.png                       
	SF_4A_SUVR_vs_degree_scatter.png                                    
	SF_4B_SUVR_vs_EC_scatter.png 
	Results.mat                                                         

Required toolboxes/utilities:

	Brain Connectivity Toolbox (BCT): https://sites.google.com/site/bctnet/
	FSLNets: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/FSLNets
	FDR Benjamini-Hockberg: https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh
	Multiple Comparison Correction: https://www.mathworks.com/matlabcentral/fileexchange/34920-mult_comp_perm_corr
	dsqrform.m
		to generate a vector of elements of the upper triangle of a square symmetric matrix
		from SymMat: Symmetric Matrix Library, written by Cheol E Han (cheolhan@gmail.com) Jan 30, 2013
	sqrform.m
		to generate a square symmetric matrix from a vector of elements of the upper triangle 
		from SymMat: Symmetric Matrix Library, written by Cheol E Han (cheolhan@gmail.com) Jan 30, 2013
