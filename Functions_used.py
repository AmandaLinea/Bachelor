#Using the functions from my PUK
from Functions_code import csv_to_two_heatmaps
from Functions_code import csv_to_histogram
from Functions_code import sequence_from_csv
#from Functions_code import clinical_var_csv_to_list_with_stability
from Functions_code import clinvar_scatterplot_ddg_and_dddg
from Clinical_function import clinical_var_csv_to_list_with_stability
from Clinical_function import clinical_var_scatterplot_ddg_and_dddg_2
from Clinical_function import clinical_var_scatterplot_ddg_and_dddg_1
from Functions_code import scatterplotting_old_and_new_predictions
from Functions_code import scatterplotting_old_vs_new_ddg_predictions
from Functions_code import scatterplotting_old_vs_new_dddg_predictions
from Functions_code import remove_wildtype_from_csv
from Functions_code import gathering_all_csv_files_into_one
from Clinical_function import scatter_marginal_quadrant_plot
from Functions_code import sequence_from_nicolas
from Functions_code import write_prismfile_from_RaSP_csv
from Functions_code import write_prismfile_from_Nicolas_csv
from Functions_code import read_prismfiles_into_dataframe
from Functions_code import big_dataframe_to_significant_dataframe
from Functions_code import scatterplot_marginals_quadrant_for_significant_dataframe
from Functions_code import max_plotted
from Functions_code import get_p_and_g_origo

'''csv file to heatmaps'''
#toO8B = csv_to_two_heatmaps(f'Old_predictions\old_pred.csv', '2PAH', 'score_ml',chainA=True, chainB=True, chainC=False, chainD=False, show_A=True, show_B=True, show_C=False, show_D=False)
#tO8B = csv_to_two_heatmaps(f'Old_predictions\old_pred.csv', '2PAH', 'score_ml_ddg_bind',chainA=True, chainB=True, chainC=False, chainD=False, show_A=True, show_B=True, show_C=False, show_D=False)
#yO8B = csv_to_two_heatmaps(f'_New_old_predictions\_New_old_pred.csv', '2PAH', 'score_ml',chainA=True, chainB=True, chainC=False, chainD=False, show_A=True, show_B=True, show_C=False, show_D=False)
#uO8B = csv_to_two_heatmaps(f'_New_old_predictions\_New_old_pred.csv', '2PAH', 'score_ml_ddg_bind',chainA=True, chainB=True, chainC=False, chainD=False, show_A=True, show_B=True, show_C=False, show_D=False)

#denheatmap1_2 = csv_to_two_heatmaps(f'_New_old_predictions\_New_old_pred.csv', '5DEN', 'score_ml',chainA=True, chainB=True, chainC=False, chainD=False, show_A=True, show_B=True, show_C=False, show_D=False)
#dendddgheatmap3_4 = csv_to_two_heatmaps(f'_New_old_predictions\_New_old_pred.csv', '5DEN', 'score_ml_ddg_bind',chainA=False, chainB=False, chainC=True, chainD=True, show_A=False, show_B=False, show_C=True, show_D=True)
#denheatmap1_2 = csv_to_two_heatmaps(f'_New_old_predictions\_New_old_pred.csv', '5DEN', 'score_ml_ddg_bind',chainA=True, chainB=True, chainC=False, chainD=False, show_A=True, show_B=True, show_C=False, show_D=False)
#dendddgheatmap3_4 = csv_to_two_heatmaps(f'_New_old_predictions\_New_old_pred.csv', '5DEN', 'score_ml',chainA=False, chainB=False, chainC=True, chainD=True, show_A=False, show_B=False, show_C=True, show_D=True)

#PAH = csv_to_two_heatmaps(f'{pdbid}_df_ml.csv', {pdbid}, 'score_ml',chainA=True, chainB=True, chainC=False, chainD=False, show_A=True, show_B=True, show_C=False, show_D=False)

'''csv file to histograms'''
#toO8Bddghistogram = csv_to_histogram(r'ddG_bind_RaSP.csv', '2o8b', 'score_ml_ddg_bind',1,8)
#a=csv_to_histogram(f'Old_predictions\old_pred.csv', '2O8B', 'score_ml_ddg_bind',-12,12, chainC=False,chainD=False)
##b=csv_to_histogram(f'Old_predictions\old_pred.csv', '2O8B', 'score_ml',-12,12,chainC=False,chainD=False)
#c=csv_to_histogram(f'_New_old_predictions\_New_old_pred.csv', '2O8B', 'score_ml_ddg_bind',-12,12, chainC=False,chainD=False)
##d=csv_to_histogram(f'_New_old_predictions\_New_old_pred.csv', '2O8B', 'score_ml',-12,12, chainC=False,chainD=False)
#e=csv_to_histogram(f'Old_predictions\old_pred.csv', '2O8B', 'score_ml_ddg_bind',1,12, chainC=False,chainD=False)
##b=csv_to_histogram(f'Old_predictions\old_pred.csv', '2O8B', 'score_ml',-12,12,chainC=False,chainD=False)
#f=csv_to_histogram(f'_New_old_predictions\_New_old_pred.csv', '2O8B', 'score_ml_ddg_bind',1,12, chainC=False,chainD=False)
#dendddghistogram = csv_to_histogram(r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv', '5DEN', 'score_ml_ddg_bind',2,8)

'''sequence from csv'''
#denseq = sequence_from_csv(r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv','5DEN')
#print(sequence_from_csv(f'predictions\All_new_df_ml.csv','2O8B','A'))

'''clinvar scatterplots'''
#mutsalpha = clinvar_scatterplot_ddg_and_dddg(r'ShinyClinVar_MSH2_2024-01-04.csv', r'ShinyClinVar_MSH2_benign_2024-01-04.csv', r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv', '2O8B',r'ShinyClinVar_MSH6_2024-01-04.csv', r'ShinyClinVar_MSH6_benign_2024-01-04.csv', r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv', '2O8B')

'''list of stabilities for clinical significant variants'''
#den_clinvar_variants = clinical_var_csv_to_list_with_stability(r'ShinyClinVar_PAH_2024-01-04.csv',r'ShinyClinVar_PAH_benign_2024-01-04 (1).csv',r'De-vigtige-filer\RaSP_of_2O8B_and_5DEN.csv', '5DEN')
#pah_clinical_var_stability = clinical_var_csv_to_list_with_stability(r'__amanda-project__data__v1__jonsson__2024.csv',"P00439", r'ddG_bind_RaSP.csv', "tetramere")

'''Nicolas clinical variant scatterplots'''
#pah_scatterplot = clinical_var_scatterplot_ddg_and_dddg_1(f'__amanda-project__data__v1__jonsson__2024.csv', "P00439", f'clean_OPAH_df_ml.csv', "OPAH")
#scatter_marginal_quadrant_plot(f'__amanda-project__data__v1__jonsson__2024.csv', "P00439", f'clean_OPAH_df_ml.csv', "OPAH", 2, 1)

'''old vs new data, in scatterplots'''
#scatterplotting_old_and_new_predictions(f'clean_2PAH_df_ml.csv', '2PAH', f'clean_OPAH_df_ml.csv', 'OPAH')
#scatterplotting_old_vs_new_ddg_predictions(f'Old_predictions\old_pred.csv', '2PAH', f'predictions\OPAH_df_ml.csv', 'OPAH')
#scatterplotting_old_vs_new_dddg_predictions(f'Old_predictions\old_pred.csv', '2PAH', f'predictions\OPAH_df_ml.csv', 'OPAH')
#
#scatterplotting_old_and_new_predictions(f'clean_2O8B_and_5DEN_df_ml.csv', '5DEN', f'clean_OPAH_df_ml.csv', 'OPAH')
#scatterplotting_old_vs_new_ddg_predictions(f'Old_predictions\old_pred.csv', '5DEN', f'predictions\OPAH_df_ml.csv', 'OPAH')
#scatterplotting_old_vs_new_dddg_predictions(f'Old_predictions\old_pred.csv', '5DEN', f'predictions\OPAH_df_ml.csv', 'OPAH')
#
#scatterplotting_old_and_new_predictions(f'clean_2O8B_and_5DEN_df_ml.csv', '2O8B', f'clean_2O8B_df_ml.csv', '2O8B')
#scatterplotting_old_vs_new_ddg_predictions(f'Old_predictions\old_pred.csv', '2O8B', f'_New_old_predictions\_New_old_pred.csv', '2O8B')
#scatterplotting_old_vs_new_dddg_predictions(f'Old_predictions\old_pred.csv', '2O8B', f'_New_old_predictions\_New_old_pred.csv', '2O8B')
#
#scatterplotting_old_and_new_predictions(f'clean_2PAH_df_ml.csv', '2PAH', f'clean_2PAHnew_df_ml.csv', '2PAH')
#scatterplotting_old_vs_new_ddg_predictions(f'Old_predictions\old_pred.csv', '2PAH', f'_New_old_predictions\_New_old_pred.csv', '2PAH')
#scatterplotting_old_vs_new_dddg_predictions(f'Old_predictions\old_pred.csv', '2PAH', f'_New_old_predictions\_New_old_pred.csv', '2PAH')
#
#scatterplotting_old_and_new_predictions(f'clean_2O8B_and_5DEN_df_ml.csv', '5DEN', f'clean_5DENnew_df_ml.csv', '5DEN')
#scatterplotting_old_vs_new_ddg_predictions(f'Old_predictions\old_pred.csv', '5DEN', f'_New_old_predictions\_New_old_pred.csv', '5DEN')
#scatterplotting_old_vs_new_dddg_predictions(f'Old_predictions\old_pred.csv', '5DEN', f'_New_old_predictions\_New_old_pred.csv', '5DEN')



'''remove wildtypes from RaSP csv file'''
#remove_wildtype_from_csv(f'All_df_ml.csv', 'df_ml_all')
#remove_wildtype_from_csv(f'RaSP_of_2O8B_and_5DEN.csv', '2O8B_and_5DEN_df_ml')
#remove_wildtype_from_csv(f'2O8B_df_ml.csv', '2O8B_df_ml')
#remove_wildtype_from_csv('New_old_predictions\_5DEN_df_ml.csv', '5DENnew_df_ml')
#remove_wildtype_from_csv('predictions\All_new_df_ml.csv', 'all_new_df_ml')

'''gathering all csv files into one'''
#gathering_all_csv_files_into_one(f'New_old_predictions','new_old_pred')

'''extracting the sequence from nicolas'''
#print(sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'O95786'))


'''Write prism files from RaSP'''
#write_prismfile_from_RaSP_csv()


'''Write prism files from nicolas'''
#write_prismfile_from_Nicolas_csv()


'''read the prismfiles into a dataframe'''
#read_prismfiles_into_dataframe()

'''Use the big dataframe with alligned variants to get a small dataframe for just the significant variants'''
#big_dataframe_to_significant_dataframe()

'''Create a scatterplot with margins for the mean stabilities'''
#scatterplot_marginals_quadrant_for_significant_dataframe(2,0.5)


'''Create a scatterplot with margins for the maximum stabilities'''
#max_plotted(2,0.5)

'''Get the information about the glycine substitution and proline variants'''
#get_p_and_g_origo(f'Old_predictions\old_pred.csv', '2PAH', f'_New_old_predictions\_New_old_pred.csv', '2PAH', 'proline_glycine_2PAH_results')
