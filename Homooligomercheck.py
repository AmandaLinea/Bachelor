import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from PrismData import PrismParser, VariantData
from Functions_code import sequence_from_csv

# This function flattens a nested list like [[0], [1], [2], [3], [4]] into just [0, 1, 2, 3, 4] 
def flatten(list):
    return [x for xs in list for x in xs]

# Converting the clinvar variant data into the same format as the RaSP csv file data
def variant_and_related_disease(variantcsvfile, uniprotnr , clinicalsigni):
    df = pd.read_csv(variantcsvfile)
    df = df[[ "uniprot", "variant", "resi" , "clinvar_signifiance"]]
    df = df[df["uniprot"] == uniprotnr]
    df = df[df["clinvar_signifiance"]== clinicalsigni]
    df = df.sort_values("resi")
    df = df[[ "variant"]]
    list_of_data = df.values.tolist()
    return list_of_data

def opah_chains_against_chains(chain1, chain2):
    df = pd.read_csv(f'clean_OPAH_df_ml.csv')
    chainA = df[df['chainid']== chain1]
    chainB = df[df['chainid']== chain2]
    chainC = df[df['chainid']== 'C']
    chainD = df[df['chainid']== 'D']

     # Step 1: Filter rows with matching 'variant' and 'chainid'
    matching_rows = chainA.merge(chainB, on=['variant'], suffixes=(f'_{chain1}', f'_{chain2}'))

    # Step 2: Extract desired columns
    score_ml_A = matching_rows[f'score_ml_ddg_bind_{chain1}']
    score_ml_B = matching_rows[f'score_ml_ddg_bind_{chain2}']

    # Create the scatter plot
    sns.scatterplot(x=score_ml_B , y=score_ml_A, data=matching_rows)

    ## Step 3: Create scatterplot using Seaborn
    #plt.figure(figsize=(8, 6))
    #sns.scatterplot(x=chainC, y=chainA, hue=df['score_ml'])
#
    ## Step 4: Fit a linear regression model using numpy
    slope, intercept = np.polyfit(score_ml_B, score_ml_A, 1)
    trendline_x = np.linspace(min(score_ml_A), max(score_ml_A), 100)
    trendline_y = slope * trendline_x + intercept
    plt.plot(trendline_x, trendline_y, color='red', label='Trendline')
#
    ## Step 5: Display R-value (calculate it separately)
    r_value = np.corrcoef(score_ml_B, score_ml_A)[0, 1]
    plt.text(0.7, 0.1, f'R = {r_value:.2f}', transform=plt.gca().transAxes, fontsize=12)

    plt.xlabel(f'OPAH chain {chain2}')
    plt.ylabel(f'OPAH chain {chain1}')
    plt.title(f'Scatterplot for DDG for OPAH chain {chain2} and {chain1} against eachother')
    plt.legend()
    plt.show()

def opah_chains_against_chains_mean_std():
    df = pd.read_csv('clean_OPAH_df_ml.csv')

    # Calculate mean and standard deviation for 'score_ml' grouped by 'variant'
    mean_scores = df.groupby('variant')['score_ml'].mean()
    std_scores = df.groupby('variant')['score_ml'].std()

    # Create a bar plot with error bars
    plt.figure(figsize=(10, 6))
    sns.scatterplot(mean_scores, x=mean_scores.index, y=mean_scores.values, hue=std_scores )
    plt.xlabel('Variant')
    plt.ylabel('Mean score_ml')
    plt.title(f'Mean score_ml with Standard Deviation for OPAH')
    plt.xticks(rotation=45)
    plt.show()

# Example usage:
#opah_chains_against_chains_mean_std()

#opah_chains_against_chains('A','D')
#opah_chains_against_chains('B','D')
#opah_chains_against_chains('C','D')
#opah_chains_against_chains('B','C')

def scatter_marginal_quadrant_plot(variantcsvpath1, uniprotnr1, RaSPdatacsvpath1, pdbid1, ddg_cutoff, dddg_cutoff):
    df = pd.read_csv(RaSPdatacsvpath1)
    df = df[df['pdbid'] == pdbid1]

    # Getting the clinvar variants as a list
    list1 = flatten(variant_and_related_disease(variantcsvpath1, uniprotnr1 ,"pathogenic"))
    list2 = flatten(variant_and_related_disease(variantcsvpath1, uniprotnr1 ,"benign"))

    # Filter rows based on variant lists
    filtered_df = df[df["variant"].isin(list1 + list2)]

    # Create the 'clinical_signi' column
    filtered_df["clinical_signi"] = filtered_df["variant"].apply(
        lambda x: "Pathogenic" if x in list1 else "Benign" 
        )

    # Calculate average values
    filtered_df["average_ddg_values"] = filtered_df.groupby("variant")["score_ml"].mean()
    filtered_df["average_dddg_values"] = filtered_df.groupby("variant")["score_ml_ddg_bind"].mean()

    # Check for identical variants
    duplicates = filtered_df.duplicated(subset=["variant"], keep=False)

    # Create quadrants based on x and y values
    #filtered_df["quadrant"] = ""
    #filtered_df.loc[(filtered_df["score_ml"] < ddg_cutoff) & (filtered_df["score_ml_ddg_bind"] > dddg_cutoff), "quadrant"] = "I"
    #filtered_df.loc[(filtered_df["score_ml"] > ddg_cutoff) & (filtered_df["score_ml_ddg_bind"] > dddg_cutoff), "quadrant"] = "II"
    #filtered_df.loc[(filtered_df["score_ml"] < ddg_cutoff) & (filtered_df["score_ml_ddg_bind"] <= dddg_cutoff), "quadrant"] = "III"
    #filtered_df.loc[(filtered_df["score_ml"] > ddg_cutoff) & (filtered_df["score_ml_ddg_bind"] <= dddg_cutoff), "quadrant"] = "IV"

    # Create quadrants based on x and y values
    filtered_df["quadrant"] = ""
    filtered_df.loc[(filtered_df["average_ddg_values"] < ddg_cutoff) & (filtered_df["average_dddg_values"] > dddg_cutoff), "quadrant"] = "I"
    filtered_df.loc[(filtered_df["average_ddg_values"] > ddg_cutoff) & (filtered_df["average_dddg_values"] > dddg_cutoff), "quadrant"] = "II"
    filtered_df.loc[(filtered_df["average_ddg_values"] < ddg_cutoff) & (filtered_df["average_dddg_values"] <= dddg_cutoff), "quadrant"] = "III"
    filtered_df.loc[(filtered_df["average_ddg_values"] > ddg_cutoff) & (filtered_df["average_dddg_values"] <= dddg_cutoff), "quadrant"] = "IV"

    dfpath = filtered_df[filtered_df["clinical_signi"]=='Pathogenic']
    dfbenign = filtered_df[filtered_df["clinical_signi"]=='Benign']
    all_I = filtered_df["quadrant"].value_counts()["I"]
    path_I = dfpath[dfpath["quadrant"]== "I"]
    benign_I = dfbenign[dfbenign["quadrant"]=="I"]

    all_II = filtered_df["quadrant"].value_counts()["II"]
    path_II = dfpath[dfpath["quadrant"]=="II"]
    benign_II = dfbenign[dfbenign["quadrant"]== "II"]

    all_III = filtered_df["quadrant"].value_counts()["III"]
    path_III = dfpath[dfpath["quadrant"]=="III"]
    benign_III = dfbenign[dfbenign["quadrant"]=="III"]

    all_IV = filtered_df["quadrant"].value_counts()["IV"]
    path_IV = dfpath[dfpath["quadrant"]=="IV"]
    benign_IV = dfbenign[dfbenign["quadrant"]=="IV"]

    all = len(filtered_df["quadrant"])
    path_all = len(dfpath)
    benign_all = len(dfbenign)

    # Create the Seaborn jointplot
    #sns.color_palette("Spectral", as_cmap=True)
    g = sns.jointplot(data=filtered_df, x="average_ddg_values", y="average_dddg_values", hue="clinical_signi", alpha=.8, palette="Spectral")
    g.refline(x=ddg_cutoff, y=dddg_cutoff, marginal=False)
    g.figure.subplots_adjust(top=.9)
    g.figure.suptitle(f'{pdbid1} Scatterplot with marginals, showing the distribution \nof pathogenic and benign clinical variants')
    plt.text(-3,2, f"I \nDots: {all_I}/{all}={round(all_I/all, 3)}\nPathogenic: {len(path_I)}/{path_all}={round(len(path_I)/path_all, 3)}\nBenign: {len(benign_I)}/{benign_all}={round(len(benign_I)/benign_all, 3)}")
    plt.text(9,2, f"II \nDots: {all_II}/{all}={round(all_II/all,3)}\nPathogenic: {len(path_II)}/{path_all}={round(len(path_II)/path_all,3)}\nBenign: {len(benign_II)}/{benign_all}={round(len(benign_II)/benign_all, 3)}")
    plt.text(-3,-2, f"III \nDots: {all_III}/{all}={round(all_III/all, 3)}\nPathogenic: {len(path_III)}/{path_all}={round(len(path_III)/path_all,3)}\nBenign: {len(benign_III)}/{benign_all}={round(len(benign_III)/benign_all,3)}")
    plt.text(9,-2, f"IV \nDots: {all_IV}/{all}={round(all_IV/all,3)}\nPathogenic: {len(path_IV)}/{path_all}={round(len(path_IV)/path_all,3)}\nBenign: {len(benign_IV)}/{benign_all}={round(len(benign_IV)/benign_all,3)}")
    plt.show()

def update_score(df):
    # Group by the 'variant' column and calculate the mean of both 'score_ml' and 'score_ml_ddg_bind'
    #grouped_df = df.groupby("variant")[["score_ml","score_ml_ddg_bind"]].mean()
    dfp2 = df[df["score_ml_ddg_bind"]>= 0]
    scoreml_df = dfp2.groupby("variant")[["score_ml"]].mean()
    scoredddg_df = dfp2.groupby("variant")[["score_ml_ddg_bind"]].mean()

    # Merge the original dataframe with the grouped dataframe
    middle_df = pd.merge(scoreml_df, scoredddg_df, on='variant')

    updated_df = pd.merge(df, middle_df, on = 'variant', suffixes=('', '_m'))
    # Replace the old 'score_ml' and 'score_ml_ddg_bind' with the calculated averages
    updated_df['score_ml'] = updated_df['score_ml_m']
    updated_df['score_ml_ddg_bind'] = updated_df['score_ml_ddg_bind_m']

    # Drop duplicate variants (keep only the first occurrence)
    updated_df.drop_duplicates(subset='variant', inplace=True)

    return updated_df


#df = pd.read_csv(f'clean_df_ml_all.csv')
## Apply the function to each dataframe segment
#dfOPAH = df[df['pdbid'] == 'OPAH']
#print(dfOPAH[dfOPAH['variant']=='V118A'])
#dfOPAH = update_score(dfOPAH)
#print(dfOPAH)

def scatter_marginal_quadrant_plot_for_all(variantcsvpath1, RaSPdatacsvpath1, ddg_cutoff, dddg_cutoff):
    df = pd.read_csv(RaSPdatacsvpath1)
    # Apply the function to each dataframe segment
    dfOPAH = df[df['pdbid'] == 'OPAH']
    df2O8B = df[df['pdbid'] == '2O8B']
    df2DN1 = df[df['pdbid'] == '2DN1']
    df2O4H = df[df['pdbid'] == '2O4H']
    df4OL0 = df[df['pdbid'] == '4OL0']
    df8G7V = df[df['pdbid'] == '8G7V']
    df1G82 = df[df['pdbid'] == '1G82']
    df1I85 = df[df['pdbid'] == '1I85']

    dfOPAH = update_score(dfOPAH)
    df2O8B = update_score(df2O8B)
    df2DN1 = update_score(df2DN1)
    df2O4H = update_score(df2O4H)
    df4OL0 = update_score(df4OL0)
    df8G7V = update_score(df8G7V)
    df1G82 = update_score(df1G82)
    df1I85 = update_score(df1I85)

    # Getting the clinically significant variants as a list
    listOPAHpath = flatten(variant_and_related_disease(variantcsvpath1, 'P00439' ,"pathogenic"))
    listOPAHbeni = flatten(variant_and_related_disease(variantcsvpath1, 'P00439' ,"benign"))
    list2O8Bpath1 = flatten(variant_and_related_disease(variantcsvpath1, 'P43246' ,"pathogenic"))
    list2O8Bbeni1 = flatten(variant_and_related_disease(variantcsvpath1, 'P43246' ,"benign"))
    list2O8Bpath2 = flatten(variant_and_related_disease(variantcsvpath1, 'P52701' ,"pathogenic"))
    list2O8Bbeni2 = flatten(variant_and_related_disease(variantcsvpath1, 'P52701' ,"benign"))
    list2DN1path1 = flatten(variant_and_related_disease(variantcsvpath1, 'P69905' ,"pathogenic"))
    list2DN1beni1 = flatten(variant_and_related_disease(variantcsvpath1, 'P69905' ,"benign"))
    list2DN1path2 = flatten(variant_and_related_disease(variantcsvpath1, 'P68871' ,"pathogenic"))
    list2DN1beni2 = flatten(variant_and_related_disease(variantcsvpath1, 'P68871' ,"benign"))
    list2O4Hpath = flatten(variant_and_related_disease(variantcsvpath1, 'P45381' ,"pathogenic"))
    list2O4Hbeni = flatten(variant_and_related_disease(variantcsvpath1, 'P45381' ,"benign"))
    list4OL0path1 = flatten(variant_and_related_disease(variantcsvpath1, 'Q9Y5L0 ' ,"pathogenic"))
    list4OL0beni1 = flatten(variant_and_related_disease(variantcsvpath1, 'Q9Y5L0 ' ,"benign"))
    list4OL0path2 = flatten(variant_and_related_disease(variantcsvpath1, 'P62826 ' ,"pathogenic"))
    list4OL0beni2 = flatten(variant_and_related_disease(variantcsvpath1, 'P62826 ' ,"benign"))
    list8G7Vpath1 = flatten(variant_and_related_disease(variantcsvpath1, 'O95786' ,"pathogenic"))
    list8G7Vbeni1 = flatten(variant_and_related_disease(variantcsvpath1, 'O95786' ,"benign"))
    list8G7Vpath2 = flatten(variant_and_related_disease(variantcsvpath1, 'Q8IUD6' ,"pathogenic"))
    list8G7Vbeni2 = flatten(variant_and_related_disease(variantcsvpath1, 'Q8IUD6' ,"benign"))
    list1G82path = flatten(variant_and_related_disease(variantcsvpath1, 'P31371' ,"pathogenic"))
    list1G82beni = flatten(variant_and_related_disease(variantcsvpath1, 'P31371' ,"benign"))
    list1I85path1 = flatten(variant_and_related_disease(variantcsvpath1, 'P16410' ,"pathogenic"))
    list1I85beni1 = flatten(variant_and_related_disease(variantcsvpath1, 'P16410' ,"benign"))
    list1I85path2 = flatten(variant_and_related_disease(variantcsvpath1, 'P42081' ,"pathogenic"))
    list1I85beni2 = flatten(variant_and_related_disease(variantcsvpath1, 'P42081' ,"benign"))

    all_lists = [
    listOPAHpath, listOPAHbeni,
    list2O8Bpath1, list2O8Bbeni1, list2O8Bpath2, list2O8Bbeni2,
    list2DN1path1, list2DN1beni1, list2DN1path2, list2DN1beni2,
    list2O4Hpath, list2O4Hbeni,
    list4OL0path1, list4OL0beni1, list4OL0path2, list4OL0beni2,
    list8G7Vpath1, list8G7Vbeni1, list8G7Vpath2, list8G7Vbeni2,
    list1G82path, list1G82beni,
    list1I85path1, list1I85beni1, list1I85path2, list1I85beni2
    ]

    patho_lists = [
    listOPAHpath, list2O8Bpath1, list2O8Bpath2,
    list2DN1path1, list2DN1path2, list2O4Hpath,
    list4OL0path1, list4OL0path2, list8G7Vpath1,
    list8G7Vpath2, list1G82path, list1I85path1,
    list1I85path2
]

    # Filter rows based on variant lists
    filtered_df = dfOPAH[dfOPAH["variant"].isin(listOPAHpath + listOPAHbeni)]
    filtered_df = pd.concat([filtered_df, df2O8B[df2O8B["variant"].isin(list2O8Bpath1+list2O8Bbeni1+list2O8Bpath2+list2O8Bbeni2)], 
                             df2DN1[df2DN1["variant"].isin(list2DN1path1 + list2DN1beni1 + list2DN1path2 + list2DN1beni2)],
                             df2O4H[df2O4H["variant"].isin(list2O4Hpath + list2O4Hbeni)],
                             df4OL0[df4OL0["variant"].isin(list4OL0path1 + list4OL0beni1 + list4OL0path2 + list4OL0beni2)],
                             df8G7V[df8G7V["variant"].isin(list8G7Vpath1 + list8G7Vbeni1 + list8G7Vpath2 + list8G7Vbeni2)],
                             df1G82[df1G82["variant"].isin(list1G82path + list1G82beni)],
                             df1I85[df1I85["variant"].isin(list1I85path1 + list1I85beni1 + list1I85path2 + list1I85beni2)]])

    filtered_df["clinical_signi"] = filtered_df["variant"].apply(
        lambda x: "Pathogenic" if any(x in patho_list for patho_list in patho_lists) else "Benign"
        )
    
    print(filtered_df)

    ## Create the 'clinical_signi' column
    #filtered_df["clinical_signi"] = filtered_df["variant"].apply(
    #    lambda x: "Pathogenic" if x in patho_list else "Benign" 
    #    )

    # Create quadrants based on x and y values
    filtered_df["quadrant"] = ""
    filtered_df.loc[(filtered_df["score_ml"] < ddg_cutoff) & (filtered_df["score_ml_ddg_bind"] > dddg_cutoff), "quadrant"] = "I"
    filtered_df.loc[(filtered_df["score_ml"] > ddg_cutoff) & (filtered_df["score_ml_ddg_bind"] > dddg_cutoff), "quadrant"] = "II"
    filtered_df.loc[(filtered_df["score_ml"] < ddg_cutoff) & (filtered_df["score_ml_ddg_bind"] <= dddg_cutoff), "quadrant"] = "III"
    filtered_df.loc[(filtered_df["score_ml"] > ddg_cutoff) & (filtered_df["score_ml_ddg_bind"] <= dddg_cutoff), "quadrant"] = "IV"


    dfpath = filtered_df[filtered_df["clinical_signi"]=='Pathogenic']
    dfbenign = filtered_df[filtered_df["clinical_signi"]=='Benign']

    all_I = filtered_df["quadrant"].value_counts()["I"]
    path_I = dfpath[dfpath["quadrant"]== "I"]
    benign_I = dfbenign[dfbenign["quadrant"]=="I"]

    all_II = filtered_df["quadrant"].value_counts()["II"]
    path_II = dfpath[dfpath["quadrant"]=="II"]
    benign_II = dfbenign[dfbenign["quadrant"]== "II"]

    all_III = filtered_df["quadrant"].value_counts()["III"]
    path_III = dfpath[dfpath["quadrant"]=="III"]
    benign_III = dfbenign[dfbenign["quadrant"]=="III"]

    all_IV = filtered_df["quadrant"].value_counts()["IV"]
    path_IV = dfpath[dfpath["quadrant"]=="IV"]
    benign_IV = dfbenign[dfbenign["quadrant"]=="IV"]

    all = len(filtered_df["quadrant"])
    path_all = len(dfpath)
    benign_all = len(dfbenign)

    # Create the Seaborn jointplot
    #sns.color_palette("Spectral", as_cmap=True)
    g = sns.jointplot(data=filtered_df, x="score_ml", y="score_ml_ddg_bind", hue="clinical_signi", alpha=1, palette="Spectral")
    g.refline(x=ddg_cutoff, y=dddg_cutoff, marginal=False)
    g.figure.subplots_adjust(top=.9)
    g.figure.suptitle(f'Scatterplot with marginals, showing the distribution \nof pathogenic and benign clinical variants')
    plt.text(-3.5,2, f"I \nDots: {all_I}/{all}={round(all_I/all, 3)}\nPathogenic: {len(path_I)}/{path_all}={round(len(path_I)/path_all, 3)}\nBenign: {len(benign_I)}/{benign_all}={round(len(benign_I)/benign_all, 3)}")
    plt.text(9.5,2, f"II \nDots: {all_II}/{all}={round(all_II/all,3)}\nPathogenic: {len(path_II)}/{path_all}={round(len(path_II)/path_all,3)}\nBenign: {len(benign_II)}/{benign_all}={round(len(benign_II)/benign_all, 3)}")
    plt.text(-3.5,-0.5, f"III \nDots: {all_III}/{all}={round(all_III/all, 3)}\nPathogenic: {len(path_III)}/{path_all}={round(len(path_III)/path_all,3)}\nBenign: {len(benign_III)}/{benign_all}={round(len(benign_III)/benign_all,3)}")
    plt.text(9.5,-0.5, f"IV \nDots: {all_IV}/{all}={round(all_IV/all,3)}\nPathogenic: {len(path_IV)}/{path_all}={round(len(path_IV)/path_all,3)}\nBenign: {len(benign_IV)}/{benign_all}={round(len(benign_IV)/benign_all,3)}")
    plt.show()
    

#scatter_marginal_quadrant_plot_for_all(f'nr2_Nicolas_data.csv', f'clean_df_ml_all.csv', 2, 1)

def write_prismfile_from_RaSP_csv():
    df = pd.read_csv(f"predictions\All_new_df_ml.csv" )
    dfb = df[df['pdbid']=='2C9S']
    dfb = dfb[dfb['chainid']=='F']
    dfb.reset_index(inplace=True)
    dfb = dfb[['variant', 'score_ml', 'score_ml_ddg_bind']]

    md = {'version':1,
                'protein':{'name':'Human SUPEROXIDE DISMUTASE', 'organism':'Homo sapiens', 'uniprot':'P00441','first_residue_number':1, 'sequence':'ATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ'},
                'method':{'name':'RaSP from a PDB'},
                'columns':{'score_ml':'DDG', 'score_ml_ddg_bind': 'DDDG'}}

    prismdata=VariantData(metadata=md, dataframe=dfb)
    prismdata.check()
    pp = PrismParser()
    pp.write("Prismfiles/prism_RaSP_P00441_F.txt", prismdata)

#write_prismfile_from_RaSP_csv()

def write_prismfile_from_Nicolas_csv():
    df = pd.read_csv(f'nr2_Nicolas_data.csv')
    dfa = df[df['uniprot'] =='P16410']
    dfa.reset_index(inplace=True)
    dfa = dfa[["variant", "clinvar_signifiance"]]

    md = {'version':1,
                'protein':{'name':'Human CYTOTOXIC T-LYMPHOCYTE-ASSOCIATED PROTEIN 4', 'organism':'Homo sapiens', 'uniprot':'P16410','first_residue_number':1, 'sequence':'MACLGFQRHKAQLNLATRTWPCTLLFFLLFIPVFCKAMHVAQPAVVLASSRGIASFVCEYASPGKATEVRVTVLRQADSQVTEVCAATYMMGNELTFLDDSICTGTSSGNQVNLTIQGLRAMDTGLYICKVELMYPPPYYLGIGNGTQIYVIDPEPCPDSDFLLWILAAVSSGLFFYSFLLTAVSLSKMLKKRSPLTTGVYVKMPPTEPECEKQFQPYFIPIN'},
                'method':{'name':'RaSP from a PDB'},
                'columns':{'clinvar_signifiance':'Clinical_significance'}}


    prismdata=VariantData(metadata=md, dataframe=dfa)
    prismdata.check()
    pp = PrismParser()
    pp.write("Prismfiles/prism_Nicolas_CTLA4_001.txt", prismdata)

def read_prismfiles_into_dataframe():
    # Read data section
    dataframe1 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_O95786.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe1.insert(0, 'uniprot', 'O95786')
    dataframe1['q_structure'] = 'Heterotetramer'

    dataframe2 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P00439.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe2.insert(0, 'uniprot', 'P00439')
    dataframe2['q_structure'] = 'Homotetramer'


    dataframe3 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P00441.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe3.insert(0, 'uniprot', 'P00441')
    dataframe3['q_structure'] = 'Homodimer'


    dataframe4 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P16410.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe4.insert(0, 'uniprot', 'P16410')
    dataframe4['q_structure'] = 'Heterotetramer'


    dataframe5 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P30566.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe5.insert(0, 'uniprot', 'P30566')
    dataframe5['q_structure'] = 'Homotetramer'


    dataframe6 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P31371.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe6.insert(0, 'uniprot', 'P31371')
    dataframe6['q_structure'] = 'Homotetramer'

    dataframe7 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P35520.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe7.insert(0, 'uniprot', 'P35520')
    dataframe7['q_structure'] = 'Homodimer'


    dataframe8 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P42081.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe8.insert(0, 'uniprot', 'P42081')
    dataframe8['q_structure'] = 'Heterotetramer'


    dataframe9 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P43246.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe9.insert(0, 'uniprot', 'P43246')
    dataframe9['q_structure'] = 'Heterodimer'

    dataframe10 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P45381.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe10.insert(0, 'uniprot', 'P45381')
    dataframe10['q_structure'] = 'Homodimer'


    dataframe11 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P52701.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe11.insert(0, 'uniprot', 'P52701')
    dataframe11['q_structure'] = 'Heterodimer'


    dataframe12 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_P62826.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe12.insert(0, 'uniprot', 'P62826')
    dataframe12['q_structure'] = 'Heterodimer'


    dataframe13 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_Q8IUD6.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe13.insert(0, 'uniprot', 'Q8IUD6')
    dataframe13['q_structure'] = 'Heterotetramer'


    dataframe14 = pd.read_csv('Prismfiles\prism_Nicolas_RaSP_Q9Y5L0.txt', delim_whitespace=True, comment='#', header=0,
                         keep_default_na=True, na_values=['Na','na'])
    dataframe14.insert(0, 'uniprot', 'Q9Y5L0')
    dataframe14['q_structure'] = 'Heterodimer'

    big_dataframe = pd.concat([dataframe1,dataframe2,dataframe3,dataframe4,dataframe5,dataframe6,dataframe7,dataframe8,dataframe9,dataframe10,dataframe11,dataframe12,dataframe13,dataframe14], ignore_index=True)

    big_dataframe.to_csv('big_dataframe.csv', index=False)

    return big_dataframe

#print(read_prismfiles_into_dataframe())

def big_dataframe_to_significant_dataframe():
    df = pd.read_csv('big_dataframe.csv', index_col=0)
    df.reset_index(drop=True)
    dfa = df[df['clinvar_signifiance_00'].notnull()]
    dfa.reset_index(inplace=True)
    
    # Calculate the average for 'score_ml' columns
    dfa['score_ml'] = dfa[['score_ml_01', 'score_ml_02', 'score_ml_03', 'score_ml_04']].mean(axis=1)

    # Calculate the average for 'score_ml_ddg_bind' columns
    dfa['score_ml_ddg_bind'] = dfa[['score_ml_ddg_bind_01', 'score_ml_ddg_bind_02', 'score_ml_ddg_bind_03', 'score_ml_ddg_bind_04']].mean(axis=1)

    # Calculate the standard deviation for 'score_ml' and 'score_ml_ddg_bind'
    dfa['std_score_ml'] = dfa[['score_ml_01', 'score_ml_02', 'score_ml_03', 'score_ml_04']].std(axis=1)
    dfa['std_score_ml_ddg_bind'] = dfa[['score_ml_ddg_bind_01', 'score_ml_ddg_bind_02', 'score_ml_ddg_bind_03', 'score_ml_ddg_bind_04']].std(axis=1)
    
    # Calculate the average for 'score_ml' columns
    dfa['score_ml_max'] = dfa[['score_ml_01', 'score_ml_02', 'score_ml_03', 'score_ml_04']].max(axis=1)

    # Calculate the average for 'score_ml_ddg_bind' columns
    dfa['score_ml_ddg_bind_max'] = dfa[['score_ml_ddg_bind_01', 'score_ml_ddg_bind_02', 'score_ml_ddg_bind_03', 'score_ml_ddg_bind_04']].max(axis=1)


    significant_dataframe=dfa[['uniprot','variant','clinvar_signifiance_00','q_structure','score_ml','score_ml_ddg_bind','std_score_ml','std_score_ml_ddg_bind', 'score_ml_max' ,'score_ml_ddg_bind_max']]

    significant_dataframe.to_csv('significant_dataframe_2.csv', index=False)
    return significant_dataframe

#big_dataframe_to_significant_dataframe()

def scatterplot_marginals_quadrant_for_significant_dataframe(ddg_cutoff,dddg_cutoff):
    df = pd.read_csv('significant_dataframe.csv')
    df = df[df['clinvar_signifiance_00']!='VUS']
    df = df[df['score_ml'].notnull()]
    df.reset_index(inplace=True)
    
    # Create quadrants based on x and y values
    df["quadrant"] = ""
    df.loc[(df["score_ml"] < ddg_cutoff) & (df["score_ml_ddg_bind"] >= dddg_cutoff), "quadrant"] = "I"
    df.loc[(df["score_ml"] >= ddg_cutoff) & (df["score_ml_ddg_bind"] >= dddg_cutoff), "quadrant"] = "II"
    df.loc[(df["score_ml"] < ddg_cutoff) & (df["score_ml_ddg_bind"] < dddg_cutoff), "quadrant"] = "III"
    df.loc[(df["score_ml"] >= ddg_cutoff) & (df["score_ml_ddg_bind"] < dddg_cutoff), "quadrant"] = "IV"


    dfpath = df[df["clinvar_signifiance_00"]=='pathogenic']
    dfbenign = df[df["clinvar_signifiance_00"]=='benign']
    #dfVUS = df[df["clinvar_signifiance_00"]=='VUS']
    dfconflict = df[df["clinvar_signifiance_00"]=='conflict']

    all_I = df["quadrant"].value_counts()["I"]
    path_I = dfpath[dfpath["quadrant"]== "I"]
    benign_I = dfbenign[dfbenign["quadrant"]=="I"]
    #vUS_I = dfVUS[dfVUS["quadrant"]== "I"]
    conflict_I = dfconflict[dfconflict["quadrant"]== "I"]

    all_II = df["quadrant"].value_counts()["II"]
    path_II = dfpath[dfpath["quadrant"]=="II"]
    benign_II = dfbenign[dfbenign["quadrant"]== "II"]
    #vUS_II = dfVUS[dfVUS["quadrant"]== "II"]
    conflict_II = dfconflict[dfconflict["quadrant"]== "II"]

    all_III = df["quadrant"].value_counts()["III"]
    path_III = dfpath[dfpath["quadrant"]=="III"]
    benign_III = dfbenign[dfbenign["quadrant"]=="III"]
    #vUS_III = dfVUS[dfVUS["quadrant"]== "III"]
    conflict_III = dfconflict[dfconflict["quadrant"]== "III"]

    all_IV = df["quadrant"].value_counts()["IV"]
    path_IV = dfpath[dfpath["quadrant"]=="IV"]
    benign_IV = dfbenign[dfbenign["quadrant"]=="IV"]
    #vUS_IV = dfVUS[dfVUS["quadrant"]== "IV"]
    conflict_IV = dfconflict[dfconflict["quadrant"]== "IV"]

    all = len(df["quadrant"])
    path_all = len(dfpath)
    benign_all = len(dfbenign)
    #VUS_all = len(dfVUS)
    conflict_all = len(dfconflict)

    # Create the Seaborn jointplot
    #sns.color_palette("Spectral", as_cmap=True)
    #markers = {"Homodimer": "^", "Homotetramer": "v", "Heterodimer":">", "Heterotetramer":"<"}
    #sns.set_theme(style="darkgrid", palette="gist_ncar_r", font="Tahoma")
    palette = ["#AA61E2", "#FFD700","#F46049", "#40E0D0",  "#BADA55"]
    sns.set_palette(palette)
    sns.set_theme(style="darkgrid", palette=palette, font="Tahoma", )
    #g = sns.jointplot(data=df, x="score_ml", y="score_ml_ddg_bind", hue="q_structure", alpha=.5, space=.8 )
    g = sns.jointplot(data=df, x="score_ml", y="score_ml_ddg_bind", hue="clinvar_signifiance_00", alpha=.5, space=.8)
    g.refline(x=ddg_cutoff, y=dddg_cutoff, marginal=False)
    g.figure.subplots_adjust(top=.9)
    sns.move_legend(g.ax_joint, "upper center", bbox_to_anchor=(.5, 1.1), ncol=3, title="Clinical significance", fontsize=12)
    #sns.move_legend(g.ax_joint, "upper center", bbox_to_anchor=(.5, 1.1), ncol=4, title="Quaternary structure", fontsize=12)
    g.figure.suptitle(f'Scatterplot with marginals, showing the distribution \nof pathogenic and benign clinical variants')
    plt.text(-3.5,2, f"I \nDots: {all_I}/{all}={round((all_I/all)*100, 2)}%\nPathogenic: {len(path_I)}/{path_all}={round((len(path_I)/path_all)*100, 2)}%\nBenign: {len(benign_I)}/{benign_all}={round((len(benign_I)/benign_all)*100, 2)}%\nConflicting: {len(conflict_I)}/{conflict_all}={round((len(conflict_I)/conflict_all)*100, 2)}%", fontsize=11)
    plt.text(9.5,2, f"II \nDots: {all_II}/{all}={round((all_II/all)*100,2)}%\nPathogenic: {len(path_II)}/{path_all}={round((len(path_II)/path_all)*100,2)}%\nBenign: {len(benign_II)}/{benign_all}={round((len(benign_II)/benign_all)*100, 2)}%\nConflicting: {len(conflict_II)}/{conflict_all}={round((len(conflict_II)/conflict_all)*100, 2)}%",fontsize=11)
    plt.text(-3.5,-1, f"III \nDots: {all_III}/{all}={round((all_III/all)*100, 2)}%\nPathogenic: {len(path_III)}/{path_all}={round((len(path_III)/path_all)*100,2)}%\nBenign: {len(benign_III)}/{benign_all}={round((len(benign_III)/benign_all)*100,2)}%\nConflicting: {len(conflict_III)}/{conflict_all}={round((len(conflict_III)/conflict_all)*100, 2)}%",fontsize=11) 
    plt.text(9.5,-1, f"IV \nDots: {all_IV}/{all}={round((all_IV/all)*100,2)}%\nPathogenic: {len(path_IV)}/{path_all}={round((len(path_IV)/path_all)*100,2)}%\nBenign: {len(benign_IV)}/{benign_all}={round((len(benign_IV)/benign_all)*100,2)}%\nConflicting: {len(conflict_IV)}/{conflict_all}={round((len(conflict_IV)/conflict_all)*100, 2)}%",fontsize=11)
    #plt.text(-3.5,2, f"I \nDots: {all_I}/{all}={round(all_I/all, 3)}\nPathogenic: {len(path_I)}/{path_all}={round(len(path_I)/path_all, 3)}\nBenign: {len(benign_I)}/{benign_all}={round(len(benign_I)/benign_all, 3)}\nConflicting: {len(conflict_I)}/{conflict_all}={round(len(conflict_I)/conflict_all, 3)} \nVUS: {len(vUS_I)}/{VUS_all}={round(len(vUS_I)/VUS_all, 3)}")
    #plt.text(9.5,2, f"II \nDots: {all_II}/{all}={round(all_II/all,3)}\nPathogenic: {len(path_II)}/{path_all}={round(len(path_II)/path_all,3)}\nBenign: {len(benign_II)}/{benign_all}={round(len(benign_II)/benign_all, 3)}\nConflicting: {len(conflict_II)}/{conflict_all}={round(len(conflict_II)/conflict_all, 3)} \nVUS: {len(vUS_II)}/{VUS_all}={round(len(vUS_II)/VUS_all, 3)}")
    #plt.text(-3.5,-0.5, f"III \nDots: {all_III}/{all}={round(all_III/all, 3)}\nPathogenic: {len(path_III)}/{path_all}={round(len(path_III)/path_all,3)}\nBenign: {len(benign_III)}/{benign_all}={round(len(benign_III)/benign_all,3)}\nConflicting: {len(conflict_III)}/{conflict_all}={round(len(conflict_III)/conflict_all, 3)} \nVUS: {len(vUS_III)}/{VUS_all}={round(len(vUS_III)/VUS_all, 3)}")
    #plt.text(9.5,-0.5, f"IV \nDots: {all_IV}/{all}={round(all_IV/all,3)}\nPathogenic: {len(path_IV)}/{path_all}={round(len(path_IV)/path_all,3)}\nBenign: {len(benign_IV)}/{benign_all}={round(len(benign_IV)/benign_all,3)}\nConflicting: {len(conflict_IV)}/{conflict_all}={round(len(conflict_IV)/conflict_all, 3)} \nVUS: {len(vUS_IV)}/{VUS_all}={round(len(vUS_IV)/VUS_all, 3)}")
    
    path_I.to_csv('pathogenic_quadrant_I.csv', index=False)

    plt.show()


scatterplot_marginals_quadrant_for_significant_dataframe(2,0.5)

def max_plotted(ddg_cutoff,dddg_cutoff):
    df = pd.read_csv('significant_dataframe_2.csv')
    df = df[df['clinvar_signifiance_00']!='VUS']
    df = df[df['score_ml'].notnull()]
    df.reset_index(inplace=True)
    
    # Create quadrants based on x and y values
    df["quadrant"] = ""
    df.loc[(df["score_ml_max"] < ddg_cutoff) & (df["score_ml_ddg_bind_max"] >= dddg_cutoff), "quadrant"] = "I"
    df.loc[(df["score_ml_max"] >= ddg_cutoff) & (df["score_ml_ddg_bind_max"] >= dddg_cutoff), "quadrant"] = "II"
    df.loc[(df["score_ml_max"] < ddg_cutoff) & (df["score_ml_ddg_bind_max"] < dddg_cutoff), "quadrant"] = "III"
    df.loc[(df["score_ml_max"] >= ddg_cutoff) & (df["score_ml_ddg_bind_max"] < dddg_cutoff), "quadrant"] = "IV"


    dfpath = df[df["clinvar_signifiance_00"]=='pathogenic']
    dfbenign = df[df["clinvar_signifiance_00"]=='benign']
    #dfVUS = df[df["clinvar_signifiance_00"]=='VUS']
    dfconflict = df[df["clinvar_signifiance_00"]=='conflict']

    all_I = df["quadrant"].value_counts()["I"]
    path_I = dfpath[dfpath["quadrant"]== "I"]
    benign_I = dfbenign[dfbenign["quadrant"]=="I"]
    #vUS_I = dfVUS[dfVUS["quadrant"]== "I"]
    conflict_I = dfconflict[dfconflict["quadrant"]== "I"]

    all_II = df["quadrant"].value_counts()["II"]
    path_II = dfpath[dfpath["quadrant"]=="II"]
    benign_II = dfbenign[dfbenign["quadrant"]== "II"]
    #vUS_II = dfVUS[dfVUS["quadrant"]== "II"]
    conflict_II = dfconflict[dfconflict["quadrant"]== "II"]

    all_III = df["quadrant"].value_counts()["III"]
    path_III = dfpath[dfpath["quadrant"]=="III"]
    benign_III = dfbenign[dfbenign["quadrant"]=="III"]
    #vUS_III = dfVUS[dfVUS["quadrant"]== "III"]
    conflict_III = dfconflict[dfconflict["quadrant"]== "III"]

    all_IV = df["quadrant"].value_counts()["IV"]
    path_IV = dfpath[dfpath["quadrant"]=="IV"]
    benign_IV = dfbenign[dfbenign["quadrant"]=="IV"]
    #vUS_IV = dfVUS[dfVUS["quadrant"]== "IV"]
    conflict_IV = dfconflict[dfconflict["quadrant"]== "IV"]

    all = len(df["quadrant"])
    path_all = len(dfpath)
    benign_all = len(dfbenign)
    #VUS_all = len(dfVUS)
    conflict_all = len(dfconflict)

    # Create the Seaborn jointplot
    #sns.color_palette("Spectral", as_cmap=True)
    #markers = {"Homodimer": "^", "Homotetramer": "v", "Heterodimer":">", "Heterotetramer":"<"}
    #sns.set_theme(style="darkgrid", palette="gist_ncar_r", font="Tahoma")
    palette = ["#AA61E2", "#FFD700","#F46049", "#40E0D0",  "#BADA55"]
    sns.set_palette(palette)
    sns.set_theme(style="darkgrid", palette=palette, font="Tahoma")
    g = sns.jointplot(data=df, x="score_ml_max", y="score_ml_ddg_bind_max", hue="q_structure", alpha=.5, space=.8 )
    #g = sns.jointplot(data=df, x="score_ml_max", y="score_ml_ddg_bind_max", hue="clinvar_signifiance_00", alpha=.5, space=.8 )
    g.refline(x=ddg_cutoff, y=dddg_cutoff, marginal=False)
    g.figure.subplots_adjust(top=.9)
    #sns.move_legend(g.ax_joint, "upper center", bbox_to_anchor=(.5, 1.1), ncol=3, title="Clinical significance", fontsize=12)
    sns.move_legend(g.ax_joint, "upper center", bbox_to_anchor=(.5, 1.1), ncol=4, title="Quaternary structure", fontsize=12)
    g.figure.suptitle(f'Scatterplot with marginals, showing the distribution \nof pathogenic and benign clinical variants for their maximum RaSP values')
    plt.text(-4,3, f"I \nDots: {all_I}/{all}={round((all_I/all)*100, 2)}%\nPathogenic: {len(path_I)}/{path_all}={round((len(path_I)/path_all)*100, 2)}%\nBenign: {len(benign_I)}/{benign_all}={round((len(benign_I)/benign_all)*100, 2)}%\nConflicting: {len(conflict_I)}/{conflict_all}={round((len(conflict_I)/conflict_all)*100, 2)}%", fontsize=11)
    plt.text(9.5,3, f"II \nDots: {all_II}/{all}={round((all_II/all)*100,2)}%\nPathogenic: {len(path_II)}/{path_all}={round((len(path_II)/path_all)*100,2)}%\nBenign: {len(benign_II)}/{benign_all}={round((len(benign_II)/benign_all)*100, 2)}%\nConflicting: {len(conflict_II)}/{conflict_all}={round((len(conflict_II)/conflict_all)*100, 2)}%",fontsize=11)
    plt.text(-4,-1, f"III \nDots: {all_III}/{all}={round((all_III/all)*100, 2)}%\nPathogenic: {len(path_III)}/{path_all}={round((len(path_III)/path_all)*100,2)}%\nBenign: {len(benign_III)}/{benign_all}={round((len(benign_III)/benign_all)*100,2)}%\nConflicting: {len(conflict_III)}/{conflict_all}={round((len(conflict_III)/conflict_all)*100, 2)}%",fontsize=11) 
    plt.text(9.5,-1, f"IV \nDots: {all_IV}/{all}={round((all_IV/all)*100,2)}%\nPathogenic: {len(path_IV)}/{path_all}={round((len(path_IV)/path_all)*100,2)}%\nBenign: {len(benign_IV)}/{benign_all}={round((len(benign_IV)/benign_all)*100,2)}%\nConflicting: {len(conflict_IV)}/{conflict_all}={round((len(conflict_IV)/conflict_all)*100, 2)}%",fontsize=11)
    
    path_I.to_csv('pathogenic_quadrant_I_max.csv', index=False)

    plt.show()

#max_plotted(2,0.5)

