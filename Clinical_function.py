#Function to read the nicholas clinical data into scatterplots 
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy

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

# A function to spit out a list of the clinical pathogenic variants and their corresponding theoretical ddg or dddg
def clinical_var_csv_to_list_with_stability(variantcsvpath, uniprotnr, RaSPdatacsvpath, pdbid):
    df = pd.read_csv(RaSPdatacsvpath)
    df = df[df['pdbid'] == pdbid]

    # Empty lists for the pathogenic variants and benign variants
    path_variant = []
    path_ddg = []
    path_dddg = []
    benign_variant = []
    benign_ddg = []
    benign_dddg = []

    # Getting the clinvar variants as a list
    list1 = flatten(variant_and_related_disease(variantcsvpath, uniprotnr,"pathogenic"))
    list2 = flatten(variant_and_related_disease(variantcsvpath, uniprotnr, "benign"))

    for word in list1:
        # Here we use list comprehension on the dataframe to filter for just the entries with our variants. This is almost instant 
        filtered = df[df['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid:
                 path_variant.append(row['variant'])
                 path_ddg.append(row['score_ml'])
                 path_dddg.append(row['score_ml_ddg_bind'])
    
    for word in list2:
        filtered = df[df['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid:
                 benign_variant.append(row['variant'])
                 benign_ddg.append(row['score_ml'])
                 benign_dddg.append(row['score_ml_ddg_bind'])

    #create empty dictionaries for the pathogenic and the benign variants
    pathogenic_dict = {}
    benign_dict = {}

    # Creating the pathogenic dictionary
    for i in range(len(path_variant)):
        pathogenic_dict[path_variant[i]] = (path_ddg[i], path_dddg[i])
    # Creating the benign dictionary
    for i in range(len(benign_variant)):
        benign_dict[benign_variant[i]] = (benign_ddg[i], benign_dddg[i])
    
    pathkeys = []
    pathvalues = []
    benignkeys = []
    benignvalues = []
    for key, value in pathogenic_dict.items():
        if value[1] >= 0.1:
            pathkeys.append(key)
            pathvalues.append(value)

    for key, value in benign_dict.items():
        if value[1] >= 0.1:
            benignkeys.append(key)
            benignvalues.append(value)

    dddg01 = pd.DataFrame({'Pathogenic variants': pathkeys, 'Pathogenic variant values' : pathvalues})
    dddgbenign01 = pd.DataFrame({'Benign variants': benignkeys, 'Benign variant values' : benignvalues})
    return print(f"The pathogenic variants of {pdbid} with dddg of more than 0.1: {dddg01} \nThe benign variants of {pdbid} with dddg of more than 0.1: {dddgbenign01}")


# A function to create scatter plots showing the clinical variants, either pathogenic or benign, for given PDB IDs
def clinical_var_scatterplot_ddg_and_dddg_2(variantcsvpath1, uniprotnr1, RaSPdatacsvpath1, pdbid1, variantcsvpath2, uniprotnr2, RaSPdatacsvpath2, pdbid2):
    df1 = pd.read_csv(RaSPdatacsvpath1)
    df1 = df1[df1['pdbid'] == pdbid1]
    df2 = pd.read_csv(RaSPdatacsvpath2)
    df2 = df2[df2['pdbid'] == pdbid2]

    # Empty lists for the pathogenic variants and benign variants
    path_ddg1 = []
    path_dddg1 = []
    benign_ddg1 = []
    benign_dddg1 = []
    path_ddg2 = []
    path_dddg2 = []
    benign_ddg2 = []
    benign_dddg2 = []

    # Getting the clinvar variants as a list
    list1 = flatten(variant_and_related_disease(variantcsvpath1, uniprotnr1 ,"pathogenic"))
    list2 = flatten(variant_and_related_disease(variantcsvpath1, uniprotnr1 ,"benign"))
    list3 = flatten(variant_and_related_disease(variantcsvpath2, uniprotnr2 ,"pathogenic"))
    list4 = flatten(variant_and_related_disease(variantcsvpath2, uniprotnr2 ,"benign"))

    for word in list1:
        # Here we use list comprehension on the dataframe to filter for just the entries with our variants. This is almost instant 
        filtered = df1[df1['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid1:
                 path_ddg1.append(row['score_ml'])
                 path_dddg1.append(row['score_ml_ddg_bind'])
    
    for word in list2:
        filtered = df1[df1['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid1:
                 benign_ddg1.append(row['score_ml'])
                 benign_dddg1.append(row['score_ml_ddg_bind'])
    
    for word in list3:
        # Here we use list comprehension on the dataframe to filter for just the entries with our variants. This is almost instant 
        filtered = df2[df2['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid2:
                 path_ddg2.append(row['score_ml'])
                 path_dddg2.append(row['score_ml_ddg_bind'])
    
    for word in list4:
        filtered = df2[df2['variant']==word]
        for _, row in filtered.iterrows():
            if row['pdbid'] == pdbid2:
                 benign_ddg2.append(row['score_ml'])
                 benign_dddg2.append(row['score_ml_ddg_bind'])
    
    figure = plt.figure(figsize = [8,8])
    plt.scatter(path_ddg1, path_dddg1, label=f'Pathogenic {pdbid1} chain A', color='red',  alpha=0.7)
    plt.scatter(benign_ddg1, benign_dddg1, label=f'Benign {pdbid1} chain A', color='lightsalmon', alpha=0.7)
    plt.scatter(path_ddg2, path_dddg2, label=f'Pathogenic {pdbid2} chain B', color='darkorchid',  alpha=0.7)
    plt.scatter(benign_ddg2, benign_dddg2, label=f'Benign {pdbid2} chain B', color='orchid', alpha=0.7)
    plt.xlabel('ddg')
    plt.ylabel('dddg')
    plt.title(f'Scatterplot of selected simple clinvar variants for {pdbid1} and {pdbid2} ddg and dddg values')
    plt.legend()
    plt.grid(True, linestyle='-', alpha=0.5)
    plt.show()
    return figure

def clinical_var_scatterplot_ddg_and_dddg_1(variantcsvpath1, uniprotnr1, RaSPdatacsvpath1, pdbid1):
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

    # Create the Seaborn jointplot
    sns.set_theme(style="whitegrid")
    sns.jointplot(data=filtered_df, x="score_ml", y="score_ml_ddg_bind", hue="clinical_signi", kind="scatter")
    plt.show()

        ## Empty lists for the pathogenic variants and benign variants
    #path_ddg1 = []
    #path_dddg1 = []
    #benign_ddg1 = []
    #benign_dddg1 = []
    #patho = []
    #benign = []
    #patho_variant = []
    #benign_variant = []
#
    #for word in list1:
    #    # Here we use list comprehension on the dataframe to filter for just the entries with our variants. This is almost instant 
    #    filtered = df1[df1['variant']==word]
    #    for _, row in filtered.iterrows():
    #        if row['pdbid'] == pdbid1:
    #             patho_variant.append(row['variant'])
    #             path_ddg1.append(row['score_ml'])
    #             path_dddg1.append(row['score_ml_ddg_bind'])
    #             patho.append('Pathogenic')
    #
    #for word in list2:
    #    filtered = df1[df1['variant']==word]
    #    for _, row in filtered.iterrows():
    #        if row['pdbid'] == pdbid1:
    #             benign_variant.append(row['variant'])
    #             benign_ddg1.append(row['score_ml'])
    #             benign_dddg1.append(row['score_ml_ddg_bind'])
    #             benign.append('Benign')
    #
    #variants = [patho_variant, benign_variant]
    #ddglist = [path_ddg1, benign_ddg1]
    #dddglist = [path_dddg1, benign_dddg1]
    #clinsig = [patho, benign]
#
    #dict = {'Variant':variants,'Clinical_significance':clinsig, 'DDG':ddglist, 'DDDG':dddglist}
#
    #df2 = pd.DataFrame(dict)
#
    ##figure = plt.figure(figsize = [8,8])
    #sns.jointplot(df2, x='DDG', y='DDDG', hue='Clinical_significance', kind='scatter', palette='rainbow')
    ##sns.jointplot(x=df2['Benign ddg'], y=df2['Benign dddg'], kind='scatter', color='m', edgecolor="limegreen", linewidth=2)
    #
    ##plt.scatter(path_ddg1, path_dddg1, label=f'Pathogenic {pdbid1} chain A', color='red',  alpha=0.7)
    ##plt.scatter(benign_ddg1, benign_dddg1, label=f'Benign {pdbid1} chain A', color='lightsalmon', alpha=0.7)
    ##plt.xlabel('ddg')
    ##plt.ylabel('dddg')
    ##plt.title(f'Scatterplot of variants for {pdbid1} ddg and dddg values')
    ##plt.legend()
    ##plt.grid(True, linestyle='-', alpha=0.5)
    #plt.show()
    ##return figure


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
    