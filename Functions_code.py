#The code from my PUK
import seaborn as sns 
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.colors
from PrismData import PrismParser, VariantData
import numpy as np 
import os
from Bio import Align

''' A function to create heatmaps showing the change in thermodynamic stability (ddg or dddg) for a PDB ID
 across all theoretic mutations across all positions in the sequence'''
def csv_to_two_heatmaps(csvfilepath, pdbid, valuecolumn, chainA=True, chainB=True, chainC=True, chainD=True, show_A =True, show_B=True, show_C=True, show_D=True):
    # Read the CSV file into a pandas dataframe
    df = pd.read_csv(f'{csvfilepath}')
    df["pos"] = df["variant"].str[1:-1].astype("int")

    # Determine the aminoacid residues in the variant column of the csv file
    df["aa"] = df["variant"].str[-1]
    
    # Filter the dataframe to only include rows where the PDB id is the right one
    df = df[df['pdbid'] == pdbid]
    if chainA==True: 
        dfA = df[df['chainid'] == 'A'] 
        # Create a pivot table with the desired columns and index
        structure_A = pd.pivot_table(dfA.dropna(subset=['aa']), values=valuecolumn, index=["aa"], columns="pos" )
        structure_A = structure_A.reindex(dfA["aa"].unique())
    if chainB==True: 
        dfB = df[df['chainid'] == 'B']
        # Create a pivot table with the desired columns and index
        structure_B = pd.pivot_table(dfB.dropna(subset=['aa']), values=valuecolumn, index=["aa"], columns="pos" )
        structure_B = structure_B.reindex(dfB["aa"].unique())
    if chainC==True: 
        dfC = df[df['chainid'] == 'C']
        # Create a pivot table with the desired columns and index
        structure_C = pd.pivot_table(dfC.dropna(subset=['aa']), values=valuecolumn, index=["aa"], columns="pos" )
        structure_C = structure_C.reindex(dfC["aa"].unique())
    if chainD==True: 
        dfD = df[df['chainid'] == 'D']
        # Create a pivot table with the desired columns and index
        structure_D = pd.pivot_table(dfD.dropna(subset=['aa']), values=valuecolumn, index=["aa"], columns="pos" )
        structure_D = structure_D.reindex(dfD["aa"].unique())
    
    # Change the names for the value column so it is more understandable in the final title
    if valuecolumn == 'score_ml':
        valuecolumn = 'DDG'
    elif valuecolumn == 'score_ml_ddg_bind':
        valuecolumn = 'DDDG'

    # Create a figure with two subplots
    fig, axs = plt.subplots(1, 2, figsize=(16,8))

    if show_A==True:
        # Plot the heatmap for the PDB ID and chosen chains
        sns.heatmap(structure_A , vmin=-2 , vmax=12, cmap='rainbow',ax=axs[0], yticklabels='auto', cbar_kws={'orientation':'horizontal'}).set_title(f'{valuecolumn} for PDB ID {pdbid} chain A')
    if show_B==True:
        # Plot the heatmap for the PDB ID and chosen chains
        sns.heatmap(structure_B , vmin=-2 , vmax=12 , cmap='rainbow',ax=axs[1], yticklabels='auto', cbar_kws={'orientation':'horizontal'}).set_title(f'{valuecolumn} for PDB ID {pdbid} chain B')
    if show_C==True:
        # Plot the heatmap for the PDB ID and chosen chains
        sns.heatmap(structure_C , vmin=-2 , vmax=12, cmap='rainbow',ax=axs[0], yticklabels='auto', cbar_kws={'orientation':'horizontal'}).set_title(f'{valuecolumn} for PDB ID {pdbid} chain C')
    if show_D==True:
        # Plot the heatmap for the PDB ID and chosen chains
        sns.heatmap(structure_D , vmin=-2 , vmax=12 , cmap='rainbow',ax=axs[1], yticklabels='auto', cbar_kws={'orientation':'horizontal'}).set_title(f'{valuecolumn} for PDB ID {pdbid} chain D')
    
    #Save the figure
    #plt.savefig(f"/content/output/Heatmap_for_{pdbid}_{valuecolumn}.png")
    plt.show()
    return fig



''' A function to create a histogram of the number of variants with the different changes in stability'''
def csv_to_histogram(csvfilepath, pdbid, valuecolumn, minval, maxval, chainA=True, chainB=True, chainC=True, chainD=True):
    # Read the CSV file and get a dataframe with the wanted pdb and stability calculations
    df = pd.read_csv(csvfilepath)
    df = df[df['pdbid'] == pdbid] 
    df = df[df[valuecolumn].between(minval, maxval)]
    df["pos"] = df["variant"].str[1:-1].astype("int")
    if chainA == True:
        dfA = df[df['chainid'] == 'A'] 
        # Create the histogram for the chain
        plt.hist(dfA[valuecolumn], bins=30, color='yellow', alpha=0.7, label='Chain A')
    if chainB ==True:
        dfB = df[df['chainid'] == 'B']
        # Create the histogram for the chain
        plt.hist(dfB[valuecolumn], bins=30, color='orange', alpha=0.4, label='Chain B')
    if chainC ==True:
        dfC = df[df['chainid'] == 'C'] 
        # Create the histogram for the chain
        plt.hist(dfC[valuecolumn], bins=30, color='purple', alpha=0.4, label='Chain C')
    if chainD==True:
        dfD = df[df['chainid'] == 'D']
        # Create the histogram for the chain
        plt.hist(dfD[valuecolumn], bins=30, color='red', alpha=0.4, label='Chain D')

    # Change the names for the value column so it is more understandable in the final title
    if valuecolumn == 'score_ml':
        valuecolumn = 'DDG'
    elif valuecolumn == 'score_ml_ddg_bind':
        valuecolumn = 'DDDG'

    # Set the title and axis labels
    plt.xlabel(f'{valuecolumn} values for {pdbid}')
    plt.ylabel('Frequency')
    plt.legend()

    return plt.show()


'''functions used to get the utilized sequence for the protein or protein sequence'''
def every_nth(list, nth):
    # Use list slicing to return elements starting from the (nth-1) index, with a step of 'nth'.
    return list[1::nth]
def sequence_from_csv(csvfilepath, pdbid, chain): 
    ## Read the CSV file into a pandas dataframe
    #df = pd.read_csv(csvfilepath)
    #df = df[df['pdbid'] == pdbid]
    #df = df[df['chainid']==chain]
    #df = df["variant"].str[0:1].astype("str")
    #aa_sequence = every_nth(df, 19)
    #aa_list = ''.join(aa_sequence)
    ##pd.set_option('display.max_rows', None)

    # Read the CSV file into a pandas dataframe
    df = pd.read_csv(csvfilepath)
    df = df[df['pdbid'] == pdbid]
    df = df[df['chainid']==chain]
    df = df['variant'].astype("str")
    aa_sequence = every_nth(df, 20)
    aa_seq_m_x =[]

    for index, val in enumerate(aa_sequence):
        num = val[1:-1]
        if index <= 2: aa_seq_m_x.append(val)
        elif str(int(num)-1) in aa_seq_m_x[-1]:
            aa_seq_m_x.append(val)
        elif str(int(num)-1) not in aa_seq_m_x[-1:-2]:
            m = int(num) - int(aa_seq_m_x[-1][1:-1])
            #mx = m+1
            for j in range(m):
                aa_seq_m_x.append(f'X{int(j)+int(num)}X')
    
    
    aa_list = ''.join([x[0] for x in aa_seq_m_x])
    #for n in aa_seq_m_x:
    #    aminoacid = n[0:1]
    #    aa.append[aminoacid]

    #df = df["variant"].str[0:1].astype("str")
    #aa_sequence = every_nth(df, 19)
    #aa_list = ''.join(aa)
    #pd.set_option('display.max_rows', None)

    return aa_list

''' This function flattens a nested list like [[0], [1], [2], [3], [4]] into just [0, 1, 2, 3, 4] '''
def flatten(list):
    return [x for xs in list for x in xs]

''' Converting the clinvar variant data into the same format as the RaSP csv file data'''
def variant_and_related_disease(clinvarcsvfile):
    #read the CSV file into a dataframe
    df = pd.read_csv(clinvarcsvfile)
    df = df[[ "ref_aa", "pos_aa", "alt_aa", "ClinicalSignificance", "PhenotypeList",]]
    df = df.sort_values("pos_aa")

    #replace all the three letter codes for aminoacids with one letter code
    df.replace("Ala", "A", inplace=True)
    df.replace("Arg", "R", inplace=True)
    df.replace("Asn", "N", inplace=True)
    df.replace("Asp", "D", inplace=True)
    df.replace("Cys", "C", inplace=True)
    df.replace("Glu", "E", inplace=True)
    df.replace("Gln", "Q", inplace=True)
    df.replace("Gly", "G", inplace=True)
    df.replace("His", "H", inplace=True)
    df.replace("Ile", "I", inplace=True)
    df.replace("Leu", "L", inplace=True)
    df.replace("Lys", "K", inplace=True) 
    df.replace("Met", "M", inplace=True) 
    df.replace("Phe", "F", inplace=True)
    df.replace("Pro", "P", inplace=True)
    df.replace("Ser", "S", inplace=True) 
    df.replace("Thr", "T", inplace=True) 
    df.replace("Trp", "W", inplace=True) 
    df.replace("Tyr", "Y", inplace=True) 
    df.replace("Val", "V", inplace=True) 
    
    #create a column for the variant
    df["variant"] = df["ref_aa"] + df["pos_aa"].astype(str) +df["alt_aa"]
    # If clinical significance and phenotype information is wanted as well, the hashtag can be removed from the next line
    #df = df[[ "variant" , "ClinicalSignificance", "PhenotypeList"]]
    df = df[[ "variant"]]
    list_of_data = df.values.tolist()
    return list_of_data

# A function to spit out a list of the clinical pathogenic variants and their corresponding theoretical ddg or dddg
def clinvar_csv_to_list_with_stability(clinvarpathogeniccsvpath, clinvarbenigncsvpath, RaSPdatacsvpath, pdbid):
    #read csv file into dataframe
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
    list1 = flatten(variant_and_related_disease(clinvarpathogeniccsvpath))
    list2 = flatten(variant_and_related_disease(clinvarbenigncsvpath))


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
    
    #empty lists for the keys and values of the dictionaries
    pathkeys = []
    pathvalues = []
    benignkeys = []
    benignvalues = []

    #appending keys and values to their rightful list, from the benign and pathogenic lists 
    for key, value in pathogenic_dict.items():
        if value[1] >= 0.1:
            pathkeys.append(key)
            pathvalues.append(value)

    for key, value in benign_dict.items():
        if value[1] >= 0.1:
            benignkeys.append(key)
            benignvalues.append(value)

    #create a new dataframe based on the lists 
    dddg01 = pd.DataFrame({'Pathogenic variants': pathkeys, 'Pathogenic variant values' : pathvalues})
    dddgbenign01 = pd.DataFrame({'Benign variants': benignkeys, 'Benign variant values' : benignvalues})
    return print(f"The pathogenic variants of {pdbid} with dddg of more than 0.1: {dddg01} \nThe benign variants of {pdbid} with dddg of more than 0.1: {dddgbenign01}")


# A function to create scatter plots showing the clinical variants, either pathogenic or benign, for given PDB IDs
def clinvar_scatterplot_ddg_and_dddg(clinvarpathogeniccsvpath1, clinvarbenigncsvpath1, RaSPdatacsvpath1, pdbid1, clinvarpathogeniccsvpath2, clinvarbenigncsvpath2, RaSPdatacsvpath2, pdbid2):
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
    list1 = flatten(variant_and_related_disease(clinvarpathogeniccsvpath1))
    list2 = flatten(variant_and_related_disease(clinvarbenigncsvpath1))
    list3 = flatten(variant_and_related_disease(clinvarpathogeniccsvpath2))
    list4 = flatten(variant_and_related_disease(clinvarbenigncsvpath2))

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
    #plt.scatter(path_ddg1, path_dddg1, label=f'Pathogenic {pdbid1} chain A', color='red',  alpha=0.7)
    #plt.scatter(benign_ddg1, benign_dddg1, label=f'Benign {pdbid1} chain A', color='lightsalmon', alpha=0.7)
    #plt.scatter(path_ddg2, path_dddg2, label=f'Pathogenic {pdbid2} chain B', color='darkorchid',  alpha=0.7)
    #plt.scatter(benign_ddg2, benign_dddg2, label=f'Benign {pdbid2} chain B', color='orchid', alpha=0.7)
    sns.jointplot(x=path_ddg1, y=path_dddg1, kind='scatter', color='m', edgecolor="skyblue", linewidth=2).set_title(f'Pathogenic variants of {pdbid1} chain A')
    sns.jointplot(x=benign_ddg1, y=benign_dddg1, kind='scatter', color='m', edgecolor="limegreen", linewidth=2).set_title(f'Benign variants of {pdbid1} chain A')
    sns.jointplot(x=path_ddg2, y=path_dddg2, kind='scatter', color='m', edgecolor="red", linewidth=2).set_title(f'Pathogenic variants of {pdbid2} chain B')
    sns.jointplot(x=benign_ddg2, y=benign_dddg2, kind='scatter', color='m', edgecolor="yellow", linewidth=2).set_title(f'Benign variants of {pdbid2} chain B')
    plt.xlabel('ddg')
    plt.ylabel('dddg')
    plt.title(f'Scatterplot of selected simple clinvar variants for {pdbid1} and {pdbid2} ddg and dddg values')
    plt.legend()
    plt.grid(True, linestyle='-', alpha=0.5)
    plt.show()
    return figure

def scatterplotting_old_and_new_predictions(oldcsv, oldpdbid, newcsv, newpdbid):
    #read the old and new dataframes into their own respective dataframes
    dfgammel = pd.read_csv(oldcsv)
    dfgammel = dfgammel[dfgammel['pdbid'] == oldpdbid]
    dfny = pd.read_csv(newcsv)
    dfny = dfny[dfny['pdbid'] == newpdbid]

    #create empty lists for the stability values to be found in each dataframe
    ddggammel = []
    ddgny = []
    dddggammel = []
    dddgny = []

    #appending all the values from the old and new dataframes, into the lists, in the rightfull order
    for row in dfgammel['score_ml']:
        ddggammel.append(row)
    
    for row in dfgammel['score_ml_ddg_bind']:
        dddggammel.append(row)
    
    for row in dfny['score_ml']:
        ddgny.append(row)
    
    for row in dfny['score_ml_ddg_bind']:
        dddgny.append(row)
    
    #Plot the scatterplot
    figure = plt.figure(figsize = [8,8])
    plt.scatter(ddggammel, dddggammel, label=f'Original RaSP data for {oldpdbid}', color='red',  alpha=0.7)
    plt.scatter(ddgny, dddgny, label=f'New RaSP data for {newpdbid}', color='purple', alpha=0.7)
    plt.xlabel('ddg')
    plt.ylabel('dddg')
    plt.title(f'Scatterplot of original {oldpdbid} RaSP data and new {newpdbid} RaSP data')
    plt.legend()
    plt.grid(True, linestyle='-', alpha=0.5)
    plt.show()
    return figure


def scatterplotting_old_vs_new_ddg_predictions(gammelcsv, gammelpdbid, nycsv, nypdbid):
    dfgammel = pd.read_csv(gammelcsv)
    dfgammel = dfgammel[dfgammel['pdbid'] == gammelpdbid]
    dfny = pd.read_csv(nycsv)
    dfny = dfny[dfny['pdbid'] == nypdbid]

    # Step 1: Filter rows with matching 'variant' and 'chainid'
    matching_rows = dfgammel.merge(dfny, on=['variant', 'chainid'], suffixes=('_old', '_new'))

    # Step 2: Extract desired columns
    score_ml_gammel = matching_rows['score_ml_old']
    score_ml_ny = matching_rows['score_ml_new']

    # Step 3: Create scatterplot using Seaborn
    cmap_blended = sns.blend_palette(["#AA61E2", "#40E0D0", "#BADA55", "#FFD700","#F46049" ], as_cmap=True)
    sns.set_theme(style="darkgrid", font="Tahoma")
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=score_ml_gammel, y=score_ml_ny, palette=cmap_blended, hue=matching_rows['score_ml_old'], legend=None)

    # Step 4: Fit a linear regression model using numpy
    slope, intercept = np.polyfit(score_ml_gammel, score_ml_ny, 1)
    trendline_x = np.linspace(min(score_ml_gammel), max(score_ml_gammel), 100)
    trendline_y = slope * trendline_x + intercept
    plt.plot(trendline_x, trendline_y, color='red', label='Trendline')

    # Step 5: Display R-value (calculate it separately)
    r_value = np.corrcoef(score_ml_gammel, score_ml_ny)[0, 1]
    plt.text(0.7, 0.1, f'R = {r_value:.2f}', transform=plt.gca().transAxes, fontsize=17)

    plt.xlabel(f'RaSP-PC data', fontsize=17)
    plt.ylabel(f'RaSP-PCAV data', fontsize=17)
    #plt.title(f'Old {gammelpdbid} DDG values against new {nypdbid} DDG values', fontsize=18)
    #plt.legend()
    plt.show()



def scatterplotting_old_vs_new_dddg_predictions(gammelcsv, gammelpdbid, nycsv, nypdbid):
    dfgammel = pd.read_csv(gammelcsv)
    dfgammel = dfgammel[dfgammel['pdbid'] == gammelpdbid]
    dfny = pd.read_csv(nycsv)
    dfny = dfny[dfny['pdbid'] == nypdbid]

    # Step 1: Filter rows with matching 'variant' and 'chainid'
    matching_rows = dfgammel.merge(dfny, on=['variant', 'chainid'], suffixes=('_old', '_new'))

    # Step 2: Extract desired columns
    score_ml_ddg_bind_gammel = matching_rows['score_ml_ddg_bind_old']
    score_ml_ddg_bind_ny = matching_rows['score_ml_ddg_bind_new']

    

    # Step 3: Create scatterplot using Seaborn
    cmap_blended = sns.blend_palette(["#AA61E2", "#40E0D0", "#BADA55", "#FFD700","#F46049" ], as_cmap=True)
    sns.set_theme(style="darkgrid", font="Tahoma")
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=score_ml_ddg_bind_gammel, y=score_ml_ddg_bind_ny, palette=cmap_blended, hue=matching_rows['score_ml_ddg_bind_old'], legend=None)

    # Step 4: Fit a linear regression model using numpy
    slope, intercept = np.polyfit(score_ml_ddg_bind_gammel, score_ml_ddg_bind_ny, 1)
    trendline_x = np.linspace(min(score_ml_ddg_bind_gammel), max(score_ml_ddg_bind_gammel), 100)
    trendline_y = slope * trendline_x + intercept
    plt.plot(trendline_x, trendline_y, color='red', label='Trendline')

    # Step 5: Display R-value (calculate it separately)
    r_value = np.corrcoef(score_ml_ddg_bind_gammel, score_ml_ddg_bind_ny)[0, 1]
    plt.text(0.7, 0.1, f'R = {r_value:.2f}', transform=plt.gca().transAxes, fontsize=17)

    plt.xlabel(f'RaSP-PC data ', fontsize=17)
    plt.ylabel(f'RaSP-PCAV data ', fontsize=17)
    #plt.title(f'Old {gammelpdbid} DDDG values against new {nypdbid} DDDG values', fontsize=18)
    #plt.legend()
    plt.show()    
    

def gathering_all_csv_files_into_one(folderpath,newname):
    #Get all files from a specific folder into one list
    all_files = os.listdir(folderpath)
    
    # Filter out non-CSV files from the read folder
    csv_files = [f for f in all_files if f.endswith('.csv')]
    
    # Create a list to hold the dataframes
    df_list = []
    
    for csv in csv_files:
        file_path = os.path.join(folderpath, csv)
        try:
            # Try reading the file using default UTF-8 encoding
            df = pd.read_csv(file_path)
            df_list.append(df)
        except UnicodeDecodeError:
            try:
                # reading the file using UTF-16 encoding with tab separator
                df = pd.read_csv(file_path, sep='\t', encoding='utf-16')
                df_list.append(df)
            except Exception as e:
                print(f"Could not read file {csv} because of error: {e}")
        except Exception as e:
            print(f"Could not read file {csv} because of error: {e}")
    
    # Concatenate all data into one DataFrame
    big_df = pd.concat(df_list, ignore_index=True)

    # Drop duplicate variants (keep only the first occurrence)
    big_df.drop_duplicates(subset=('pdbid','chainid','variant'), inplace=True)
    
    # Save the final result to a new CSV file
    big_df.to_csv(os.path.join(folderpath, f'{newname}.csv'), index=False)

def remove_wildtype_from_csv(csvfile, name):
    df = pd.read_csv(csvfile)
    # Extract the first and last characters from the "variant" column
    df["first_char"] = df["variant"].str[0]
    df["last_char"] = df["variant"].str[-1]

    # Filter out rows where the first and last characters are identical
    df = df[df["first_char"] != df["last_char"]]

    # Drop the temporary columns
    df.drop(columns=["first_char", "last_char"], inplace=True)

    # Save the cleaned dataframe to a new CSV file
    df.to_csv(f"clean_{name}.csv", index=False)

def sequence_from_nicolas(csvfile, uniprotnr):
    #read the csv file into a dataframe
    df = pd.read_csv(csvfile)
    df = df[df['uniprot']== uniprotnr]
    #get the second column of the file since it contains the entire sequence
    sequence = df['sequence'].iloc[2]
    return sequence 


def write_prismfile_from_RaSP_csv():
    #read csv into dataframe, and pick the chosen columns of the dataframe
    df = pd.read_csv(f"predictions\All_new_df_ml.csv" )
    dfb = df[df['pdbid']=='2C9S']
    dfb = dfb[dfb['chainid']=='F']
    dfb.reset_index(inplace=True)
    dfb = dfb[['variant', 'score_ml', 'score_ml_ddg_bind']]

    #write the meta data for the prism file
    md = {'version':1,
                'protein':{'name':'Human SUPEROXIDE DISMUTASE', 'organism':'Homo sapiens', 'uniprot':'P00441','first_residue_number':1, 'sequence':'ATKAVCVLKGDGPVQGIINFEQKESNGPVKVWGSIKGLTEGLHGFHVHEFGDNTAGCTSAGPHFNPLSRKHGGPKDEERHVGDLGNVTADKDGVADVSIEDSVISLSGDHCIIGRTLVVHEKADDLGKGGNEESTKTGNAGSRLACGVIGIAQ'},
                'method':{'name':'RaSP from a PDB'},
                'columns':{'score_ml':'DDG', 'score_ml_ddg_bind': 'DDDG'}}

    #Call the use of the VariantData module in prismdata 
    prismdata=VariantData(metadata=md, dataframe=dfb)
    prismdata.check()

    #parse and align the data
    pp = PrismParser()
    #write the prismfile into a directory with prism, the method behind the data and the uniprot ID 
    pp.write("Prismfiles/prism_RaSP_P00441_F.txt", prismdata)



def write_prismfile_from_Nicolas_csv():
    #read csv into dataframe, and pick the chosen columns of the dataframe
    df = pd.read_csv(f'nr2_Nicolas_data.csv')
    dfa = df[df['uniprot'] =='P16410']
    dfa.reset_index(inplace=True)
    dfa = dfa[["variant", "clinvar_signifiance"]]

    #write the meta data for the prism file
    md = {'version':1,
                'protein':{'name':'Human CYTOTOXIC T-LYMPHOCYTE-ASSOCIATED PROTEIN 4', 'organism':'Homo sapiens', 'uniprot':'P16410','first_residue_number':1, 'sequence':'MACLGFQRHKAQLNLATRTWPCTLLFFLLFIPVFCKAMHVAQPAVVLASSRGIASFVCEYASPGKATEVRVTVLRQADSQVTEVCAATYMMGNELTFLDDSICTGTSSGNQVNLTIQGLRAMDTGLYICKVELMYPPPYYLGIGNGTQIYVIDPEPCPDSDFLLWILAAVSSGLFFYSFLLTAVSLSKMLKKRSPLTTGVYVKMPPTEPECEKQFQPYFIPIN'},
                'method':{'name':'RaSP from a PDB'},
                'columns':{'clinvar_signifiance':'Clinical_significance'}}

    #Call the use of the VariantData module in prismdata
    prismdata=VariantData(metadata=md, dataframe=dfa)
    prismdata.check()

    #parse and align the data
    pp = PrismParser()
    #write the prismfile into a directory with prism, the method behind the data and the uniprot ID
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

    #concatenate the data into one dataframe
    big_dataframe = pd.concat([dataframe1,dataframe2,dataframe3,dataframe4,dataframe5,dataframe6,dataframe7,dataframe8,dataframe9,dataframe10,dataframe11,dataframe12,dataframe13,dataframe14], ignore_index=True)
    #get the dataframe as a csv file
    big_dataframe.to_csv('big_dataframe.csv', index=False)

    return big_dataframe


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

    #make a new dataframe based on the columns that explain what we want to show
    significant_dataframe=dfa[['uniprot','variant','clinvar_signifiance_00','q_structure','score_ml','score_ml_ddg_bind','std_score_ml','std_score_ml_ddg_bind', 'score_ml_max' ,'score_ml_ddg_bind_max']]

    #save the new dataframe to a csv file
    significant_dataframe.to_csv('significant_dataframe_2.csv', index=False)
    return significant_dataframe


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

    #make specific clinical significances into their own variable
    dfpath = df[df["clinvar_signifiance_00"]=='pathogenic']
    dfbenign = df[df["clinvar_signifiance_00"]=='benign']
    #dfVUS = df[df["clinvar_signifiance_00"]=='VUS']
    dfconflict = df[df["clinvar_signifiance_00"]=='conflict']

    #get the values from the rows in every quadrant
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
    #g = sns.jointplot(data=df, x="score_ml", y="score_ml_ddg_bind", hue="clinvar_signifiance_00", alpha=.5, space=.8)
    g = sns.jointplot(data=df, x="score_ml", y="score_ml_ddg_bind", hue="q_structure", alpha=.5, space=.8 )
    g.refline(x=ddg_cutoff, y=dddg_cutoff, marginal=False)
    g.figure.subplots_adjust(top=.9)
    #sns.move_legend(g.ax_joint, "upper center", bbox_to_anchor=(.5, 1.1), ncol=3, title="Clinical significance", fontsize=12)
    sns.move_legend(g.ax_joint, "upper center", bbox_to_anchor=(.5, 1.1), ncol=4, title="Quaternary structure", fontsize=12)
    #g.figure.suptitle(f'Figure 11: Scatterplot with marginals, showing the distribution of clinical variants\n across mean \u0394\u0394G and mean \u0394\u0394\u0394G values')
    g.figure.suptitle(f'Figure 12: Scatterplot with marginals, showing the distribution of clinical variants´ quaternary structure\n across mean \u0394\u0394G and mean \u0394\u0394\u0394G values')
    plt.text(-4,2, f"I \nDots: {all_I}/{all}={round((all_I/all)*100, 2)}%\nPathogenic: {len(path_I)}/{path_all}={round((len(path_I)/path_all)*100, 2)}%\nBenign: {len(benign_I)}/{benign_all}={round((len(benign_I)/benign_all)*100, 2)}%\nConflicting: {len(conflict_I)}/{conflict_all}={round((len(conflict_I)/conflict_all)*100, 2)}%", fontsize=11)
    plt.text(9.5,2, f"II \nDots: {all_II}/{all}={round((all_II/all)*100,2)}%\nPathogenic: {len(path_II)}/{path_all}={round((len(path_II)/path_all)*100,2)}%\nBenign: {len(benign_II)}/{benign_all}={round((len(benign_II)/benign_all)*100, 2)}%\nConflicting: {len(conflict_II)}/{conflict_all}={round((len(conflict_II)/conflict_all)*100, 2)}%",fontsize=11)
    plt.text(-4,-1, f"III \nDots: {all_III}/{all}={round((all_III/all)*100, 2)}%\nPathogenic: {len(path_III)}/{path_all}={round((len(path_III)/path_all)*100,2)}%\nBenign: {len(benign_III)}/{benign_all}={round((len(benign_III)/benign_all)*100,2)}%\nConflicting: {len(conflict_III)}/{conflict_all}={round((len(conflict_III)/conflict_all)*100, 2)}%",fontsize=11) 
    plt.text(9.5,-1, f"IV \nDots: {all_IV}/{all}={round((all_IV/all)*100,2)}%\nPathogenic: {len(path_IV)}/{path_all}={round((len(path_IV)/path_all)*100,2)}%\nBenign: {len(benign_IV)}/{benign_all}={round((len(benign_IV)/benign_all)*100,2)}%\nConflicting: {len(conflict_IV)}/{conflict_all}={round((len(conflict_IV)/conflict_all)*100, 2)}%",fontsize=11)
    #plt.text(-3.5,2, f"I \nDots: {all_I}/{all}={round(all_I/all, 3)}\nPathogenic: {len(path_I)}/{path_all}={round(len(path_I)/path_all, 3)}\nBenign: {len(benign_I)}/{benign_all}={round(len(benign_I)/benign_all, 3)}\nConflicting: {len(conflict_I)}/{conflict_all}={round(len(conflict_I)/conflict_all, 3)} \nVUS: {len(vUS_I)}/{VUS_all}={round(len(vUS_I)/VUS_all, 3)}")
    #plt.text(9.5,2, f"II \nDots: {all_II}/{all}={round(all_II/all,3)}\nPathogenic: {len(path_II)}/{path_all}={round(len(path_II)/path_all,3)}\nBenign: {len(benign_II)}/{benign_all}={round(len(benign_II)/benign_all, 3)}\nConflicting: {len(conflict_II)}/{conflict_all}={round(len(conflict_II)/conflict_all, 3)} \nVUS: {len(vUS_II)}/{VUS_all}={round(len(vUS_II)/VUS_all, 3)}")
    #plt.text(-3.5,-0.5, f"III \nDots: {all_III}/{all}={round(all_III/all, 3)}\nPathogenic: {len(path_III)}/{path_all}={round(len(path_III)/path_all,3)}\nBenign: {len(benign_III)}/{benign_all}={round(len(benign_III)/benign_all,3)}\nConflicting: {len(conflict_III)}/{conflict_all}={round(len(conflict_III)/conflict_all, 3)} \nVUS: {len(vUS_III)}/{VUS_all}={round(len(vUS_III)/VUS_all, 3)}")
    #plt.text(9.5,-0.5, f"IV \nDots: {all_IV}/{all}={round(all_IV/all,3)}\nPathogenic: {len(path_IV)}/{path_all}={round(len(path_IV)/path_all,3)}\nBenign: {len(benign_IV)}/{benign_all}={round(len(benign_IV)/benign_all,3)}\nConflicting: {len(conflict_IV)}/{conflict_all}={round(len(conflict_IV)/conflict_all, 3)} \nVUS: {len(vUS_IV)}/{VUS_all}={round(len(vUS_IV)/VUS_all, 3)}")
    plt.xlabel(f'\u0394\u0394G ', fontsize=17)
    plt.ylabel(f'\u0394\u0394\u0394G ', fontsize=17)
    path_I.to_csv('pathogenic_quadrant_I.csv', index=False)

    plt.show()


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

    #make specific clinical significances into their own variable
    dfpath = df[df["clinvar_signifiance_00"]=='pathogenic']
    dfbenign = df[df["clinvar_signifiance_00"]=='benign']
    #dfVUS = df[df["clinvar_signifiance_00"]=='VUS']
    dfconflict = df[df["clinvar_signifiance_00"]=='conflict']

    #get the values from the rows in every quadrant
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
    palette = ["#AA61E2", "#FFD700","#F46049", "#40E0D0",  "#BADA55"]
    sns.set_palette(palette)
    sns.set_theme(style="darkgrid", palette=palette, font="Tahoma")
    g = sns.jointplot(data=df, x="score_ml_max", y="score_ml_ddg_bind_max", hue="clinvar_signifiance_00", alpha=.5, space=.8 )
    #g = sns.jointplot(data=df, x="score_ml_max", y="score_ml_ddg_bind_max", hue="q_structure", alpha=.5, space=.8 )
    g.refline(x=ddg_cutoff, y=dddg_cutoff, marginal=False)
    g.figure.subplots_adjust(top=.9)
    sns.move_legend(g.ax_joint, "upper center", bbox_to_anchor=(.5, 1.1), ncol=3, title="Clinical significance", fontsize=12)
    #sns.move_legend(g.ax_joint, "upper center", bbox_to_anchor=(.5, 1.1), ncol=4, title="Quaternary structure", fontsize=12)
    g.figure.suptitle(f'Figure 13: Scatterplot with marginals, showing the distribution of clinical variants\n across max \u0394\u0394G and max \u0394\u0394\u0394G values')
    #g.figure.suptitle(f'Figure 14: Scatterplot with marginals, showing the distribution of clinical variants´ quaternary structure\n across max \u0394\u0394G and max \u0394\u0394\u0394G values')
    plt.text(-3.5,3, f"I \nDots: {all_I}/{all}={round((all_I/all)*100, 2)}%\nPathogenic: {len(path_I)}/{path_all}={round((len(path_I)/path_all)*100, 2)}%\nBenign: {len(benign_I)}/{benign_all}={round((len(benign_I)/benign_all)*100, 2)}%\nConflicting: {len(conflict_I)}/{conflict_all}={round((len(conflict_I)/conflict_all)*100, 2)}%", fontsize=11)
    plt.text(9.5,3, f"II \nDots: {all_II}/{all}={round((all_II/all)*100,2)}%\nPathogenic: {len(path_II)}/{path_all}={round((len(path_II)/path_all)*100,2)}%\nBenign: {len(benign_II)}/{benign_all}={round((len(benign_II)/benign_all)*100, 2)}%\nConflicting: {len(conflict_II)}/{conflict_all}={round((len(conflict_II)/conflict_all)*100, 2)}%",fontsize=11)
    plt.text(-3.5,-1, f"III \nDots: {all_III}/{all}={round((all_III/all)*100, 2)}%\nPathogenic: {len(path_III)}/{path_all}={round((len(path_III)/path_all)*100,2)}%\nBenign: {len(benign_III)}/{benign_all}={round((len(benign_III)/benign_all)*100,2)}%\nConflicting: {len(conflict_III)}/{conflict_all}={round((len(conflict_III)/conflict_all)*100, 2)}%",fontsize=11) 
    plt.text(9.5,-1, f"IV \nDots: {all_IV}/{all}={round((all_IV/all)*100,2)}%\nPathogenic: {len(path_IV)}/{path_all}={round((len(path_IV)/path_all)*100,2)}%\nBenign: {len(benign_IV)}/{benign_all}={round((len(benign_IV)/benign_all)*100,2)}%\nConflicting: {len(conflict_IV)}/{conflict_all}={round((len(conflict_IV)/conflict_all)*100, 2)}%",fontsize=11)
    plt.xlabel(f'\u0394\u0394G ', fontsize=17)
    plt.ylabel(f'\u0394\u0394\u0394G ', fontsize=17)

    path_I.to_csv('pathogenic_quadrant_I_max.csv', index=False)

    plt.show()


def get_p_and_g_origo(gammelcsv, gammelpdbid, nycsv, nypdbid, nameofnewfile):
    #Read the old and new csv files into their own dataframes
    dfgammel = pd.read_csv(gammelcsv)
    dfgammel = dfgammel[dfgammel['pdbid'] == gammelpdbid]
    dfny = pd.read_csv(nycsv)
    dfny = dfny[dfny['pdbid'] == nypdbid]

    #Filter rows with matching 'variant' and 'chainid'
    matching_rows = dfgammel.merge(dfny, on=['variant', 'chainid'], suffixes=('_old', '_new'))

    #get only the values around origo
    old_origo=matching_rows[matching_rows["score_ml_ddg_bind_old"].between(-0.2, 0.2)]
    new_origo=matching_rows[matching_rows["score_ml_ddg_bind_new"].between(-0.2, 0.2)]

    #collect the two dataframes with data around origo and save the dataframe as a csv file
    origo_data = pd.concat([old_origo, new_origo], ignore_index=True).drop_duplicates()
    origo_data.to_csv('dddg_origo_data.csv', index=False)

    #read the newly made csv file into it's own dataframe
    origo = pd.read_csv('dddg_origo_data.csv')

    #pick out the mutation in the variants column, and the wild type amino acid in the variants column
    origo["mut"] = origo["variant"].str[-1]
    origo["wt"] = origo["variant"].str[0]
    #Choose data from: the mutations that had proline, and choose wildtypes that had glycine 
    proline = origo[origo["mut"]=='P']
    glycine = origo[origo["wt"]=='G']

    #concatenate the proline and glycine dataframes and save them to their own csv file
    pg = pd.concat([proline, glycine], ignore_index=True)
    pg.to_csv(f'{nameofnewfile}.csv', index=False)


 



