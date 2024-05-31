#The code from my PUK
import seaborn as sns 
import pandas as pd 
import matplotlib.pyplot as plt
import matplotlib.colors
import numpy as np 
import os
from Bio import Align

# A function to create heatmaps showing the change in thermodynamic stability (ddg or dddg) for a PDB ID
# across all theoretic mutations across all positions in the sequence
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
    #plt.savefig(f"/content/output/Heatmap_for_{pdbid}_{valuecolumn}.png")
    plt.show()
    return fig



# A function to create a histogram of the number of variants with the different changes in stability
def csv_to_histogram(csvfilepath, pdbid, valuecolumn, minval, maxval, chainA=True, chainB=True, chainC=True, chainD=True):
    # Read the CSV file
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

# functions used to get the utilized sequence for the protein or protein sequence

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

# This function flattens a nested list like [[0], [1], [2], [3], [4]] into just [0, 1, 2, 3, 4] 
def flatten(list):
    return [x for xs in list for x in xs]

# Converting the clinvar variant data into the same format as the RaSP csv file data
def variant_and_related_disease(clinvarcsvfile):
    df = pd.read_csv(clinvarcsvfile)
    df = df[[ "ref_aa", "pos_aa", "alt_aa", "ClinicalSignificance", "PhenotypeList",]]
    df = df.sort_values("pos_aa")
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
    
    df["variant"] = df["ref_aa"] + df["pos_aa"].astype(str) +df["alt_aa"]
    # If clinical significance and phenotype information is wanted as well, the hashtag can be removed from the next line
    #df = df[[ "variant" , "ClinicalSignificance", "PhenotypeList"]]
    df = df[[ "variant"]]
    list_of_data = df.values.tolist()
    return list_of_data

# A function to spit out a list of the clinical pathogenic variants and their corresponding theoretical ddg or dddg
def clinvar_csv_to_list_with_stability(clinvarpathogeniccsvpath, clinvarbenigncsvpath, RaSPdatacsvpath, pdbid):
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
    dfgammel = pd.read_csv(oldcsv)
    dfgammel = dfgammel[dfgammel['pdbid'] == oldpdbid]
    dfny = pd.read_csv(newcsv)
    dfny = dfny[dfny['pdbid'] == newpdbid]

    ddggammel = []
    ddgny = []
    dddggammel = []
    dddgny = []

    for row in dfgammel['score_ml']:
        ddggammel.append(row)
    
    for row in dfgammel['score_ml_ddg_bind']:
        dddggammel.append(row)
    
    for row in dfny['score_ml']:
        ddgny.append(row)
    
    for row in dfny['score_ml_ddg_bind']:
        dddgny.append(row)
    
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

    cmap_blended = sns.blend_palette(["#AA61E2", "#40E0D0", "#BADA55", "#FFD700","#F46049" ], as_cmap=True)
    sns.set_theme(style="darkgrid", font="Tahoma")

    # Step 3: Create scatterplot using Seaborn
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=score_ml_gammel, y=score_ml_ny, palette=cmap_blended, hue=matching_rows['score_ml_old'])

    # Step 4: Fit a linear regression model using numpy
    slope, intercept = np.polyfit(score_ml_gammel, score_ml_ny, 1)
    trendline_x = np.linspace(min(score_ml_gammel), max(score_ml_gammel), 100)
    trendline_y = slope * trendline_x + intercept
    plt.plot(trendline_x, trendline_y, color='red', label='Trendline')

    # Step 5: Display R-value (calculate it separately)
    r_value = np.corrcoef(score_ml_gammel, score_ml_ny)[0, 1]
    plt.text(0.7, 0.1, f'R = {r_value:.2f}', transform=plt.gca().transAxes, fontsize=17)

    plt.xlabel(f'Old RaSP data for {gammelpdbid}', fontsize=17)
    plt.ylabel(f'New RaSP data for {nypdbid}', fontsize=17)
    plt.title(f'Old {gammelpdbid} DDG values against new {nypdbid} DDG values', fontsize=18)
    plt.legend()
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

    
    cmap_blended = sns.blend_palette(["#AA61E2", "#40E0D0", "#BADA55", "#FFD700","#F46049" ], as_cmap=True)
    sns.set_theme(style="darkgrid", font="Tahoma")

    # Step 3: Create scatterplot using Seaborn
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=score_ml_ddg_bind_gammel, y=score_ml_ddg_bind_ny, palette=cmap_blended, hue=matching_rows['score_ml_ddg_bind_old'])

    # Step 4: Fit a linear regression model using numpy
    slope, intercept = np.polyfit(score_ml_ddg_bind_gammel, score_ml_ddg_bind_ny, 1)
    trendline_x = np.linspace(min(score_ml_ddg_bind_gammel), max(score_ml_ddg_bind_gammel), 100)
    trendline_y = slope * trendline_x + intercept
    plt.plot(trendline_x, trendline_y, color='red', label='Trendline')

    # Step 5: Display R-value (calculate it separately)
    r_value = np.corrcoef(score_ml_ddg_bind_gammel, score_ml_ddg_bind_ny)[0, 1]
    plt.text(0.7, 0.1, f'R = {r_value:.2f}', transform=plt.gca().transAxes, fontsize=17)

    plt.xlabel(f'Old RaSP data for {gammelpdbid}', fontsize=17)
    plt.ylabel(f'New RaSP data for {nypdbid}', fontsize=17)
    plt.title(f'Old {gammelpdbid} DDDG values against new {nypdbid} DDDG values', fontsize=18)
    plt.legend()
    plt.show()  


def gathering_all_csv_files_into_one(folderpath,newname):
    
    all_files = os.listdir(folderpath)
    
    # Filter out non-CSV files
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
    df = pd.read_csv(csvfile)
    df = df[df['uniprot']== uniprotnr]
    sequence = df['sequence'].iloc[2]
    return sequence 

def align_data():
    ''' i need to align my sequences to the sequences of Nicolas' datasheet.
     to do this i will have to extract the sequence from my RaSP dataframe, i will do this using the sequence_from_csv() funtion
     i will then have to extract the sequence from Nicolas' datasheet, which i can do by just extracting the first sequence of a 
     uniprot-id, using the sequence_from_nicolas function. The sequences of the rightfull protein chains should then be aligned, 
     and i am gonna use Biopythons module pairwise2 to do a global pairwise alignment. From then it is then necessary to fix any 
     numbering of the mutants that could have been changed by this alignment. '''
    
    # my PDB file chain sequences
    opah_a = sequence_from_csv(f'clean_df_ml_all.csv', 'OPAH_A')
    opah_b = sequence_from_csv(f'clean_df_ml_all.csv', 'OPAH_B')
    opah_c = sequence_from_csv(f'clean_df_ml_all.csv', 'OPAH_C')
    opah_d = sequence_from_csv(f'clean_df_ml_all.csv', 'OPAH_D')
    TOo8b_a = sequence_from_csv(f'clean_df_ml_all.csv', '2O8B_A')
    TOo8b_b = sequence_from_csv(f'clean_df_ml_all.csv', '2O8B_B')
    TOdn1_a = sequence_from_csv(f'clean_df_ml_all.csv', '2DN1_A')
    TOdn1_b = sequence_from_csv(f'clean_df_ml_all.csv', '2DN1_B')
    TOo4h_a = sequence_from_csv(f'clean_df_ml_all.csv', '2O4H_A')
    TOo4h_b = sequence_from_csv(f'clean_df_ml_all.csv', '2O4H_B')
    FIREol0_a = sequence_from_csv(f'clean_df_ml_all.csv', '4OL0_A')
    FIREol0_b = sequence_from_csv(f'clean_df_ml_all.csv', '4OL0_B')
    OTTEg7v_a = sequence_from_csv(f'clean_df_ml_all.csv', '8G7V_A')
    OTTEg7v_b = sequence_from_csv(f'clean_df_ml_all.csv', '8G7V_B')
    OTTEg7v_c = sequence_from_csv(f'clean_df_ml_all.csv', '8G7V_C')
    OTTEg7v_d = sequence_from_csv(f'clean_df_ml_all.csv', '8G7V_D')
    ETg82_a = sequence_from_csv(f'clean_df_ml_all.csv', '1G82_A')
    ETg82_b = sequence_from_csv(f'clean_df_ml_all.csv', '1G82_B')
    ETg82_c = sequence_from_csv(f'clean_df_ml_all.csv', '1G82_C')
    ETg82_d = sequence_from_csv(f'clean_df_ml_all.csv', '1G82_D')
    ETi85_a = sequence_from_csv(f'clean_df_ml_all.csv', '1I85_A')
    ETi85_b = sequence_from_csv(f'clean_df_ml_all.csv', '1I85_B')
    ETi85_c = sequence_from_csv(f'clean_df_ml_all.csv', '1I85_C')
    ETi85_d = sequence_from_csv(f'clean_df_ml_all.csv', '1I85_D')

    # my uniprot number chain sequences
    P00439 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'P00439')
    P43246 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'P43246')
    P52701 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'P52701')
    P69905 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'P69905')
    P68871 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'P68871')
    P45381 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'P45381')
    Q9Y5L0 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'Q9Y5L0')
    P62826 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'P62826')
    O95786 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'O95786')
    Q8IUD6 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'Q8IUD6')
    P31371 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'P31371')
    P16410 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'P16410')
    P42081 = sequence_from_nicolas(f'nr2_Nicolas_data.csv', 'P42081')

    # Create a PairwiseAligner object
    aligner = Align.PairwiseAligner()

    # Perform pairwise alignment
    def pairwise_align(seq1, seq2):
        alignment = aligner.align(seq1, seq2)
        best_alignment = max(alignment, key=lambda x: x.score)
        return best_alignment

    # Create a dataframe to store the results
    #df = pd.DataFrame(columns=["pdbid", "uniprot", "alignment_score"])
    columns =  ["pdbid", "uniprot", "alignment_score"]
    results = []

    # Example: Align opah_a with P00439
    alignment1 = pairwise_align(opah_a, P00439)
    results.append({"pdbid": "OPAH_A", "uniprot": "P00439", "alignment_score": alignment1.score})
    alignment2 = pairwise_align(opah_b, P00439)
    results.append({"pdbid": "OPAH_B", "uniprot": "P00439", "alignment_score": alignment2.score})
    alignment3 = pairwise_align(opah_c, P00439)
    results.append({"pdbid": "OPAH_C", "uniprot": "P00439", "alignment_score": alignment3.score})
    alignment4 = pairwise_align(opah_d, P00439)
    results.append({"pdbid": "OPAH_D", "uniprot": "P00439", "alignment_score": alignment4.score})

    alignment5 = pairwise_align(TOo8b_a, P43246)
    results.append({"pdbid": "2O8B_A", "uniprot": "P43246", "alignment_score": alignment5.score})
    alignment6 = pairwise_align(TOo8b_b, P52701)
    results.append({"pdbid": "2O8B_B", "uniprot": "P52701", "alignment_score": alignment6.score})
    alignment7 = pairwise_align(TOo8b_a, P52701)
    results.append({"pdbid": "2O8B_A", "uniprot": "P52701", "alignment_score": alignment7.score})
    alignment8 = pairwise_align(TOo8b_b, P43246 )
    results.append({"pdbid": "2O8B_B", "uniprot": "P43246 ", "alignment_score": alignment8.score})

    alignment9 = pairwise_align(TOdn1_a, P69905)
    results.append({"pdbid": "2DN1_A", "uniprot": "P69905", "alignment_score": alignment9.score})
    alignment10 = pairwise_align(TOdn1_b, P68871)
    results.append({"pdbid": "2DN1_B", "uniprot": "P68871", "alignment_score": alignment10.score})
    alignment11 = pairwise_align(TOdn1_a, P68871)
    results.append({"pdbid": "2DN1_A", "uniprot": "P68871", "alignment_score": alignment11.score})
    alignment12 = pairwise_align(TOdn1_b, P69905)
    results.append({"pdbid": "2DN1_B", "uniprot": "P69905", "alignment_score": alignment12.score})

    alignment13 = pairwise_align(TOo4h_a, P45381)
    results.append({"pdbid": "2O4H_A", "uniprot": "P45381", "alignment_score": alignment13.score})
    alignment14 = pairwise_align(TOo4h_b, P45381)
    results.append({"pdbid": "2O4H_B", "uniprot": "P45381", "alignment_score": alignment14.score})

    alignment15 = pairwise_align(FIREol0_a, Q9Y5L0)
    results.append({"pdbid": "4OL0_A", "uniprot": "Q9Y5L0", "alignment_score": alignment15.score})
    alignment16 = pairwise_align(FIREol0_b, P62826)
    results.append({"pdbid": "4OL0_B", "uniprot": "P62826", "alignment_score": alignment16.score})
    alignment17 = pairwise_align(FIREol0_a, P62826)
    results.append({"pdbid": "4OL0_A", "uniprot": "P62826", "alignment_score": alignment17.score})
    alignment18 = pairwise_align(FIREol0_b, Q9Y5L0)
    results.append({"pdbid": "4OL0_B", "uniprot": "Q9Y5L0", "alignment_score": alignment18.score})

    alignment19 = pairwise_align(OTTEg7v_a, O95786)
    results.append({"pdbid": "8G7V_A", "uniprot": "O95786", "alignment_score": alignment19.score})
    alignment20 = pairwise_align(OTTEg7v_b, Q8IUD6)
    results.append({"pdbid": "8G7V_B", "uniprot": "Q8IUD6", "alignment_score": alignment20.score})
    alignment21 = pairwise_align(OTTEg7v_c, O95786)
    results.append({"pdbid": "8G7V_C", "uniprot": "O95786", "alignment_score": alignment21.score})
    alignment22 = pairwise_align(OTTEg7v_d, Q8IUD6)
    results.append({"pdbid": "8G7V_D", "uniprot": "Q8IUD6", "alignment_score": alignment22.score})
    alignment23 = pairwise_align(OTTEg7v_a, Q8IUD6)
    results.append({"pdbid": "8G7V_A", "uniprot": "Q8IUD6", "alignment_score": alignment23.score})
    alignment24 = pairwise_align(OTTEg7v_b, O95786)
    results.append({"pdbid": "8G7V_B", "uniprot": "O95786", "alignment_score": alignment24.score})
    alignment25 = pairwise_align(OTTEg7v_c, Q8IUD6)
    results.append({"pdbid": "8G7V_C", "uniprot": "Q8IUD6", "alignment_score": alignment25.score})
    alignment26 = pairwise_align(OTTEg7v_d, O95786)
    results.append({"pdbid": "8G7V_D", "uniprot": "O95786", "alignment_score": alignment26.score})

    alignment27 = pairwise_align(ETg82_a, P31371)
    results.append({"pdbid": "1G82_A", "uniprot": "P31371", "alignment_score": alignment27.score})
    alignment28 = pairwise_align(ETg82_b, P31371)
    results.append({"pdbid": "1G82_B", "uniprot": "P31371", "alignment_score": alignment28.score})
    alignment29 = pairwise_align(ETg82_c, P31371)
    results.append({"pdbid": "1G82_C", "uniprot": "P31371", "alignment_score": alignment29.score})
    alignment30 = pairwise_align(ETg82_d, P31371)
    results.append({"pdbid": "1G82_D", "uniprot": "P31371", "alignment_score": alignment30.score})

    alignment31 = pairwise_align(ETi85_a, P16410)
    results.append({"pdbid": "1I85_A", "uniprot": "P16410", "alignment_score": alignment31.score})
    alignment32 = pairwise_align(ETi85_b, P42081)
    results.append({"pdbid": "1I85_B", "uniprot": "P42081", "alignment_score": alignment32.score})
    alignment33 = pairwise_align(ETi85_c, P16410)
    results.append({"pdbid": "1I85_C", "uniprot": "P16410", "alignment_score": alignment33.score})
    alignment34 = pairwise_align(ETi85_d, P42081)
    results.append({"pdbid": "1I85_D", "uniprot": "P42081", "alignment_score": alignment34.score})
    alignment35 = pairwise_align(ETi85_a, P42081)
    results.append({"pdbid": "1I85_A", "uniprot": "P42081", "alignment_score": alignment35.score})
    alignment36 = pairwise_align(ETi85_b, P16410)
    results.append({"pdbid": "1I85_B", "uniprot": "P16410", "alignment_score": alignment36.score})
    alignment37 = pairwise_align(ETi85_c, P42081)
    results.append({"pdbid": "1I85_C", "uniprot": "P42081", "alignment_score": alignment37.score})
    alignment38 = pairwise_align(ETi85_d, P16410)
    results.append({"pdbid": "1I85_D", "uniprot": "P16410", "alignment_score": alignment38.score})

    # Create the dataframe
    alignment_df = pd.DataFrame(results)

    # Save to CSV or any other format you prefer
    alignment_df.to_csv("pairwise_alignments.csv", index=False)
    return print(alignment_df)




 



