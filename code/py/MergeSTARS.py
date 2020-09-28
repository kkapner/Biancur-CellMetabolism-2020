import pandas as pd
import os
import numpy as np
import sys

def merge_stars(data_folder):
    """
    This function merges the data from the negative and positive direction STARS data by taking the data corresponding to the values with the highest p-value.

    :param data_folder: full path to folder containing combined STARS data (postivie and negative directions appended together)
    :return: saves the merged data over the original data so that all genes appear once and have one corresponding STARS score.
    """
    files = os.listdir(data_folder)

    for f in files:
        if not f.startswith(".") and not os.path.isdir(os.path.join(data_folder, f)):
            path = os.path.join(data_folder, f)
            all_data = pd.read_csv(path, sep = "\t")
            merged_data = pd.DataFrame()
            print(f"On file {f}")
            for gene in all_data['Gene.Symbol'].unique():
                gene_data = all_data.loc[all_data['Gene.Symbol'] == gene]

                max_index = gene_data["p.value"].idxmax()
                min_index = gene_data["p.value"].idxmin()
                
                # Taking greater p-value unless one of them is signifcant with
                # an FDR of less than 0.25
                if all_data.iloc[min_index]["p.value"] < 0.05 and all_data.iloc[min_index]["FDR"] < 0.25:
                    merged_data = merged_data.append(all_data.iloc[min_index])

                else:
                    merged_data = merged_data.append(all_data.iloc[max_index])


            merged_data = merged_data[["Gene.Symbol", "STARS.Score", "p.value", "FDR", "q.value"]].sort_values("STARS.Score", ascending = False).reset_index(drop = True)
            merged_data.to_csv(os.path.join(data_folder, f), sep = "\t", index = False)

    return

def read_csv(csv_file):
    """
    This function reads in any number of csv files (list or one) and returns a dictionary of all of the data in the individual csv files. The key is the name of the file and the value is a pandas data frame.

    :param csv_file: single or list of full file paths to csv files
    :return: dictionary of pandas data frames with keys as file names and values as the corresponding data from the csv file
    """
    if type(csv_file) == str:
        return pd.read_csv(csv_file, sep = "\t")

    elif len(csv_file) > 1:
        data = {}
        for file in list(csv_file):
            fname = os.path.basename(file)
            data[fname] = pd.read_csv(file, sep = "\t")

        return data
    
    elif csv_file:
        return pd.read_csv(csv_file[0], sep = "\t")

    else:
        print("Issue with file locations.")
        sys.exit()

def data_extract(dataframe, keys):
    cols = dataframe.columns.tolist()
    gene_col = [x for x in cols if "gene" in x.lower()]
    gene_col = gene_col[0]
    return dataframe.loc[dataframe[gene_col].isin(keys)].sort_values(by = [gene_col]).set_index(gene_col)

def add_lfc(data_path, reference_path):
    """
    This function takes in the STARS files with a STARS score, p.value, fdr, and q.value for each gene and then appends the median LFC from the reference file provided by reference_path.

    :param data_path: full file path to the location containing the STARS files.
    :param reference_path: full file path to the reference log fold change data.

    :return: a copy of each of the STARS files with the appended median LFC from the reference_file.
    """
    stars_files = os.listdir(data_path)
    stars_files = [os.path.join(data_path, file) for file in stars_files if (
        not file.startswith(".") and not os.path.isdir(os.path.join(data_path, file)))]
    data = read_csv(stars_files)
    normalized_lfc = read_csv(reference_path)

    for key in data:
        merged_data = data[key].sort_values(by="Gene.Symbol")
        m_d_genes = merged_data["Gene.Symbol"].tolist()
        filtered_data = data_extract(normalized_lfc, m_d_genes)
        sample_name = key[key.rfind("-")+1:key.rfind(".")]
        median_lfc = filtered_data.groupby(
            "GeneID").agg({sample_name: np.median})
        merged_data[sample_name] = median_lfc.values
        merged_data.to_csv(os.path.join(data_path, f"STARS-with-Median-LFC-{sample_name}.csv"), 
                           index=False)

    return

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Merge STARS data')
    parser.add_argument('DataPath', metavar='Y', type=str, 
                        help='The full file path to the merged STARS data')
    parser.add_argument('ReferenceData', metavar='Y', type=str,
                        help='The full file path to the reference log fold change data to extract the median LFC from')
    args = parser.parse_args()

    merge_stars(args.DataPath)
    add_lfc(args.DataPath, args.ReferenceData)
