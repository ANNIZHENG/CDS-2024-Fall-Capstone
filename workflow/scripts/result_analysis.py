import pandas as pd

def find_opposite_sign(file1, file2, file3):
    """
    Reads three CSV files, merges data, and identifies cases where beta estimates
    for two enhancers are of the same sign, but the combined interaction is of the opposite sign.

    Args:
    - file1 (str): Path to the first CSV file.
    - file2 (str): Path to the second CSV file.
    - file3 (str): Path to the third CSV file.

    Returns:
    - pd.DataFrame: A DataFrame containing cases of opposite sign interaction.
    """
    enhancer1 = pd.read_csv(file1, header=0, names=["gene_1", "enhancer1_1", "enhancer2_1", "intercept_1", 
                                                    "beta.estimate_1", "beta.pvalue_1", "bootstrap.pvalue_1"])
    enhancer2 = pd.read_csv(file2, header=0, names=["gene_2", "enhancer1_2", "enhancer2_2", "intercept_2", 
                                                    "beta.estimate_2", "beta.pvalue_2", "bootstrap.pvalue_2"])
    interaction = pd.read_csv(file3, header=0, names=["gene", "enhancers", "intercept", 
                                                      "beta.estimate", "beta.pvalue", "bootstrap.pvalue"])

    interaction[['enhancer1', 'enhancer2']] = interaction['enhancers'].str.split('_', expand=True)
    interaction = interaction.drop(columns=['enhancers'])

    merged_df = pd.merge(enhancer1, enhancer2, left_on=['gene_1', 'enhancer1_1', 'enhancer2_1'], 
                         right_on=['gene_2', 'enhancer1_2', 'enhancer2_2'])
    merged_df = pd.merge(merged_df, interaction, left_on=['gene_1', 'enhancer1_1', 'enhancer2_1'], 
                         right_on=['gene', 'enhancer1', 'enhancer2'])

    merged_df = merged_df[['gene_1', 'enhancer1_1', 'enhancer2_1', 'intercept_1',
                           'beta.estimate_1', 'beta.pvalue_1', 'bootstrap.pvalue_1', 'intercept_2', 
                           'beta.estimate_2', 'beta.pvalue_2', 'bootstrap.pvalue_2', 'intercept',
                           'beta.estimate', 'beta.pvalue', 'bootstrap.pvalue']]
    merged_df.columns = ['gene', 'enhancer1', 'enhancer2', 'intercept_enhancer1', 'beta.estimate_enhancer1',
                         'beta.pvalue_enhancer1', 'bootstrap.pvalue_enhancer1', 'intercept_enhancer2', 
                         'beta.estimate_enhancer2', 'beta.pvalue_enhancer2', 'bootstrap.pvalue_enhancer2', 
                         'intercept_both', 'beta.estimate_both', 'beta.pvalue_both', 'bootstrap.pvalue_both']

    significant = merged_df[(merged_df['beta.pvalue_enhancer1'] >= 0.05) & 
                            (merged_df['beta.pvalue_enhancer2'] >= 0.05) & 
                            (merged_df['beta.pvalue_both'] >= 0.05)]

    betas = significant[['gene', 'enhancer1', 'enhancer2', 'beta.estimate_enhancer1',
                          'beta.estimate_enhancer2', 'beta.estimate_both']]
    betas['opposite_sign_interaction'] = (
        ((betas['beta.estimate_enhancer1'] > 0) & (betas['beta.estimate_enhancer2'] > 0) & 
         (betas['beta.estimate_both'] < 0)) |
        ((betas['beta.estimate_enhancer1'] < 0) & (betas['beta.estimate_enhancer2'] < 0) & 
         (betas['beta.estimate_both'] > 0))
    )

    opposite_sign = betas[betas['opposite_sign_interaction']]
    
    return opposite_sign

# example
if __name__ == "__main__":
    file1 = '/Users/huajingru/Desktop/Fall_2024/Capstone/CDS-2024-Fall-Capstone/workflow/scripts/model_results/model_results_cells_01.csv'
    file2 = '/Users/huajingru/Desktop/Fall_2024/Capstone/CDS-2024-Fall-Capstone/workflow/scripts/model_results/model_results_cells_10.csv'
    file3 = '/Users/huajingru/Desktop/Fall_2024/capstone store/hpc_model_results/CD8-naive/model_results_cells_00.csv'
    
    result = find_opposite_sign(file1, file2, file3)
    print(result)