#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 08-02-2023
# Authors:   Sara Vega (saravg99@gmail.com)
#            Alejandra Gonzalez (gonmola@hotmail.es)
# ---



import pandas as pd
from os.path import join, exists, dirname, basename
import matplotlib.pyplot as plt
import seaborn as sns
from utils import rm_prefs, check_dir
from IPython.display import display
import utils
import re
import itertools
import numpy as np




def plot_qualitative(mockname, mock_df, stat = 'fscore', levels = [1,2,3,4,5,6,7],
col_palette='colorblind', group="database", outdir = None, save = True, return_df = False):

    """
    Plot qualitative scores in a lineplot, by level.

    Arguments
    ---------
        mockname: str
            Name of the mock community (for plot title)
        mock_df: dataframe
            Dataframe with qualitative evaluation scores
        stat: str, default = 'fscore'
            Qualitative stat to plot (must be a column of mock_df).
        levels: list, default = [1,2,3,4,5,6,7]
            List with taxonomic levels to plot.
        col_palette: default = 'colorblind'
            Color palette, format compatible with seaborn.
        group: str, default = 'database'
            Grouping feature to plot (must be a column of mock_df).
        outdir: str, optional
            Output directory. Only required if save = True
        save: bool, default = True
            If True, save the plot as png. Requires setting `outdir` to a valid path
        return_df: bool, default = False
            If True, return the dataframe, else, return the plot

    Returns
    -------
    fig
        Line plot of the qualitative scores by level


    Raises
    ------
    ValueError
        If save = True and no output directory is provided.

    See also
    --------
    confusion_matrix_scores() Computes qualitative evaluation scores from a confusion matrix
    """

    if save and not outdir:
        raise ValueError("Please provide a output directory to save the figure.")

    if save:
        # Output
        check_dir(outdir)
        outfile = join(outdir, f'{group}_comparison_{stat}.png')

    stat_df_all = pd.DataFrame()
    l7_order = []

    for lvl in levels:
        try:
            stat_df = mock_df[lvl].groupby(group)[[stat]].mean().sort_values(by=[stat], ascending=False)
            stat_df = stat_df.reset_index()
            stat_df['level']= lvl
        except Exception as e:
            display(mock_df[lvl])
            raise e

        if not stat_df_all.empty:
            stat_df_all = pd.concat([stat_df_all, stat_df], ignore_index=True)
        else:
            stat_df_all = stat_df

        # Order the plot according to stat at last level
        if lvl== levels[-1]:
            l7_order = list(stat_df[group])


    fig = plt.figure(figsize=(11,8), dpi=300)
    sns.lineplot(data=stat_df_all, x='level', y=stat, hue=group, hue_order = l7_order,
                markers=True, legend='brief', style=group, dashes=False,
                palette=col_palette, markersize=10, linewidth=2.5)
    plt.legend(loc='best', bbox_to_anchor=(1, 1))
    plt.title(f"{group.upper()}-{mockname}: Line plot showing the average {stat.upper()} at each level. Legend is descending ordered at level {levels[-1]}")
    plt.ylim(0,1.05)
    if save: plt.savefig(outfile, bbox_inches='tight')
    plt.show()

    if return_df:
        return stat_df_all
    else:
        return fig



def modify_taxabarplots(md, comp_tab, results_dir, mockname, level=7, force=False):


    """Create taxa barplot adjusting taxonomy to match expected taxonomy.

    Clustered taxonomies from observed taxa barplots are replaced by the expected
    taxonomies if both match.

    Arguments
    ---------
        comp_tab: dataframe
            Dataframe containing the expected and observed taxonomies for each feature.
            Expected taxonomy column must be named 'exp_tax'
        md: dataframe
            Metadata conaining classification parameters for each observation.
            Index must be the observed ids.
            It must contain the columns: database, pgram, conf
        results_dir: str
            Output directory
        mockname: str
            Name of the mock community that is being analysed
        level: int, default = 7
            Level at which to collapse taxonomy.
        force = bool, default = False
            If True, overwrite already existing results

    Returns
    -------
    None

    See also
    --------
    merge_taxaplots() Merge taxaplots from different classifications

    """


    from copy import copy

    comp_tab_mod = copy(comp_tab)
    obs_rows = list(md.index)

    # Obtain the abundancy columns for each mock_community sample
    exp_col_index = comp_tab.columns.get_loc("exp_tax")
    mock_samples = comp_tab.iloc[:, 1:exp_col_index]
    n_samples = len(mock_samples.columns)
    mock_sample_list = list(mock_samples.columns)
    comp_tab = comp_tab.fillna("Unknown")


    # Collapse exp taxonomy at selected level
    comp_tab_mod.loc[:, 'exp_tax'] = comp_tab_mod.loc[:, 'exp_tax'].str.split(
        '; ', expand=False).apply(lambda x: '; '.join(x[:level]))



    for obs_name in obs_rows:
        # Collapse observed taxoxnomies
        comp_tab_mod.loc[:, obs_name] = comp_tab_mod.loc[:, obs_name].str.split(
        '; ', expand=False).apply(lambda x: '; '.join(x[:level]))

        db = md.loc[obs_name, 'database']
        params = f"0.001::{md.loc[obs_name, 'pgram']}:{md.loc[obs_name, 'conf']}"

        tbp_out = join(results_dir, db, 'naive-bayes', params,
                       f'taxa-barplot-L{level}_{db}_{params}_{mockname}_mod.tsv')


        if not exists(tbp_out) or force:
            for row_i in comp_tab_mod.index:
                obs = comp_tab_mod.loc[row_i, obs_name]
                exp = comp_tab_mod.loc[row_i, 'exp_tax']

                status = eval_classification_expclust(exp, obs, level)

                if status == 'match':
                    comp_tab_mod.loc[row_i, obs_name] = exp

            # Generate Taxabarplot format
            tbp = comp_tab_mod.loc[:, (obs_name, *mock_sample_list)]

            # Rename columns
            tbp.rename(columns={obs_name:"#OTU ID"}, inplace = True)
            sample = 0
            for mocksample in mock_sample_list:
                tbp.rename(columns={mocksample: f"{mockname}0{sample}S01"},
                          inplace=True)
                sample += 1

            # Collape same taxonomies
            tbp = tbp.groupby(by=["#OTU ID"]).sum().reset_index()
            # Add header
            tbp.columns = pd.MultiIndex.from_tuples(
                    zip(['# Constructed from biom file', '', '', '', '', ''], tbp.columns))
            # Download
            tbp.to_csv(tbp_out, sep="\t", index=False)


def merge_taxaplots(mockname, obs_ele, exp, outdir, level=7, format_taxonomy = False, force=False):

    """
    Merge taxaplots from different classifications of the same features.

    Arguments
    ---------
        mockname: str
            Name of the mock community being analysed
        obs_ele: list of array-like
            List of [path, dataframe] for each taxa barplot to be merged.
            Dataframe taxonomy column must be named 'Taxon'
        exp: dataframe
            Dataframe with the expected taxonomy and its expected abundance
            Taxonomy column must be named 'Taxon'
        outdir: str
            Output directory
        level: int, default = 7
            Taxonomical level to collapse the taxonomies
        format_taxonomy: bool, default = False
            If True, replace taxonomy with its synonym from NCBI taxonomy database
            If False, do not change provided taxonomies
        force: bool, default = False
            If True, overwrite already existing results

    Returns
    -------
        metadata_df: dataframe
            Classification parameters for each observation in the merged taxa barplot
        merged_df: dataframe
            Merged taxa barplot


    See also
    --------
    modify_taxabarplots() Create taxabarplots taking into account expected taxonomy



    """


    from ete3 import NCBITaxa
    import update_taxonomy

    # Counters for colnames
    counter_e = 0
    counter_s = 0

    # Predefined vars
    metadata = []
    rm_cols = ['Kingdom', 'Phylum', 'Class', 'Order','Family', 'Genus', 'Species']

    # Output files
    merged_df_out = join(outdir, f"merged_taxaplot_L{level}_{mockname}.csv")
    metadata_df_out = join(outdir, f"merged_metadata_taxaplot_L{level}_{mockname}.csv")


    if not exists(merged_df_out) or not exists(metadata_df_out) or force:
        # Expected abundancies
        ## Collapse at some Level. Use Kindogm, Phylum,..., col df to collapse
        exp['Taxon'] = exp[rm_cols[:level]].agg('; '.join, axis = 1)
        ## Delete clade cols ( Kindogm, Phylum,...,) and collapse
        exp = exp.drop(columns=rm_cols)
        exp = exp.groupby(by=['Taxon']).sum().reset_index()
        ## Rename columns and save metadata
        for col in exp.columns:
            if col == 'Taxon':
                continue
            new_colname = f"E{counter_e}"
            exp = exp.rename(columns={col: new_colname})

            # sampleID, mocksampleID, database, pgram, confidence
            md_row = (new_colname, col, "Expected", None, None)
            metadata.append(md_row)

            counter_e = counter_e + 1
        ## Create merge dataframe
        merge_df = exp.copy()

        # Observed abundancies
        for item in obs_ele:
            taxbp_path = item[0]
            taxbp_df = item[1]

            # taxonomy assignment params
            path_split = taxbp_path.split('/')
            database = path_split[2]
            clf_params = path_split[4].split(':')
            pgram = clf_params[2]
            confidence = clf_params[3]

            # Format Taxon according to NCBI nomenclature
            if format_taxonomy:
                taxons_up = []

                for taxon in list(taxbp_df['Taxon']):
                    if taxon:
                        new_tax = update_taxonomy.update_taxonomy_ncbi(
                            taxonomy = taxon,
                            join_taxa = True,
                            remove_patterns = True,
                            full_update = True
                            )
                    else:
                        new_tax = taxon

                    new_tax = new_tax.split('; ')[:level]
                    taxons_up.append('; '.join(new_tax))
                # Assign new taxa
                taxbp_df['Taxon'] = taxons_up

            # Collapse same taxonomies
            ## Remove extra cols
            taxbp_df.drop(columns=rm_cols[:level])
            taxbp_df = taxbp_df.groupby(by=['Taxon']).sum().reset_index()

            # Rename columns and save metadata
            for col in taxbp_df.columns:
                if col == 'Taxon':
                    continue
                new_colname = f"S{counter_s}"
                taxbp_df = taxbp_df.rename(columns={col: new_colname})

                # sampleID, mocksampleID, database, pgram, confidence
                md_row = (new_colname, col, database, pgram, confidence)
                metadata.append(md_row)

                counter_s = counter_s + 1

            # merge df
            merge_df = pd.merge(merge_df, taxbp_df, on="Taxon", how='outer')

        # Create metadata dataframe
        colnames =  ["sampleID", "mocksampleID", "database", "pgram", "confidence"]
        metadata_df = pd.DataFrame(metadata, columns=colnames)
        metadata_df.to_csv(metadata_df_out, index=False)

        # Save merged df
        merged_df = merge_df.fillna(0)
        merged_df.to_csv(merged_df_out, index=False)
    else:
        metadata_df = pd.read_csv(metadata_df_out, header=0)
        merged_df = pd.read_csv(merged_df_out, header=0)

    return metadata_df, merged_df



def per_sequence_comparison(obs_tsv_tax, seq_comp_tab, results_dir, mockname,
    format_taxonomy=True, force=False):

    """
    Merge results of several taxonomy assignments of the same features
    with their expected taxonomies and their abundances.
    Also creates a metadata file with the classification details of each observation.

    Arguments
    ------------
        obs_tsv_tax: list
            List of paths of the taxonomy tables obtained in the taxonomy assignments
        seq_comp_tab: dataframe
            Dataframe containing the expected taxonomy and abundance of each feature.
        results_dir: str
            Output directory
        mockname: str
            Name of the mock community being analysed
        format_taxonomy: bool, default = True
            If True, replace taxonomy with its synonym from NCBI taxonomy database
            If False, do not change provided taxonomies
        force: bool, default = False
            If True, overwrite already existing results

    Returns
    -----------
        mockname: str
            Name of the mock community being analysed
        seq_comp_tab: dataframe
            Dataframe with expected taxonomy, expected abundancies and
            observed taxonomies for each feature.
        md: dataframe
            Metadata of the observations.
    """


    import update_taxonomy

    check_dir(results_dir)

    # Metadata list
    md = {}
    counter = 0

    # Output files
    if format_taxonomy:
        md_out = join(results_dir, f"metadata_per_sequence_{mockname}.csv")
        df_out = join(results_dir, f"comparison_taxon_per_sequence_{mockname}.csv")
        # print(md_out)
    else:
        md_out = join(results_dir, f"metadata_per_sequence_NF_{mockname}.csv")
        df_out = join(results_dir, f"comparison_taxon_per_sequence_NF_{mockname}.csv")

    if not exists(md_out) or not exists(df_out) or force:
        for obs_tsv in obs_tsv_tax:
            # Metadata info
            params = obs_tsv.split("/")
            db_name = params[2]
            pgram = params[4].split(":")[2]
            conf = params[4].split(":")[3]

            # print(db_name, pgram, conf)

            # if mockname == 'mockrobiota':
            #     print(obs_tsv)
            #     print(db_name)

            # Load dataframe
            df = pd.read_csv(obs_tsv, sep="\t", header=0)
            # Delete confidence column
            df = df[['Feature ID', 'Taxon']]
            
            # display(df)


            # Unify taxa
            if format_taxonomy:
                taxons = list(df['Taxon'])

                taxons_up = []
                for tax in taxons:

                    # print(f"tax = {tax}")
                    

                    # if taxa comes from GTDB: remove sufixes like _A before formatting
                    # also add space between levels if necessary
                    if 'gtdb' in db_name:
                        pattern = re.compile(r"_[A-Z]([^a-zA-Z])")
                        tax = pattern.sub(r'\1', tax)
                        pattern2 = re.compile(r";([a-z])")
                        tax = pattern2.sub(r'; \1', tax)

                        # print(f"tax_filt = {tax}")


                    new_tax = update_taxonomy.update_taxonomy_ncbi(
                        taxonomy = tax,
                        join_taxa = True,
                        remove_patterns = True,
                        full_update = True
                    )
                    taxons_up.append(new_tax)

                    # print(f"tax_up = {new_tax}")
                    # print(f"tax_up_list = {taxons_up}")
                    # raise

                    # not_found = [""]*7

                    # taxons_up = update_taxa_from_taxa(
                    #     taxa = tax,
                    #     taxa_updated_list = taxons_up, not_found = not_found)
                    #
                    # lowest_level, lowest_rank = update_taxonomy.obtain_lowest_taxa_level_from_taxa(tax)
                    # tax_up = update_taxonomy.get_seven_lineage(lowest_level)
                    #
                    # # If taxa found in NCBI, update, else, keep original and search the next level
                    # if tax_up:
                    #     taxons_up.append(tax_up)
                    # else:
                    #     # tax = tax.replace('d__', 'k__')
                    #     while not tax_up and lowest_rank > 1:
                    #         try:
                    #             not_found[lowest_rank-1] = lowest_level[0].replace(" ", "_")
                    #             # Update lowest_rank
                    #             lowest_rank = lowest_rank - 1
                    #             # Update_lowest_taxa
                    #             lowest_level = [utils.rm_prefs(tax.split('; '))[lowest_rank-1]]
                    #             tax_up = update_taxonomy.get_seven_lineage(
                    #                 lowest_level, add_ranks = not_found)
                    #         except Exception as e:
                    #             print(f"{lowest_rank}")
                    #             print(tax)
                    #             raise e
                    #     if lowest_rank == 0:
                    #         tax = tax.replace('d__', 'k__')
                    #         taxons_up.append(tax)
                    #     else:
                    #         taxons_up.append(tax_up)

                # Update new taxonomy
                df['Taxon'] = taxons_up


            # Change dataframe according to metadata id
            ID = f"obs_0{counter}"
            md[ID] = {"database": db_name, "pgram": pgram, "conf": conf}
            df = df.rename(columns = {'Taxon': ID})

            # Merge observed column
            seq_comp_tab = pd.merge(seq_comp_tab, df, on="Feature ID")


            counter = counter + 1

        md = pd.DataFrame.from_dict(md)
        md = md.transpose()

        # Save files
        md.to_csv(md_out, index=True)
        seq_comp_tab.to_csv(df_out, index=False)
    else:
        md = pd.read_csv(md_out, header=0, index_col=0)
        seq_comp_tab = pd.read_csv(df_out, header=0)

    return mockname, seq_comp_tab, md





def eval_classification_expclust(exp_taxon, obs_taxon, level):
    """
    Compares and expected taxonomy with an observed taxonomy.
    It is able to evaluate clustered taxonomies.

    Arguments
    ---------
        exp_taxon: str
            Expected taxonomy
        obs_taxon: str
            Observed taxonomy
        level:
            Level at which the taxonomies should be compared

    Return
    ------
    str
        Result of the comparison. It can be: match, misclassification, overclassification, underclassification


    """
    # if obs_taxon == np.nan :
    #     return 'misclassification'

    exp_taxon = utils.rm_prefs(exp_taxon.split("; ")[0:level],replace_unknowns=True)
    obs_taxon = utils.rm_prefs(obs_taxon.split("; ")[0:level],replace_unknowns=True)
    # matches counter
    matches = 0
    # compare level by level
    for lvl, (exp, obs) in enumerate(zip(exp_taxon, obs_taxon)):
        match = False
        if lvl != 6:
            exp = exp_taxon[lvl].split('-')
            obs = obs_taxon[lvl].split('-')
            if set(exp) & set(obs):
                match = True
        else:
            # species stuff
            exp = [re.split('_|-', sps) for sps in exp_taxon[lvl].split(':')]
            obs = [re.split('_|-', sps) for sps in obs_taxon[lvl].split(':')]
            if any([len(set(exp_sp)&set(obs_sp)) >= 2 for exp_sp, obs_sp in itertools.product(exp, obs)]):
                match = True

        if match:
            matches += 1
        elif (exp == ['Unknown']):
            return 'overclassification'
        elif (obs == ['Unknown']):
            return 'underclassification'
        else:
            return 'misclassification'


    if matches == len(exp_taxon) and matches == len(obs_taxon):
        return 'match'
    elif (len(exp_taxon)<len(obs_taxon)):
            return 'overclassification'
    elif (len(obs_taxon)<len(exp_taxon)):
            return 'underclassification'
    else:
        raise Exception(f'Not able to evaluate classification.\nExpected: {exp_taxon}\nObserved: {obs_taxon}')




def display_confusion_matrix(confusion_matrix, level = 7, pattern = '', out_img_pref = None,
    fig_size = (22, 17), labelsize = 10, only_mismatches = True, params = None, save = True):

    """ Display a confusion matrix

    Arguments
    ---------
        confusion_matrix: dataframe
            Confusion matrix
        level: int, default = 7
            Level to collapse taxonomies
        pattern: str, optional
            Display only expected taxonomies containing the pattern.
            Is not compatible with `only_mismatches = True`
        out_img_pref: str, optional
            File prefix. Required if save = True.
        fig_size: tuple, default = (22, 17)
            Figure size
        labelsize: int, default = 10
            Label size
        only_mismatches: bool, default = True
            If True, only display mismatches. Not compatible with `pattern`.
        params: str, optional
            String to show on the title
        save: bool, default = True
            If True, save the confusion matrix as an image. Requires `out_img_pref`


    Returns
    --------
    None


    See also
    ---------
    confusion_matrix_scores() Creates a confusion matrix

    """


    if pattern:
        sub = confusion_matrix.filter(regex=pattern).replace(0, np.nan)
        display(sub.dropna())
    else:
        if only_mismatches:
            # Display the mismatches
            # diag_vals = np.diagonal(confusion_matrix)
            # no_matches = [ i for i, val in enumerate(diag_vals) if val == 0]
            # cm_no_matches = confusion_matrix.iloc[no_matches]
            # confusion_matrix = cm_no_matches

    #         display()

            # For mismatches, evaluate if it is a real mistmatches
            # or if contain clusters
            for obs_row in confusion_matrix.index:
                for exp_col in confusion_matrix.columns:
                    value = confusion_matrix.loc[obs_row, exp_col]


                    if value == 0:
                        continue
                    else:
    #                     print(f"{value=}")
    #                     print(f"{exp_col=}")
    #                     print(f"{obs_row=}")
                        comparison = eval_classification_expclust(exp_col, obs_row, level)
    #                     print(f"{comparison=}")
                        # If there is a match, set value to 0
                        if comparison == 'match':
                            confusion_matrix.loc[obs_row, exp_col] = 0
            # Delete empty rows and columns
            confusion_matrix = confusion_matrix.loc[~(confusion_matrix==0).all(axis=1)]
            confusion_matrix = confusion_matrix.loc[:, confusion_matrix.any()]

            if save: out_img = f"{out_img_pref}_MM.png"
        else:
            if save: out_img = f"{out_img_pref}.png"

        if not confusion_matrix.empty:

            # setting the dimensions of the plot
            fig, ax = plt.subplots(figsize=fig_size)
            sns.set(font_scale=0.9)
            sns.heatmap(confusion_matrix, annot=True, annot_kws={"size": 12},
                       cmap="Blues", ax = ax)
            plt.tick_params(axis='both', labelsize=labelsize)
            plt.title(params)
            plt.yticks(rotation=0)
            plt.xlabel("Expected")
            plt.ylabel('Observed')

            if save: plt.savefig(out_img, bbox_inches='tight')

        else:
            print("NO MATRIX (No mismatches)")



def get_best_configurations(mockname, per_sequence_score_obsid_df, group, levels=[1,2,3,4,5,6,7],
    pval_thr=0.05, out_all = True, hyphotesis = 'two-sided'):

    """
    Display best classifications across levels + common best in all levels.

    Arguments
    ---------
        mockname: str
            Name of the mock community being analysed
        per_sequence_score_obsid_df: dict
            Dataframes with the evaluation scores for each observation. Format: ``{level: df}``
        group: str
            Grouping variable to compare. Must be a column in `per_sequence_score_obsid_df`
        levels: list, default = [1,2,3,4,5,6,7]
            Taxonomic levels to evaluate
        pval_thr: float, default = 0.01
            Significance threshold for the statistical tests.
        out_all: bool, default = True
            If True, display results for all steps and levels
            If False, display only the common best configurations across levels
        hypothesis: {'two-sided', 'greater', 'less'}
            Defines the alternative hypothesis. The following options are available (default is ‘two-sided’):
                - 'two-sided': the means of the distributions underlying the samples are unequal.
                - 'less': the mean of the distribution underlying the first sample is less than the mean of the distribution underlying the second sample.
                - 'greater': the mean of the distribution underlying the first sample is greater than the mean of the distribution underlying the second sample.

    Returns
    --------
    dataframe
        Dataframe with the common best configurations across levels


    """


    from scipy.stats import shapiro, levene
    from statistics import mean
    import numpy as np

    import pingouin as pg

    best_classifiers = []
    best_fscores_classifiers = {}
    print_message_once = False

    if out_all: print(f"####\n## {mockname}\n####")
    for lvl, df in per_sequence_score_obsid_df.items():
        if lvl in levels:
            # display(df)
            if out_all: print(f"~~~~\n~~ LEVEL {lvl}\n~~~~")
            df = df.sort_values(by='fscore', ascending=False).reset_index(drop=True)
            df = df.astype({'all_fscores':'string'})
            df = df[~df['all_fscores'].str.contains('nan')].reset_index(drop=True)

            # performs statistics to decide which classifiers are not performig worst than the others
            # if  group == 'database':
            # make dataframe with fscores per obs_id
            groups_dict = { df.loc[i, group] : df.loc[i, 'all_fscores'].split(',') for i in df.index}
            groups_dict = { key: list(map(float, val)) for key, val in groups_dict.items()}
            groups_df = pd.DataFrame(groups_dict)
            groups_df = pd.melt(groups_df)
            # groups_df --> dataframe with 2 colums: `value` (fscore) and `variable` (obs_id)
            if len(groups_df.index) == 5:
                if print_message_once:
                       pass
                else:
                    print_message_once = True
                    print('Watch out! Only one group. No comparisons, only performing average across samples')

                if out_all:
                    print("Your data:")
                    display(groups_df)
                    print(f"AssertionError: Only one group. At least two are needed for comparisons. Therefore, your best configuration in this level is the only one")

                best_confs = list(df[group])
                fscores_best_confs = df[df[group].isin(best_confs)]
                fscores_best_confs = fscores_best_confs[[group, 'fscore', 'precision', 'recall']]
                best_fscores_classifiers[lvl] = fscores_best_confs.set_index(group).T.to_dict('list')
                best_classifiers.append(set(best_confs))
                continue
            # Get pairs to compare
            # 1. get the group name of the configuration with best mean fscore (e.g. get ID name)
            best_conf = df.loc[0, group]



            # 2. get all the comparisons to be made between the best configuration and the others
            comparisons = [(best_conf, obs_conf) for obs_conf in df[group] if best_conf != obs_conf]

            # create list of tuples for storing p-values and create the dataframe
            ## group1, group2, pvalue
            comp_pvals = []

            # Check Normality
            res_norm = pg.normality(data=groups_df, dv='value', group='variable')
            if all(list(res_norm['normal'])):
                if out_all: print("Normality checked")

                # Check Homogenity of variances
                res_homo = pg.homoscedasticity(
                        data=groups_df, dv='value', group='variable', method='levene')

                if all(list(res_homo['equal_var'])):
                    if out_all: print("Homogenity of variances checked")

                    # perform t_test with homoscedasticity one-sided cause we
                    # are only interested in knowing if the observed configuration
                    # (obs_conf) is not worse than the best configuration (best_conf).
                    for pair in comparisons:
                        t_df = pg.ttest(
                            groups_df.loc[groups_df['variable']== pair[0], 'value'],
                            groups_df.loc[groups_df['variable']== pair[1], 'value'],
                            correction = False, alternative = hyphotesis)
                        comp_pvals.append((pair[0], pair[1], list(t_df['p-val'])[0] ))
                else:
                    # perform t_test without homoscedasticity one-sided cause we
                    # are only interested in knowing if the observed configuration
                    # (obs_conf) is not worse than the best configuration (best_conf).
                    for pair in comparisons:
                        t_df = pg.ttest(
                            groups_df.loc[groups_df['variable']== pair[0], 'value'],
                            groups_df.loc[groups_df['variable']== pair[1], 'value'],
                            correction = True, alternative = hyphotesis)
                        comp_pvals.append((pair[0], pair[1], list(t_df['p-val'])[0]))

            else:
                # perform mann withney one-sided cause we are only interested in
                # knowing if the observed configuration (obs_conf) is not worse than
                # the best configuration (best_conf).
                for pair in comparisons:
                    try:
                        man_df = pg.mwu(
                            groups_df.loc[groups_df['variable']== pair[0], 'value'],
                            groups_df.loc[groups_df['variable']==pair[1], 'value'],
                            alternative = hyphotesis)
                        man_pval = list(man_df['p-val'])[0]

                    except ValueError as e:
                        if out_all: print(f"For {pair}: {str(e)}. Pvalue set to 1")
                        man_pval = 1


                    comp_pvals.append((pair[0], pair[1], man_pval ))

            # Store the results for all comparison in dataframe
            comp_pvals_df = pd.DataFrame(comp_pvals,
                columns = ['group1', 'group2', 'pval'])
            comp_pvals_df['reject'] = np.where(
                comp_pvals_df['pval'] < pval_thr, True, False)

            # do pvalue correction
            try:
                corr_res = pg.multicomp(comp_pvals_df['pval'], method = 'fdr_bh',
                alpha=pval_thr)
            except Exception as e:
                display(comp_pvals_df)
                display(groups_df)
                raise e
            comp_pvals_df['padjust'] = corr_res[1]
            comp_pvals_df['reject_adjust'] = corr_res[0]

            if out_all:
                display(comp_pvals_df.loc[comp_pvals_df['reject_adjust'] == False])

            # Only keep the configurations that show no difference with the best one
            best_confs = [best_conf]
            best_confs.extend(list(comp_pvals_df.loc[comp_pvals_df['reject_adjust']==False, 'group2']))
            fscores_best_confs = df[df[group].isin(best_confs)]
            fscores_best_confs = fscores_best_confs[[group, 'fscore', 'precision', 'recall']]

            best_fscores_classifiers[lvl] = fscores_best_confs.set_index(group).T.to_dict('list')
            # For itgdbmock and mockrobiota
            # else:
            #     # Stablish the threshold from the highest fscore
            #     fscore_thr = list(df['fscore'])[0] - window_threshold
            #     df = df[df['fscore'] > fscore_thr]
            #     best_confs = list(df[group])
            #     fscores_best_confs = df[df[group].isin(best_confs)]
            #     # display(fscores_best_confs)
            #     fscores_best_confs = fscores_best_confs[[group, 'fscore', 'precision', 'recall']]
            #
            #     best_fscores_classifiers[lvl] = fscores_best_confs.set_index(group).T.to_dict('list')

            # BEST CLASSIFIERS for each level
            best_classifiers.append(set(best_confs))

            if out_all: print(f"~~~~\n~~ BEST CLASSIFIERS FOR {mockname} level {lvl}\n~~~~")

            if group == "ID":
                if out_all:
                    display(df[df[group].isin(best_confs)][[group, 'database', 'pgram', 'conf', 'fscore', 'precision', 'recall']])
            else:
                if out_all:
                    display(df[df[group].isin(best_confs)][[group, 'fscore', 'precision', 'recall']])


    common_best_classifiers = set.intersection(*best_classifiers)

    print(f"~~~~\n~~ COMMON BEST CLASSIFIERS FOR {mockname}\n~~~~")

    if group == "ID":
        df_common_best = df[df[group].isin(common_best_classifiers)][[group, 'database', 'pgram', 'conf']]
    else:
        df_common_best = df[df[group].isin(common_best_classifiers)][[group]]

    common_best_confs = list(df_common_best[group])

    # Calculate the average fscore across levels and order df_common_best by it
    mean_common_best_confs = []
    mean_precision = []
    mean_recall = []
    for conf in common_best_confs:
        all_common_best = []
        all_precision = []
        all_recall = []
        for lvl in levels:
            # display(best_fscores_classifiers[lvl][conf])
            all_common_best.append(best_fscores_classifiers[lvl][conf][0])
            all_precision.append(best_fscores_classifiers[lvl][conf][1])
            all_recall.append(best_fscores_classifiers[lvl][conf][2])

        mean_common_best_confs.append(mean(all_common_best))
        mean_precision.append(mean(all_precision))
        mean_recall.append(mean(all_recall))

    df_common_best['mean_fscores'] = mean_common_best_confs
    df_common_best['mean_precision'] = mean_precision
    df_common_best['mean_recall'] = mean_recall

    df_common_best = df_common_best.sort_values(by=['mean_fscores'], ascending=False)

    display(df_common_best)

    return df_common_best



def plot_best_configs(mockname, mock_df, best_ids, outdir = None, stat = 'fscore',
                    levels = [1,2,3,4,5,6,7], col_palette = 'colorblind', save = True):

    """
    Plot evaluation scores for the desired observations.


    Arguments
    ---------
        mockname: str
            Name of the mock community being analysed
        mock_df: dict
            Dictionary of dataframes by level, containing the evaluation scores for each
            observation. Format: ``{level: df}``
        best_ids: list
            List of ids of the observations to be plotted
        outdir: str, optional
            Output directory, required if save = True
        levels: list, default = [1,2,3,4,5,6,7]
            List with taxonomic levels to plot.
        col_palette: default = 'colorblind'
            Color palette, format compatible with seaborn.
        save: bool, default = True
            If True, save the plot as png. Requires setting `outdir` to a valid path

    Returns
    -------
    fig
        Line plot of the qualitative evaluation scores by level


    Raises
    ------
    ValueError
        If save = True and no output directory is provided.

    """

    # mock_df --> per_sequence_score_obsid[v_region][mockname]

    if save and not outdir:
        raise ValueError("Please provide an output directory.")

    if save:
        check_dir(outdir)
        outfile = join(outdir, f'best_configs_{stat}.png')

    stat_df_all = pd.DataFrame()
    l7_order = []

    for lvl in levels:
        # filter only best_ids
        stat_df = mock_df[lvl][mock_df[lvl]['ID'].isin(best_ids)]
        stat_df = stat_df.sort_values(by=[stat], ascending=False)
        stat_df['level'] = lvl
        stat_df['config'] = stat_df.database.astype(str) + ':' + stat_df.pgram.astype(str) + ':' + stat_df.conf.astype(str)


        if not stat_df_all.empty:
            stat_df_all = pd.concat([stat_df_all, stat_df], ignore_index=True)
        else:
            stat_df_all = stat_df

        # Order the plot according to stat at last level
        if lvl== levels[-1]:
            l7_order = list(stat_df['config'])


    fig = plt.figure(figsize=(11,8), dpi=80)
    sns.lineplot(data=stat_df_all, x='level', y=stat, hue='config', hue_order = l7_order,
                markers=True, legend='brief', style='config', dashes=False,
                palette=col_palette, markersize=10, linewidth=2.5)
    plt.legend(loc='best', bbox_to_anchor=(1, 1))
    plt.title(f"{mockname}: The best classifiers {stat.upper()} at each level. Legend is descending ordered at level {levels[-1]}")
    plt.ylim(0,1.05)
    if save: plt.savefig(outfile, bbox_inches='tight')
    plt.show()

    return fig


def confusion_matrix_scores(exp, obs, abun, level=7):
    """Compute per sequence evaluation scores with a confusion matrix
    Creates a confusion matrix of expected vs observed taxonomies. For each
    taxonomy, computes the following metrics: true positives (TP), true negatives (TN),
    false positives (FP), false negatives (FN), precision, recall, fscore, accuracy, specificity.

    Arguments
    ---------
        exp: list
            Expected taxonomies of the features
        obs: list
            Observed taxonomies of the features
        abun: list
            Abundancy of each feature
        level: int, default = 7
            Taxonomical level to evaluate

    Returns
    --------
        confusion_matrix: dataframe
            Confusion matrix. Rows are observed taxa and columns are expected taxa
        scores: dataframe
            Computed metrics for each expected taxonomy


    """

    exp = ["; ".join(e.split("; ")[0:level]) for e in exp ]
    obs = ["; ".join(o.split("; ")[0:level]) for o in obs]

    names = list(set(exp).union(set(obs)))


    # create the confusion matrix
    df = pd.DataFrame(0, index=names, columns=names)

    for o, e, a in zip(obs, exp, abun):
        if eval_classification_expclust(e, o, level) == 'match':
            df.loc[e, e] += a
        else:
            df.loc[o, e] +=a

    # remove expected not present in sample (columns all 0):
    df = df.loc[:, (df != 0).any(axis=0)]

    # compute precision, recall and fscore for each taxonomy in expected:

    cols = ['TP', 'TN', 'FP', 'FN', 'precision', 'recall', 'fscore', 'accuracy', 'specificity']
    scores = pd.DataFrame(0, index=df.columns, columns=cols)

    for tax in df.columns:

        scores.loc[tax, 'TP'] = TP = df.loc[tax, tax]
        # false positives: obs = tax, exp != tax
        scores.loc[tax, 'FP'] = FP = df.loc[tax, df.columns !=tax].sum()
        # false negatives: the other way round
        scores.loc[tax, 'FN'] = FN = df.loc[df.index != tax, tax].sum()
        # true negatives
        scores.loc[tax, 'TN'] = TN = df.loc[df.index != tax, df.columns != tax].to_numpy().sum()

        #scores
        scores.loc[tax, 'precision'] = pre = TP/(TP+FP)
        scores.loc[tax, 'recall'] = rec =  TP/(TP+FN)
        scores.loc[tax, 'fscore'] = (2*pre*rec)/(pre+rec)
        scores.loc[tax, 'accuracy'] = (TP+TN)/(TP+TN+FP+FN)
        scores.loc[tax, 'specificity'] = TN/(TN+FP)


    # replace nans with 0
    scores = scores.fillna(0)


    return df, scores


def per_sequence_scores(tab, md, tax_abun, mockname, outdir, level = 7, force = False):

    """
    Computes performance scores using a confusion matrix.
    Scores for each sample and observation are the wheighted average of each taxa.


    Arguments
    ---------
        tab: dataframe
            Table containing feature id, expected taxonomy, feature abundancies and observed taxonomies
            Expected taxonomy column must be named 'exp_tax'.
        md: dataframe
            Metadata with the classification parameters of each observation
        tax_abun: dataframe
            Expected relative abundance for each taxonomy in each sample.
            This will be used as wheight for computing the wheighted average of the scores
        mockname: str
            Name of the mock community being analysed. Sample names must start with this name
        outdir: str
            Output directory
        level: int, default = 7
            Taxonomical level to evaluate
        force: bool, default = False
            If True, overwrite previous results


    Returns
    -------
        mockname: str
            Name of the mock community
        level: int
            Taxonomical level that has been evaluated
        scores: dataframe
            Dataframe with the scores for each observation. Scores are the average of all samples.
            It includes precision, recall, fscore, accuracy, specificity, all_fscores.
            All_fscores contains the fscores for each sample.


    See also
    --------
    per_sequence_comparison() Creates a dataframe with expected taxonomy, abundancies and observed taxonomies for each feature
    confusion_matrix_scores() Creates the confusion matrix and computes the scores per taxonomy
    """

    #reference: https://towardsdatascience.com/confusion-matrix-for-your-multi-class-machine-learning-model-ff9aa3bf7826

    output = f"per_sequence_scores_L{level}_per_obsid_{mockname}_v2.csv"
    outpath = join(outdir, output)

    if not exists(outpath) or force:

        samples = [col for col in tab.columns if col.startswith(mockname)]
        new_tab = pd.DataFrame.to_dict(md, orient = 'index')

        for obs_id in md.index:
            precision = []
            recall = []
            fscore = []
            accuracy = []
            specificity = []

            for sample in samples:
                exp = list(tab['exp_tax'])
                obs = list(tab[obs_id])
                abun = list(tab[sample])

                confusion_matrix, scores = confusion_matrix_scores(exp, obs, abun, level)

                scores = scores.join(tax_abun[sample])

                precision.append((scores['precision']*scores[sample]).sum())
                recall.append((scores['recall']*scores[sample]).sum())
                fscore.append((scores['fscore']*scores[sample]).sum())
                accuracy.append((scores['accuracy']*scores[sample]).sum())
                specificity.append((scores['specificity']*scores[sample]).sum())


            new_tab[obs_id]['precision'] = np.mean(precision)
            new_tab[obs_id]['recall'] = np.mean(recall)
            new_tab[obs_id]['fscore'] = np.mean(fscore)
            new_tab[obs_id]['accuracy'] = np.mean(accuracy)
            new_tab[obs_id]['specificity'] = np.mean(specificity)
            new_tab[obs_id]['all_fscores'] = ', '.join(map(str, fscore))


        new_tab_df = pd.DataFrame.from_dict(new_tab, orient='index').reset_index()
        new_tab_df = new_tab_df.rename(columns = {'index':'ID'})
        new_tab_df.to_csv(outpath, index=False)

    else:
        new_tab_df = pd.read_csv(outpath, header=0)

    return mockname, level, new_tab_df
