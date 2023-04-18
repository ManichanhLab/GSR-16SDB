#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 03-02-2023
# Authors:   Sara Vega (saravg99@gmail.com)
#            Alejandra Gonzalez (gonmola@hotmail.es)
# ---


###########
# IMPORTS #
###########

from ete3 import NCBITaxa
import re
import sys
import utils

# NCBI db
ncbi = NCBITaxa()

# Define 7-level taxa + formating prefixes
seven_lineage = {'superkingdom':'k__', 'phylum':'p__', 'class':'c__', 'order':'o__',
    'family':'f__', 'genus':'g__', 'species':'s__'}

seven_lineage_list = ['superkingdom', 'phylum', 'class', 'order',
    'family', 'genus', 'species']


#############
# FUNCTIONS #
#############

def get_taxa_from_specie_taxid(sp_taxid):
    """
    Takes a specie taxID and returns a list with all the taxonomy levels
    
    Arguments
    ----------
        sp_taxid: list
            List containing the specie taxID
    
    Returns
    -------
    list
        Seven-lineage taxonomy (i.e., Kingdom, Phylum,..., Specie) 
    
    """

    lineage = ncbi.get_lineage(sp_taxid) # contains the original order of taxids
    names = ncbi.get_taxid_translator(lineage)
    # Get only 7-lineage (taxids)
    ranks = ncbi.get_rank(lineage)
    # invert the ranks dictionary
    inv_ranks = {v: k for k, v in ranks.items()}

    value = []

    for rank in seven_lineage_list:
        if rank in inv_ranks:
            id = inv_ranks[rank]
            value.append(names[id])
        else:
            value.append("")

    # replace spaces by _
    value = list(map(lambda x: x.replace(" ", "_"), value))

    return value



def obtain_lowest_taxa_level(taxonomy, return_rank = False, remove_patterns = True):
    """
    Get the lowest taxa available level of a taxonomy
    
    Arguments
    ----------
        taxonomy: list or str
            Taxonomy to parse
        return_rank: bool, default = False
            If True, the rank of the lowest taxa found is returned.
        remove_patterns: bool, default = True
            Whether to remove or not unknown patterns (aka Unknown, uncultured, etc)
        
    Returns
    -------
    str
        Lowest taxa level
    
    """
        
    if type(taxonomy) is str:
        # get taxonomy ranks list
        taxonomy = taxonomy.split("; ")

    # clean the names (aka remove k__ etc)
    clean_taxa = utils.rm_prefs(taxonomy, replace_unknowns=True)

    # remove empty levels
    clean_taxa = list(filter(None, clean_taxa))


    if remove_patterns:
        # remove levels == Unknown (aka empty levels)
        clean_taxa = list(filter('Unknown'.__ne__, clean_taxa))

        # discard levels with rare patters (e.g. uncultured, metagenome, unidentified,...)
        rare_pattern_file = "/mnt/synology/DATABASES/QIIME2/SSU/gsp_99/db_evaluation/discard_patterns.txt"

        ## Open the file and extract the patterns
        with open(rare_pattern_file, 'r') as fh:
            file = fh.read()
            patterns = file.split()
        ## Filter taxa with the patterns
        clean_taxa = [tax for tax in clean_taxa if not any(pat in tax for pat in patterns)]


    # Obtain the lowest level (last element of the list)
    if clean_taxa:
        lowest_taxa_level = clean_taxa[-1].replace("_", " ")
        if return_rank:
            lowest_rank = utils.rm_prefs(taxonomy).index(clean_taxa[-1]) +1
            return lowest_taxa_level, lowest_rank
        else:
            return lowest_taxa_level


def update_taxonomy_ncbi(taxonomy, join_taxa = True, remove_patterns = True, full_update = True):
    """
    Get the lowest taxa available level of a taxonomy
    
    Arguments
    ----------
        taxonomy: list or str
        join_taxa: bool, default = True
            If True, returns a string with prefixes, if False, returns a list without prefixes.
        remove_patterns: bool, default = True
            Whether to remove or not unknown patterns (aka Unknown, uncultured, etc)
        full_update: bool, default = True
            wWether to take into account the lowest taxa level or check all the levels for synonyms
        
    Returns
    -------
    str
        Updated taxonomy
    
    """

    lowest_taxa_level, lowest_rank = obtain_lowest_taxa_level(taxonomy, return_rank = True, remove_patterns = remove_patterns)
#     print(lowest_taxa_level)
#     print(lowest_rank)
    taxa_up = get_seven_lineage(lowest_taxa_level, join_taxa = join_taxa)

    if taxa_up or not full_update:
        return taxa_up

    else:
        not_found = [""]*7
        if type(taxonomy) is str:
            taxonomy = taxonomy.split('; ')
        taxonomy = utils.rm_prefs(taxonomy)
        # taxonomy now is a list without prefixes

        while not taxa_up and lowest_rank > 1:
            try:
                not_found[lowest_rank-1] = lowest_taxa_level.replace(" ", "_")
                # Update lowest_rank
                lowest_rank = lowest_rank - 1
                # Update_lowest_taxa
                lowest_taxa_level = taxonomy[lowest_rank-1].replace("_", " ")
                taxa_up = get_seven_lineage(
                    lowest_taxa_level, add_ranks = not_found, join_taxa = join_taxa)
                # print(lowest_taxa_level)
                # print(f"{lowest_rank}")
                # print(not_found)
                # print(taxa_up)

            except Exception as e:
                print(f"{lowest_rank}")
                print(not_found)
                print(taxa)
                raise e
        if lowest_rank == 1 and not taxa_up:
            # no rank was found, return nothing
            return ""
        else:
            return taxa_up




def get_seven_lineage(lowest_taxa_level, add_ranks = None, join_taxa = True):  # join_taxa = True
    """
    Get the taxanomy from a single taxa name.
    
    Arguments
    ----------
        lowest_taxa_level: str or list
            Taxa name (e.g, specie, genus,...)
        add_ranks: list (optional)
            list with 7 elements, contaning taxa ranks to add. 
            Ranks to conserve must be empty strings.
            e.g. if we want to add the species level: add_ranks = ['', '', '', '', '', '', 'species']
        join_taxa: bool, default = True
            If True, returns a string with prefixes, if False, returns a list without prefixes.
            
    
    Returns
    -------
    str or list
        Seven-lineage taxa

    """
    # make sure that input lowest taxa level is a List
    # if it is a str, convert to list
    if type(lowest_taxa_level) is str:
        lowest_taxa_level = [lowest_taxa_level]

    # Get taxid from lowest taxa level
    name2taxid = ncbi.get_name_translator(lowest_taxa_level)
    # If name found
    if name2taxid:
        # Get full lineage (taxids)
        id = list(name2taxid.values())[0][0]
        lineage = ncbi.get_lineage(id) # contains the original order of taxids
        names = ncbi.get_taxid_translator(lineage)
        # Get only 7-lineage (taxids)
        ranks = ncbi.get_rank(lineage)
        # invert the ranks dictionary
        inv_ranks = {v: k for k, v in ranks.items()}

        value = []

        # join the names with the prefixes
        for rank in seven_lineage_list:
            if rank in inv_ranks:
                id = inv_ranks[rank]
                #value.append(f"{seven_lineage[rank]}{names[id]}")
                value.append(names[id])
            else:
                # value.append(seven_lineage[rank])
                value.append("")


        # replace spaces by _
        value = list(map(lambda x: x.replace(" ", "_"), value))

        if add_ranks:
            value = ["".join(item) for item in zip(value, add_ranks)]

        if join_taxa:
            # join all names in a string
            taxa  = utils.join_taxa_lineages(value)
        else:
            taxa = value

    # If name not found:
    else:
        taxa = ""

    return taxa
