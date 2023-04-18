#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 08-02-2023
# Authors:   Sara Vega (saravg99@gmail.com)
#            Alejandra Gonzalez (gonmola@hotmail.es)
# ---


from os.path import join, basename, dirname, exists, splitext
import os
from re import sub

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from qiime2 import Artifact
import qiime2.plugins.feature_classifier.actions as feature_classifier_actions

import utils

def integrate_db(base_db, candidate_db):
    """
    For each entry in candidate db look if the taxonomy is present in base_db.

    - NO: add entry to the database
    - YES: if here is already a taxonomy with this sequence
        - YES: skip
        - NO: check if seq is a substring of other sequences
            - YES: add the new entry to the database if the candidate seq is larger than the base seq.
            - NO: add the new entry to the database


    Arguments
    ----------
        base_db: dict
            Dictionary with taxa and seq keys.
        candidate_db: dict
            Dictionary with taxa and seq keys.

    Returns
    -------
    dict
        Merged database

    See Also
    --------
    utils.load_db_from_dict()

    """

    # INFO
    initial_len = len(base_db)
    print(f"Initial len: {initial_len}")
    new_seqs = 0
    new_taxa = 0
    n_substr = 0

    # Reverse base_db
    base_db_rev = reverse_base_db(base_db)  # base_db_rev = {taxa1: [seq1, seq2]}

    # for each entry in gg
    for ID, item in candidate_db.items():
        tax_str = "-".join(item['taxa'])

        # if taxonomy in base_db
        if tax_str in base_db_rev:
            candidate_seq = item['seq']

            # check if seq is not already
            if not candidate_seq in base_db_rev[tax_str]:

                # check if seq is a substring of other sequences
                for i in range(len(base_db_rev[tax_str])):
                    found_substr = False

                    seq = base_db_rev[tax_str][i]
                    shortest_seq = min(candidate_seq, seq)
                    largest_seq = max(candidate_seq, seq)

                    if shortest_seq in largest_seq:
                        # keep largest seq
                        base_id = [key for key, item in base_db.items() if item['seq'] == seq][0]

                        base_db[base_id]['seq'] = largest_seq
                        base_db_rev[tax_str][i] = largest_seq

                        found_substr = True
                        n_substr = n_substr + 1
                        break

                # if no substring found
                if not found_substr:
                    # add entry to the database
                    base_db[ID] = item
                    base_db_rev[tax_str].append(candidate_seq)

                    # INFO
                    new_seqs = new_seqs + 1

        # If not in base, add entry
        else:
            # add entry to the database
            base_db[ID] = item
            base_db_rev.setdefault(tax_str, []).append(item['seq'])

            # INFO
            new_taxa = new_taxa + 1
            new_seqs = new_seqs + 1

    # INFO
    final_len = len(base_db)
    print(f"Final len: {final_len}")
    print(f"Number of new taxonomies added: {new_taxa}")
    print(f"Number of new sequences added: {new_seqs}")
    print(f"Number of new substr found: {n_substr}")

    return(base_db)

def reverse_base_db(base_db):
    """
    Reverse base database dictionary to facilitate search:
        {taxa: [tax\_str list], seq: SeqObject} -> {tax\_str:[Seq1, Seq2,...]}

    Arguments
    ----------
        base_db: dict
            Dictionary with taxa and seq keys.

    Returns
    -------
    dict
        Reversed database

    See Also
    --------
    integrate_db() Integrates two databases

    """

    base_db_rev = {}

    for ID, item in base_db.items():
        tax_str = "-".join(item['taxa'])

        base_db_rev.setdefault(tax_str, []).append(item['seq'])

    return base_db_rev


def extract_region(seqs, prefix_db, threads=1, force=False, output_formats = ['qza', 'fasta'],
    f_primer = "GTGCCAGCMGCCGCGGTAA", r_primer="GGACTACHVGGGTWTCTAAT", min_length = 200,
    max_length=350):

    """
    Extract region of fasta file, using the primers provided, using QIIME2 api

    Arguments
    ---------
        seqs : str
            Path of the seq file of the database. qza and fasta formats are accepted
        prefix_db: str
            Output prefix name
        output_format: list, default = ['qza', 'fasta']
            Formats to output the extracted sequences. If two formats are provided, the two
            files are generated. Available formats: 'qza', 'fasta'.
        f_primer: str, default = "GTGCCAGCMGCCGCGGTAA"
            Forward primer. Default primer correspond to V4 region of 16S gene.
        r_primer: str, default = "GGACTACHVGGGTWTCTAAT"
            Reverse primer. Default primer correspond to V4 region of 16S gene.
        min_length: int, default = 200
            Minimum length of the extracted region. Extracted sequences shorter than this
            will be discarded.
        max_length: int, default = 350
            Maximum length of the extracted region. Extracted sequences longer than this
            will be discarded.
        threads: int, default = 1
            Number of processes to use.
        force: bool, default = False
            If True, overwrite already existing results


    Returns
    -------
    dict
        Dictionary with output formats and path(s) of the extracted sequence.

    """

    if seqs.endswith('qza'):
        seqs_ar = Artifact.load(seqs)
    elif seqs.endswith('fasta'):
        seqs_ar = Artifact.import_data('FeatureData[Sequence]', seqs)

    seqs_db_out = {outfmt: f"{prefix_db}_seqs.{outfmt}" for outfmt in output_formats}

    if not all(exists(file) for file in seqs_db_out.values()) or force:
        seqs_region_art,  = feature_classifier_actions.extract_reads(
            sequences = seqs_ar,
            f_primer = f_primer,
            r_primer = r_primer,
            min_length = min_length,
            max_length = max_length,
            n_jobs = threads)

        for outfmt in output_formats:
            if outfmt == 'qza':
                # Save to qza
                seqs_region_art.save(f"{prefix_db}_seqs")
            elif outfmt == 'fasta':
                seqs_region_art.save(f"{prefix_db}_seqs")
                utils.unqza(f"{prefix_db}_seqs.qza", f"{prefix_db}_seqs.fasta")
            else:
                raise ValueError(f"Output format \"{outfmt}\" not available")

    return seqs_db_out


def dereplicate(taxa_in, seqs_in, taxa_out, seqs_out, force=False):
    """
    Retain identical sequence records that have differing taxonomies,
    filtering sequence record with the same taxonomy and sequence. It uses
    the the RESCRIPt plugin [1]_.

    Generate qza artifacts, as well as taxonomy txt and sequence fasta files.

    Arguments
    ---------
        taxa_in : str
            Path of qza taxonomy file
        seqs_in: str
            Path of qza sequence file
        taxa_out: str
            prefix of taxonomy file output
        seqs_out: str
            prefix of sequence file output
        force: bool, default = False
            If True, overwrite already existing results

    Returns
    -------
    None

    References
    ----------
    .. [1] Michael S Robeson II, Devon R O'Rourke, Benjamin D Kaehler, Michal Ziemski, Matthew R Dillon, Jeffrey T Foster, Nicholas A Bokulich. 2021. “RESCRIPt: Reproducible sequence taxonomy reference database management.” PLoS Computational Biology 17 (11): e1009581. doi: 10.1371/journal.pcbi.1009581.

    """

    if not exists(taxa_out) or force:
        cmd_dereplicate = f"qiime rescript dereplicate \
        --i-sequences {seqs_in} \
        --i-taxa {taxa_in} \
        --p-rank-handles 'disable' \
        --p-mode 'uniq' \
        --o-dereplicated-sequences {seqs_out}\
        --o-dereplicated-taxa {taxa_out}"

        os.system(cmd_dereplicate)

    # Decompress qza
    seqs_out_fasta = f"{os.path.splitext(seqs_out)[0]}.fasta"
    taxa_out_txt = f"{os.path.splitext(taxa_out)[0]}.txt"

    if not exists(seqs_out_fasta) or not exists(taxa_out_txt) or force:
        utils.unqza(seqs_out, seqs_out_fasta)
        utils.unqza(taxa_out, taxa_out_txt)
