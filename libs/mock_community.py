#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 08-02-2023
# Authors:   Sara Vega (saravg99@gmail.com)
#            Alejandra Gonzalez (gonmola@hotmail.es)
# ---


from qiime2 import Artifact
import qiime2.plugins.demux.actions as demux_actions
import qiime2.plugins.dada2.actions as dada2_actions
import qiime2.plugins.taxa.actions as taxa_actions


from Bio import SeqIO

from copy import copy, deepcopy
import pandas as pd
import os
from os.path import join, basename, dirname, exists, splitext
import random
import subprocess
import gzip
import shutil

# For parallelization
from joblib import Parallel, delayed

# Utils
from utils import check_dir, unqza
import utils



def eval_ref(abundance_tab, ref_tax, ref_seqs):
    
    """
    Evaluate the reference file according to abundance table.

    - Keep the taxa only present in abundance table
    - Check if the abundance of the taxa in the reference file match with 
    the abudance table

    Arguments
    ----------
        abundance_tab: dict
            Desired taxonomy abundances. Dictionary with format {taxa: abundance}
        ref_tax: dataframe
            Reference database taxonomy.
        ref_seqs: dict
            Reference database sequences. Dictionary with format {id: Bio.Seq object}

    Returns
    -------
    dict 
        Dictionary containing the taxa entry and its database seqs {taxa: [Bioseq1, Bioseq2,...]}
        with abundancies matching the input abundance.
    
    """


    random.seed(123)
    # STEP1: get this {taxa:[seq1, seq2, ...]}

    # Filtrate ref taxa to contain only tax present in abundance_tab
    sub_tax = ref_tax[ref_tax['Taxon'].isin(set(abundance_tab.keys()))]
    sub_tax = sub_tax.set_index(sub_tax.columns[0])
    # Filtrate ref seqs to contain only tax present in abundance tab
    ## Agrupate feature ids according to taxa ({taxa: [id1, id2,...]})
    tax_ids = sub_tax.groupby(by=["Taxon"]).groups
    tax_ids = {key:list(value) for key, value in tax_ids.items()}
    # Agrupate seqs according to taxa ({taxa: [seq1, seq2,...]})
    tax_seqs= {taxa:[ ref_seqs[val] for val in ids ] for taxa, ids in tax_ids.items()}
    # Add description (taxa) to SeqRecord object
    for taxa in tax_seqs.keys():
        for sq in tax_seqs[taxa]:
            sq.description = taxa

    # STEP2: check abundances and adjust them
    env_seqs = {}

    for tax, abun in abundance_tab.items():
        if tax in set(tax_seqs.keys()):
            # Multiply abudance to avoid low numbers
            abun = abun * 100
            # If there are less sequences than required abundance for a taxon
            if abun > len(tax_seqs[tax]):
                # take N random sequences from seqs with replacement
                mylist = random.choices(tax_seqs[tax], k=round(abun))
                env_seqs[tax] = [deepcopy(item) for item in mylist]
            else:
                # take N random sequences from seqs without replacement
                env_seqs[tax] = random.sample(tax_seqs[tax], round(abun))


    return env_seqs




def simulate_env(abun_table, db, db_name, out_dir, force=False):
    """
    Generate taxonomy (tsv) and sequence (FASTA) files from abundancy_table.
    The abundance of FASTA sequences is directly proportional to the abundances
    especified in the input table.
    It generates as many samples as specified in the abundance table. 

    Arguments
    ---------
        abun_table: dataframe
            Dataframe containing taxonomical abudancies for each sample. First column must be named 'Taxon',
            next columns must be named after the sample names.
        db: dict
            Reference database files in dictionary format as follows: ``{db_name:{taxa:taxafile, seq:fastafile}}``
        db_name: str
            Reference database name. It must match the database key in the `db` object
        out_dir: str
            Directory where to store the output files
        force: bool
            Whether to overwritte or not already existing files

    Returns
    -------
        env_out_seqs_dict: dict
            Dictionary with the fasta file paths for each simulated sample.
        new_abun: dataframe
            Dataframe with the exact resulting environment abundances. 

    See Also
    --------
    eval_ref() Evaluate abundance file and correct sequence abundancies. 

    """

    # 1. Adapt input and output files
    env_out_seqs_dict = {}
    env_new_abun = {}
    ## Input files
    for col in abun_table.columns[1:]:
        abundance_tab = dict(zip(abun_table['Taxon'], abun_table[col]))
        # Remove taxa with abundance 0
        abundance_tab = {key:value for key, value in abundance_tab.items() if value != 0}

        env_new_abun[col] = {}

        # Read files
        ref_tax  = pd.read_csv(db[db_name]['taxa'], sep = '\t') # df
        ref_seqs = {rec.id : rec for rec in SeqIO.parse(db[db_name]['seq'], "fasta")} # dictionary

        ## Output files
        ### Sequences of the created environment
        env_out_seqs  = join(out_dir,f"env-seqs_{col}.fasta")
        ### Sequence identifiers of the environment for each taxa
        env_out_ids = join(out_dir,f"env-identifiers_{col}.tsv")
        ### Recalculated abundances after creating the environment
        env_abundance = join(out_dir, f"env_abundance.tsv")

        if not os.path.exists(env_out_seqs) or not os.path.exists(env_out_ids) or force:
            #  Evaluate reference files and check abundances
            env_seqs = eval_ref(abundance_tab, ref_tax, ref_seqs)

            # Assign identifiers to each sequence and downloand environment files
            i=0
            ## tax_seqs in fasta format
            with open(env_out_seqs, "w") as out_seqs:
                ## tax_ids in tsv format
                with open(env_out_ids, "w") as out_taxa:
                    # Add header
                    out_taxa.write(f"Taxa\tAbundance\tIDs\n")
                    for tax, seqs in env_seqs.items():
                        id_list = []
                        for seq in env_seqs[tax]:
                            ID = f"id0{i}"
                            seq.id = ID
                            id_list.append(ID)

                            i = i+1

                        # Taxa file
                        if id_list:
                            num_seqs = len(id_list)
                            env_new_abun[col][tax] = num_seqs
                            id_line = ",".join(id_list)
                            out_taxa.write(f"{tax}\t{num_seqs}\t{id_line}\n")
                        # Fasta File
                        SeqIO.write(env_seqs[tax], out_seqs, "fasta")
        else:
            env_out_ids_df = pd.read_csv(env_out_ids, header = 0, sep = '\t')
            env_new_abun[col] = {tax: num_seqs for tax, num_seqs in zip(env_out_ids_df.Taxa, env_out_ids_df.Abundance)}

        env_out_seqs_dict[col] = env_out_seqs

    if not exists(env_abundance) or force:
        new_abun = pd.DataFrame.from_dict(env_new_abun).fillna(0)
        new_abun = new_abun.reset_index()
        new_abun = new_abun.rename(columns = {'index': 'Taxon'})
        new_abun.to_csv(env_abundance, sep="\t", index=False)
    else:
        new_abun = pd.read_csv(env_abundance, sep="\t")

    return env_out_seqs_dict, new_abun



def simulate_repseqs(abun_table, db, db_name, out_dir, force=False):

    """ Simulate environment with the desired taxonomical abundances per sample

    Given a taxonomical abundancy table and a reference database, generates the following environment files:
            
        - env_seqs.fasta and env_taxa.txt: environment sequences and its taxonomy
        - env_tab.tsv: feature table with the sequence abundancies per sample
        - env_abundance: final environment abundancy table

    The files ``env_seqs.fasta`` and ``env_tab.tsv`` are ready to be imported into qiime2 artifacts, 
    and to be used as a substitute of the DADA2 output. 


    Arguments
    ---------
        abun_table: dataframe
            Dataframe containing the desired taxonomical abudancies for each sample. 
            First column must be named 'Taxon', next columns must be named after the sample names
        db: dict
            Reference database files in dictionary format as follows: ``{db_name:{taxa:taxafile, seq:fastafile}}``
        db_name: str
            Reference database name. It must match the database key in the `db` object
        out_dir: str
            Directory where to store the output files
        force: bool
            Whether to overwritte or not already existing files

    Returns
    -------
        abun_tab: dataframe
            Environment abundancy table: environment abundance of each taxa in each sample

        feature_tab: dataframe
            Environment feature table: environtment abundance of each feature in each sample

    """


    # initialize variables
    env_ids = set() #ids to be included in the environment "repseqs"
    env_tab = {} # dictionary to store the abundance of each id in each sample
    env_abun = {} # dictionary to store the abundance of each taxonomy in each sample


    # output files

    env_tab_file = join(out_dir, 'env_tab.tsv')
    env_abun_file = join(out_dir, 'env_abundance.tsv')
    env_taxa_file = join(out_dir, 'env_taxa.txt')
    env_seqs_file = join(out_dir, 'env_seqs.fasta')



    if any(not exists(file) for file in [env_tab_file, env_abun_file, env_taxa_file, env_seqs_file]) or force:
        # Read reference files
        ref_tax  = pd.read_csv(db[db_name]['taxa'], sep = '\t') # df
        ref_seqs = {rec.id : rec for rec in SeqIO.parse(db[db_name]['seq'], "fasta")} # dictionary

        for col in abun_table.columns[1:]:
            # new entry for this sample
            env_tab[col] = {}

            # subset abundance table and turn into a dictionary
            abundance_tab = dict(zip(abun_table['Taxon'], abun_table[col]))
            # Remove taxa with abundance 0
            abundance_tab = {key:value for key, value in abundance_tab.items() if value != 0}
            # get sequences according to abundance
            env_seqs = eval_ref(abundance_tab, ref_tax, ref_seqs)

            # get list of all sequences included in environment
            all_ids = [seq.id for tax, seqs in env_seqs.items() for seq in seqs]
            # add ids to env_ids set
            env_ids.update(all_ids)
            # add ids + its abundance in the env_tab table
            env_tab[col] = {ID: all_ids.count(ID) for ID in all_ids}

            # add taxonomy + number of sequences (abundance) to the env_abun table
            env_abun[col] = {tax: len(seqs) for tax, seqs in env_seqs.items()}


        # create output files #

        # load database, subset and download

        ref_db = utils.load_db_from_files(taxa_path = db[db_name]['taxa'], seqs_path = db[db_name]['seq'])
        env = {key: items for key, items in ref_db.items() if key in env_ids}
        utils.download_db_from_dict(env,
                                    taxa_out_file= env_taxa_file,
                                    seqs_out_file = env_seqs_file)
        print(f"Saved {env_taxa_file}")
        print(f"Saved {env_seqs_file}")


        # download feature table

        env_feature_tab = pd.DataFrame.from_dict(env_tab).fillna(0)
        env_feature_tab = pd.DataFrame.from_dict(env_feature_tab).fillna(0)
        env_feature_tab = env_feature_tab.reset_index()
        env_feature_tab = env_feature_tab.rename(columns = {'index': 'ID'})
        env_feature_tab.to_csv(env_tab_file, sep="\t", index=False)
        print(f"Saved {env_tab_file}")
#         display(env_feature_tab)

        # download abundance table
        env_abun_tab = pd.DataFrame.from_dict(env_abun).fillna(0)
        env_abun_tab = pd.DataFrame.from_dict(env_abun_tab).fillna(0)
        env_abun_tab = env_abun_tab.reset_index()
        env_abun_tab = env_abun_tab.rename(columns = {'index': 'Taxon'})
        env_abun_tab.to_csv(env_abun_file, sep="\t", index=False)
        print(f"Saved {env_abun_file}")
#         display(env_abun_tab)

    else:
        env_abun_tab = pd.read_csv(env_abun_file, sep="\t")
        env_feature_tab = pd.read_csv(env_tab_file, sep="\t")

    return env_abun_tab, env_feature_tab

def art_simulation(input_fasta, output_prefix, dir_out, art_path, seq_method = "MSv1", length = 200, count = 1,
                   n_samples = 1, mlen = 400, paired = False, seed=123, minQ=25,
                   force=False):
    """
    Simulate Miseq Illumina Reads with ART simulation tools [1]_.
    The function of amplicon simulation is used. More information `here <https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm>`_ .


    Arguments
    ---------
        input_fasta: str
            Reference fasta file
        output_prefix: str
            Output file prefix
        dir_out: str
            Output directory
        art_path: str
            Path to the art software binaries directory
        seq_method: str, default = 'MSv1'
            Sequencing method to simulate
        length: int, default = 200
            Length of the reads to simulate
        count: int, default = 1
            Number of reads per sequence in the fasta file.
        n_samples: int, default = 1
            Number of samples to generate
        mlen: int, default = 400
            Minimum total read length, only applicable if paired = True
        paired: bool, default = False
            If False, simulate single end reads
            If True, simulate paired end reads
        minQ: int, default = 25
            Minimal read quality to simulate


    Returns
    -------
    dict
        Dictionary of fasta paths for each sample {sample:fasta_path}

    References
    ----------
    .. [1] Weichun Huang, Leping Li, Jason R Myers, and Gabor T Marth. ART: a next-generation sequencing read simulator, Bioinformatics (2012) 28 (4): 593-594

    See Also
    --------
    eval_ref()
    
    """


    # Output
    if paired:
        dir_out_prefix = join(dir_out, "paired")
    else:
        dir_out_prefix= join(dir_out, "single")
    check_dir(dir_out_prefix)

    # Simulation
    commands = []

    outfiles_art = []
    # Generate command read generation for each sample
    for i in range(1, n_samples+1):
        # Define outname and commands
        outname = join(dir_out_prefix, f"{output_prefix}S0{i}_sample_L001_R")

        if paired:
            existing_files = [outname + '1.fq']
            existing_files.append(outname + '2.fq')

        else:
            outname = outname + "1"
            existing_files = [outname + '1.fq']

        art_outname = f"{outname}.fq"

        if not all(exists(f) for f in existing_files) or force:
            cmd = f"{art_path}art_illumina -ss {seq_method} -amp -na -nf 0 -i {input_fasta} -l {length} -c {count} -rs {seed} -o {outname} --minQ {minQ}"

            if paired:
                ## add -p parameter and mean size of fragments
                cmd = cmd + f" -p -m {mlen}"


            cmd_split = cmd.split(" ")
            commands.append(cmd_split)

        # Change seed for every sample
        seed = seed + 1

        outfiles_art.append(art_outname)
    # Run all the commands in parallel
    if commands:
        Parallel(n_jobs=10)(delayed(subprocess.run)(cmd_split) for cmd_split in commands);

    return outfiles_art


def adapt_art_output_to_qiime(fileprefix, paired=False, force=False):
    """
    Compress and change name of fq files from art simulation tools.

    In this case, there is one fastq file for each sample (Casava 1.8 format).

    The file name will include the sample identifier and should look like
    vagi_mock01_sample_L001_R1_001.fastq.gz. The underscore-separated fields in
    this file name are:

    - mock01 : the sample identifier (variable),
    - sample : the barcode sequence or a barcode identifier (constant),
    - L001 : the lane number (constant),
    - R1 : the direction of the read (e.g. R1 or R2)
    - 001 : the set number (constant).
    
    Find more info `here <https://docs.qiime2.org/2022.8/tutorials/importing/>`_


    Parameters
    ----------
        fileprefix: str
            file prefix
        paired: bool, default = False
            If reads are single end, set to False. If reads are paired end, set to True

    Returns
    -------
    str
        Adapted file path       


    """
    filenames = []
    if paired:
        filenames.append(fileprefix.replace('_R', '_R1'))
        filenames.append(fileprefix.replace('_R', '_R2'))
    else:
        filenames.append(fileprefix)

    for filename in filenames:
        # Create output directory for fastq.gz files
        art_out_prefix = dirname(filename)
        fq_out_dir = join(art_out_prefix, "fastq")
        check_dir(fq_out_dir)


        # compress fq files and rename it
        fname = basename(filename)
        ## Compressed file ('fq.gz')
        fileout_gz = join(fq_out_dir, f"{fname}.gz")
        ## Compressed file with modified name
        only_name = fileout_gz[:-6] # get file name without extension '.fq.gz'
        new_name = f"{only_name}_001.fastq.gz"


        if not exists(new_name) or force:
            # Compress file to .gz extension
            if not exists(fileout_gz):
                with open(filename, 'rb') as f_in:
                    with gzip.open(fileout_gz, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)
            # Rename
            os.rename(fileout_gz, new_name)


    return new_name

def qiime_demultiplex(out_dir, fq_out_dir, paired=False, force=False):
    """
    Merges samples generated with art simulation tools using qiime2.
    More information in `here <https://docs.qiime2.org/2022.8/tutorials/moving-pictures-usage/>`_.


    Arguments
    ---------
        out_dir: str
            Output directory
        fq_out_dir : str
            Input directory, where the fastq files are located
        paired: bool, default = False
            If reads are single end, set to False. If reads are paired end, set to True
        force: bool, default = False
            If True, overwrites already existing results. 

    Returns
    -------
        demux: qiime2 artifact
            Qiime2 artifact with the demultiplexed sequences
        demux_viz: qiime visualization
            Qiime2 visualization of the demux results
        qiime_out: str
            Output directory where the resulting files are located
    
    """

    # Outputs and parameter
    if paired:
        qiime_out = join(out_dir, "paired/")
        demux_out = join(qiime_out, "demux-paired-end.qza")
        end_type = 'SampleData[PairedEndSequencesWithQuality]'
    else:
        qiime_out = join(out_dir, "single/")
        demux_out = join(qiime_out, "demux-single-end.qza")
        end_type = 'SampleData[SequencesWithQuality]'

    check_dir(qiime_out)

    if not exists(demux_out) or force:
        # Demultiplex
        cmd = f"qiime tools import --type {end_type} --input-path {fq_out_dir} --input-format CasavaOneEightSingleLanePerSampleDirFmt --output-path {demux_out}"
        print(cmd+'\n')
        os.system(cmd)

    # Load results in qiime artifact
    demux = Artifact.load(demux_out)
    demux_viz, = demux_actions.summarize(
        data=demux)
    demux_viz.save(splitext(demux_out)[0])


    return demux, demux_viz, qiime_out

def denoise_dada2(qiime_out,mockname, demux, paired = False,
                  trim_left = 0, trunc_len=0, trunc_len_f = 0, trunc_len_r = 0,
                  trim_left_f = 0, trim_left_r = 0, n_reads_learn = 1000000,
                  threads = 1, force = False):


    """ Dada2 denoising of Illumina reads using qiime2 implementation.
        More information in `here <https://docs.qiime2.org/2022.8/tutorials/moving-pictures-usage/>`_.

    Arguments
    ---------
        qiime_out: str
            Output directory
        mockname: str
            Label for the output files
        demux: qiime2 artifact
            Qiime2 artifact with demultiplexed sequences
        paired: bool, default = False
            If reads are single end, set to False. If reads are paired end, set to True
        trim_left: int, default = 0
            Position at which sequences should be trimmed due to
            low quality. This trims the 5' end of the of the
            input sequences, which will be the bases that were
            sequenced in the first cycles.
            Only applies if paired = False
        trunc_len: int, default = 0
            Position at which sequences should be truncated due
            to decrease in quality. This truncates the 3' end of
            the of the input sequences, which will be the bases
            that were sequenced in the last cycles. Reads that
            are shorter than this value will be discarded. If 0
            is provided, no truncation or length filtering will
            be performed.
            Only applies if paired = False
        trunc_len_f: int, default = 0
            Like trunc_len for the forward reads, only applies if paired = True
        trunc_len_r: int, default = 0
            Like trunc_len for the reverse reads, only applies if paired = True
        trim_left_f: int, default = 0
            Like trim_left for the forward reads, only applies if paired = True
        trim_left_r: int, default = 0
            Like trim_left for the reverse reads, only applies if paired = True
        n_reads_learn: int, default = 1000000
            The number of reads to use when training the error
            model. Smaller numbers will result in a shorter run
            time but a less reliable error model.
        threads: int, default = 1
            Number of processes to use
        force: bool, default = False
            If True, overwrite existing results.   

    Returns
    --------
        table: qiime2 artifact
            Feature table
        rep_seqs: qiime2 artifact
            Representative sequences
        stats: qiime2 artifact
            Stats of the denoising process


    See also
    --------
    qiime_demultiplex() 

    """


    # Outputs
    table_out = join(qiime_out, f"table_{mockname}")
    rep_seqs_out = join(qiime_out, f"rep_seqs_{mockname}")
    stats_out =  join(qiime_out, f"stats_{mockname}")

    if not os.path.exists(f"{table_out}.qza") or not os.path.exists(f"{rep_seqs_out}.qza") or not os.path.exists(f"{stats_out}.qza") or force:
        if paired:
            table, rep_seqs, stats = dada2_actions.denoise_paired(
            demultiplexed_seqs=demux,
            trunc_len_f = trunc_len_f,
            trunc_len_r=trunc_len_r,
            trim_left_f = trim_left_f,
            trim_left_r = trim_left_r,
            n_threads=threads,
            n_reads_learn = n_reads_learn
        )
        else:

            table, rep_seqs, stats = dada2_actions.denoise_single(
                demultiplexed_seqs=demux,
                trim_left=trim_left,
                trunc_len=trunc_len,
                n_threads=threads
            )
        # Save files
        ## TABLE
        table.save(table_out)
        # in biom and tsv format
        unqza(input = f"{table_out}.qza",
            output = f"{table_out}.tsv", output_type='feature_tab')
        ## REP-SEQS
        rep_seqs.save(rep_seqs_out)
        unqza(f"{rep_seqs_out}.qza", f"{rep_seqs_out}.fasta" )
        ## STATS
        stats.save(stats_out)
    else:
        table = Artifact.load(f"{table_out}.qza")
        rep_seqs = Artifact.load(f"{rep_seqs_out}.qza")
        stats = Artifact.load(f"{stats_out}.qza")

    return table, rep_seqs, stats



def generate_taxa_barplots(tax_artifact, tab_artifact, level, output, force = False):

    """ Generate taxa barplot tables from a taxonomy and a feature table artifacts. 
    It has the capability to collapse taxonomies at the desired level

    Arguments
    ---------
        tax_artifact: qiime2 artifact
            Taxonomy
        tab_artifact: qiime2 artifact
            Feature table 
        level: int
            Taxonomical level to collapse taxonomy
        output: str
            Output tsv file path
        force: bool, default = False
            If True, overwrite already existing results


    Return
    -------
    str
        Output file path

    """

    add_species = False

    # Output
    base_out = splitext(output)[0]
    check_dir(dirname(output))

    if not exists(output) or force:
        try:
            # Generate TaxBarplot
            tab_L7, = taxa_actions.collapse(
                table=tab_artifact, taxonomy=tax_artifact, level=level)
        except ValueError:
            tab_L7, = taxa_actions.collapse(
                table=tab_artifact, taxonomy=tax_artifact, level=level-1)

            # if error is risen, add an extra level to the taxonomies
            add_species = True

        tab_L7.save(base_out)

        tsv_tab = utils.unqza(input = f"{base_out}.qza",
                       output = output, output_type = 'feature_tab')

        if add_species:
            tax_df = pd.read_csv(tsv_tab, header = 1, sep = "\t")
            tax_df['#OTU ID'] =  tax_df['#OTU ID']+';__'
            tax_df.columns = pd.MultiIndex.from_tuples(
                zip(['# Constructed from biom file', '', '', '', '', ''], tax_df.columns))
            tax_df.to_csv(tsv_tab, sep= "\t", index=False)

        return tsv_tab
    else:
        return output
