#! /usr/bin/python
# coding: utf-8

# ---
# Last update: 03-02-2023
# Authors:   Sara Vega (saravg99@gmail.com)
#            Alejandra Gonzalez (gonmola@hotmail.es)
# ---

import os
from os.path import join, basename, exists, dirname, splitext
from re import sub
import copy

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



def get_palette_c25():

    """
    Return
    ------
    list
        List of 25 hex color codes
    """

    palette_c25 = ['#1c86ee', '#E31A1C', '#008b00', '#6A3D9A', '#FF7F00', "black",
                  '#ffd700', '#7ec0ee','#FB9A99', '#90ee90', '#CAB2D6', '#FDBF6F',
                  '#b3b3b3', '#eee685', '#800000', '#ff83fa','#ff1493', '#0000ff',
                  '#36648b', '#00ced1', '#00ff00', '#8b8b00', '#cdcd00', '#8b4500',
                   'brown']

    return palette_c25



def get_vregion_primers():
    """Get forward and reverse primers to extract 16S variable regions.

    Primers were extracted from this `article <https://doi.org/10.1128/msphere.01202-20>`_

    Return
    ------
    dict
        Dictionary with forward and reverse primers for 16S variable regions [1]_.

    References
    ----------
    .. [1] Abellan-Schneyder I, Matchado MS, Reitmeier S, Sommer A, Sewald Z, Baumbach J, List M, Neuhaus K. 2021. "Primer, pipelines, parameters: issues in 16S rRNA gene sequencing." mSphere 6:e01202-20., https://doi.org/10.1128/msphere.01202-20


    """
    vregion_dict = {
        'V1-V2': {"f":"AGAGTTTGATYMTGGCTCAG", "r": "GCTGCCTCCCGTAGGAGT"},
        'V1-V3': {"f":"AGAGTTTGATYMTGGCTCAG", "r": "ATTACCGCGGCTGCTGG"},
        'V3-V4': {"f":"CCTACGGGNGGCWGCAG", "r": "GACTACHVGGGTATCTAATCC"},
        'V3-V5': {"f": "CCTACGGGAGGCAGCAG", "r": "CCGTCAATTCMTTTRAGT"},
        'V4': {"f":"GTGCCAGCMGCCGCGGTAA", "r": "GGACTACHVGGGTWTCTAAT"},
        'V4-V5': {"f":"GTGCCAGCMGCCGCGGTAA", "r": "GAATTAAACCACATGCTC"},
        'V6-V8': {"f":"GAATTGACGGGGGCCCGCACAAG", "r": "CGGTGTGTACAAGGCCCGGGAACG"},
        'V7-V9': {"f":"CAACGAGCGCAACCCT", "r": "TACGGYTACCTTGTTACGACTT"}
    }

    return vregion_dict


def get_vregion_lengths():
    """
    Get minimum and maximum length of the v regions

    Return
    ------
    dict
        Dictionary with the min and max length of V regions

    """
    vregion_dict = {
        'V1-V2': {'min': None, 'max': None},
        'V1-V3': {'min': 400, 'max': 600},
        'V3-V4': {'min': 400, 'max': 500},
        'V3-V5': {'min': 450, 'max': 600},
        'V4': {'min': 200, 'max': 350},
        'V4-V5': {'min': None, 'max': None},
        'V6-V8': {'min': None, 'max': None},
        'V7-V9': {'min': None, 'max': None}
    }

    return vregion_dict






def check_dir(path):
    """
    Create ``path`` directories if they do not exist.

    Arguments
    ----------
        path : str
            Path of directories to check and create.

    Returns
    -------
    None

    """

    if not os.path.exists(path):
        os.makedirs(path)




def format_taxonomy_of_taxabarplots(table_df):
    """
    Modify QIIME taxaplots (pandas dataframe) so that the taxonomy syntax matches
    the syntax database. Specifically:

    - add a space between the separator (';') of the clades

    - add the prefix (e.g. g\_\_; s\_\_...) before an unspecified clade

    Note: Taxon must be in the index

    Arguments
    ----------
        table_df: dataframe-like
            Dataframe to format

    Returns
    ----------
    dataframe
        Input dataframe with formatted taxonomy

    """

    # modify taxonomy to match the databases
    table_df.index = table_df.index.str.replace(';__$',';s__', regex=True)
    table_df.index = table_df.index.str.replace(';__;s__$',';g__;s__', regex=True)
    table_df.index = table_df.index.str.replace(';__;g__;s__$',';f__;g__;s__', regex=True)
    table_df.index = table_df.index.str.replace(';__;f__;g__;s__$',';o__;f__;g__;s__', regex=True)
    table_df.index = table_df.index.str.replace(';__;o__;f__;g__;s__$',';c__;o__;f__;g__;s__', regex=True)
    table_df.index = table_df.index.str.replace(';__;c__;o__;f__;g__;s__$',';p__;c__;o__;f__;g__;s__', regex=True)
    table_df.index = table_df.index.str.replace(';','; ', regex=True)

    table_df = table_df.reset_index()
    table_df.rename(columns={"#OTU ID":"Taxon"}, inplace=True)

    return table_df


def get_keys_from_value(d, val):
    """
    Search for a value `val` in a dictionary `d` and return a list of the
    keys that point to that value

    Arguments
    ----------
        d: dictionary
        val: value in dictionary

    Returns
    -------
    list
        List of the keys that point to that value

    """
    return [k for k,v in d.items() if v == val]


def unqza(input, output, output_type = None):
    from biom import load_table

    """
    Extract file from a .qza file

    Arguments
    ----------
        input: str
            .qza file to uncompress
        output: str
            output file.
        output_type: {None, 'feature_tab'}
            If None (default), qza is decompressed to default format
            If 'feature_tab', qza is decompressed to biom format and
                then converted to tsv.


    Returns
    -------
    str
        Output path

    """


    unqza = "/mnt/synology/SHARED_SCRIPTS/extract_qza.sh"


    if output_type:
        if output_type == 'feature_tab':
            base_input = splitext(input)[0]
            # To biom format
            os.system(f"{unqza} {input} {base_input}.biom")
            biom_tab = load_table(f"{base_input}.biom")
            # To tsv format
            if not output.endswith('.tsv'):
                print("BE CAREFUL! YOU HAVE CHOOSE FEATURE TAB BUT YOUR OUTPUT NAME DO NOT CONTAIN THE TSV EXTENSION!")

            with open(output, 'w') as f:
                biom_tab.to_tsv(direct_io=f)
        else:
            print("No output type supported.")
    else:
        os.system(f"{unqza} {input} {output}")

    return output


def file_to_qza(input_file, input_type = None):
    """
    Converts input file to .qza file.
    If no ``input_type`` given, tries to infer from file extension.

    Arguments
    ----------
        input_file: str
            File to import
        input_type: {None, 'taxa', 'seq', 'feature_tab', 'qza'}
            If None (default), it is inferred from file extension

    Returns
    -------
    str
        Output path

    Raises
    ------
    ValueError
        If the provided ``input_type`` is not available and/or it cannot be inferred from the file extension


    See Also
    --------
    import_qiime_artifact: Import a file into a qiime2 artifact object.
    """

    from qiime2 import Artifact

    output = f"{splitext(input_file)[0]}.qza"

    artifact = import_qiime_artifact(input_file, input_type)
    artifact.save(output)

    return output


def import_qiime_artifact(input_file, input_type = None):

    """Imports a file into a qiime2 artifact object.

    Arguments
    ----------
        input_file: str
            File to import
        input_type: {None, 'taxa', 'seq', 'feature_tab', 'qza'}
            If None (default), it is inferred from file extension

    Returns
    -------
    Qiime2 artifact object

    Raises
    ------
    ValueError
        If the provided ``input_type`` is not available and/or it cannot be inferred from the file extension

    """

    from qiime2 import Artifact

    formats = ["taxa", "seq", "feature_tab", "qza"]
    file_ext = {
        '.qza': 'qza',
        '.tsv': 'feature_tab',
        '.txt': 'taxa',
        '.fasta': 'seq'
    }

    if input_type not in formats:
        try:
            input_type = file_ext[splitext(input_file)[1]]
        except KeyError:
            raise ValueError(f'"{input_type}" is not a valid type. Valid types: {formats}')

    if input_type == 'qza':
        artifact = Artifact.load(input_file)

    elif input_type ==  'feature_tab':
        biom_file = f"{splitext(input_file)[0]}.biom"
        os.system(f"biom convert -i {input_file} -o {biom_file} --table-type='Table' --to-hdf5")
        artifact = Artifact.import_data('FeatureTable[Frequency]',biom_file)

    elif input_type == 'taxa':
        artifact = Artifact.import_data('FeatureData[Taxonomy]',input_file)

    elif input_type == 'seq':
        artifact = Artifact.import_data('FeatureData[Sequence]', input_file)

    return artifact





def load_db_from_files(taxa_path=None, seqs_path=None):

    """
    Reads taxonomy file and sequence file from a database and stores it in
    a dictionary with the following format:

    .. code-block:: python

        {
            ID1: {'taxa': [rank1, rank2, rank3, ...], 'seq': Seq_object },
            ID2: ...
        }


    It also deletes entries with uncompleted informations (i.e. only taxa or seq)

    Arguments
    ---------
        taxa_path, seqs_path : str
            Paths of the taxa and seqs files of the database

    Returns
    -------
    dict
        Dictionary with database information


    See Also
    --------
    download_db_from_dict : Download dictionary database

    """

    def load_taxa(taxa_path, database):
        with open(taxa_path, "r") as taxafile:
            for line in taxafile:
                l_split = line.strip().split("\t")
                # taxa id
                ID = l_split[0].strip()

                if ID == 'Feature ID':
                    continue

                # Split Seven-lineages
                ftax = l_split[1].strip()
                ftax_split = rm_prefs(ftax.split(";"))

                database.setdefault(ID, {})
                database[ID]['taxa'] = ftax_split
        return database

    def load_seqs(seqs_path, database):
        with open(seqs_path, "r") as seqsfile:
            for record in SeqIO.parse(seqsfile, "fasta"):
                ID = record.id
                seq= record.seq

                database.setdefault(ID, {})
                database[ID]['seq'] = seq
        return database


    database = {}

    if taxa_path:
        database = load_taxa(taxa_path, database)
    if seqs_path:
        database = load_seqs(seqs_path, database)
    if taxa_path and seqs_path:
        database = delete_uncompleted_entries(database)

    return database



def delete_uncompleted_entries(database):
    """
    Delete entries without sequence or taxa information.

    Arguments
    ---------
        database: dict
            dictionary obtained from function `load_db_from_files()` or with the same format

    Returns
    -------
    dict
        Modified database dictionary

    See Also
    --------
    load_db_from_files : Load database into a dictionary
    download_db_from_dict : Download dictionary database

    """
    n_deleted = 0

    database_copy = copy.deepcopy(database)

    for ID, item in database.items():

        if len(item.keys()) != 2:
            # delete entry
            database_copy.pop(ID)

            # Info
            n_deleted = n_deleted + 1


    database = database_copy

    # Info
    if n_deleted > 0:
        print(f"{n_deleted} entries without sequence or taxa information were deleted")

    return database


def join_taxa_lineages(taxa):

    """
    Join list of clades into formatted taxa string, adding rank prefixes

    Arguments
    ---------
        taxa: list
            List of ranks to join

    Returns
    -------
    str

    See Also
    --------
    rm_prefs : Remove prefixes from a list of clades

    Examples
    -------
    >>> import utils
    >>> taxa = ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacteriales', 'Enterobacteriaceae', 'Salmonella', 'Salmonella_enterica']
    >>> utils.join_taxa_lineages(taxa)
    'k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Enterobacteriales; f__Enterobacteriaceae; g__Salmonella; s__Salmonella_enterica'


    """

    prefixes = ["k__", "p__","c__","o__", "f__", "g__", "s__"]
    taxa_prefs = [''.join(item) for item in list(zip(prefixes, taxa))]
    taxa_full = '; '.join(taxa_prefs)
    return taxa_full



def download_db_from_dict(d, taxa_out_file=None, seqs_out_file=None):

    """
    Downloads database dictionary, creating a taxonomy file and/or a sequences file.
    The function will download only the provided files, e.g. if only the taxa file is provided,
    it will only download the database taxonomy.

    Arguments
    ---------
        d : dict
            database dictionary, as obtained by the function `load_db_from_files()`.
        taxa_out_file: str
            Taxonomy file output path
        seqs_out_file: str
            Sequence file output path

    Returns
    -------
    None

    Raises
    --------
    Exception
        If missing both file paths or missing the required taxa/sequence information in the database dictionary.


    See Also
    --------
    load_db_from_files : Load database files into a dictionary

    """


    records = []


    if taxa_out_file:
        with open(taxa_out_file, 'w') as taxa_out:
            # Add header to taxa file
            taxa_out.write("Feature ID\tTaxon\n")

            for ID, item in d.items():
                if taxa_out_file and 'taxa' in item:
                    taxa_full = join_taxa_lineages(taxa = item['taxa'])
                    taxa_out.write(f"{ID}\t{taxa_full}\n")
                elif taxa_out_file and 'taxa' not in item:
                    raise Exception(f"This dictionary doesn't contain taxa information. ID = {ID}")

                # Seqs file
                if seqs_out_file and 'seq' in item:
                    records.append(SeqRecord(item['seq'], id=ID, description=""))
                elif seqs_out_file and 'seq' not in item:
                    raise Exception(f"This dictionary doesn't contain sequence information. ID = {ID}")

    elif seqs_out_file:
        for ID, item in d.items():
            if seqs_out_file and 'seq' in item:
                    records.append(SeqRecord(item['seq'], id=ID, description=""))
            elif seqs_out_file and 'seq' not in item:
                raise Exception(f"This dictionary doesn't contain sequence information. ID = {ID}")


    else:
        raise Exception(f'No path provided.')

    # Write Seqs file
    if records:
        SeqIO.write(records, seqs_out_file,'fasta')



def rm_prefs(it, replace_unknowns = False):
    """
    Removes taxonomy prefixes from all the elements of the list

    Arguments
    ---------
        it: list or set or pandas dataframe column
            Taxonomy clades to modify
        replace_unknowns: bool, default = True
            If True, replaces undefined taxon levels with 'Unknown'.
            If False (default), those levels will become empty strings.

    Returns
    -------
    list
        Taxa without prefixes


    See Also
    --------
    join_taxa_lineages : join list of clades into formatted taxa string

    Examples
    --------
    >>> import utils
    >>> taxa = ['k__Bacteria', 'p__Proteobacteria', 'c__Gammaproteobacteria', 'o__Enterobacteriales', 'f__Enterobacteriaceae', 'g__Salmonella', 's__']
    >>> utils.rm_prefs(taxa)
    ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacteriales', 'Enterobacteriaceae', 'Salmonella', '']

    >>> utils.rm_prefs(taxa, replace_unknowns = True)
    ['Bacteria', 'Proteobacteria', 'Gammaproteobacteria', 'Enterobacteriales', 'Enterobacteriaceae', 'Salmonella', 'Unknown']




    """
    result = [sub(r"^[kpcofgsd]__", "", item.strip().strip('"')) for item in it]

    if replace_unknowns:
        result = [i if i!= "" else "Unknown" for i in result]


    return result


def get_ids_from_taxafile(taxapath):

    """
    Read taxonomy file and extract the entry IDs

    Arguments
    ---------
        taxapath: str
            Path of the taxonomy file to read

    Returns
    -------
    list
        List of IDs contained in the taxa file
    """

    idlist = []
    with open(taxapath, "r") as taxafile:
        for line in taxafile:
            l_split = line.strip().split("\t")
            # taxa id
            ID = l_split[0].strip()

            if ID == 'Feature ID':
                continue
            else:
                idlist.append(ID)

    return idlist

def filter_patterns(df, exclude_dbs=None, exclude_pgram=None, exclude_conf=None):
    """
    Filter databases, pgrams or confidence levels from per_sequence_score_obsid dataframes.
    If no optional argument is provided, no filtering is applied.

    Arguments
    ---------
        df: dataframe (required)
            dataframe to be filtered
        exclude_dbs: list (optional)
            patterns of databases to be excluded
        exclude_pgram: list (optional)
            pgrams to be excluded
        exclude_conf: list (optional)
            confidence levels to be excluded

    Returns
    -------
    list
        List of IDs contained in the taxa file
    """
    if exclude_dbs:
        df = df[
            df['database'].apply(lambda db: all(pat not in db for pat in exclude_dbs)
                                )].reset_index(drop=True)
    if exclude_pgram:
        df = df[
            df['pgram'].apply(lambda pg: all(pgram not in pg for pgram in exclude_pgram)
                             )].reset_index(drop=True)
    if exclude_conf:
        df = df[
            df['conf'].apply(lambda conf: all(
                str(confidence) not in conf for confidence in exclude_conf)
                            )].reset_index(drop=True)

    return df

################################################################################
## Benchmarking
################################################################################

import time

class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""

class Timer:
    timers = {}

    def __init__(
        self,
        name=None,
        text="Elapsed time: {} {:0.4f} seconds",
        logger=print
    ):

        self._start_time = None
        self.name = name
        self.text = text
        self.logger = logger


        # Add new named timers to dictionary of timers
        if name:
            self.timers.setdefault(name, 0)

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None

        if self.logger:
            self.logger(self.text.format(self.name, elapsed_time))

        if self.name:
            self.timers[self.name] += elapsed_time

        return elapsed_time

    def reset_timers(self):
        self.timers = {}
