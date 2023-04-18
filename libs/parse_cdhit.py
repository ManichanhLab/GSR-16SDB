#!/usr/bin/env python
# coding: utf-8

# ---
# Last update: 08-02-2023
# Authors:   Sara Vega (saravg99@gmail.com)
#            Alejandra Gonzalez (gonmola@hotmail.es)
# ---


import re
import networkx as nx
import matplotlib.pyplot as plt
from os.path import join, basename, splitext, dirname, exists
import pandas as pd
import multiprocessing as mp



def read_cluster(clusterfile):
    """
    Parse cluster file from CD-HIT output
    
    Arguments
    ----------
        clusterfile: str
            Path of the CD-HIT cluster file
    
    Returns
    -------
    list of lists
        List containing the cluster number and the list of all 
        the components of the clusters. The cluster centroid is
        the first element of the component list.
    
    """
        

    # Get all of the sequence IDs from the cluster file
    clusters = []
    counter = -1
    with open(clusterfile, "r") as cfh:
        for line in cfh:
            if line.startswith(">Cluster"):
                counter = counter + 1
                clusters.append(list())
                # id = line.strip()
                # clusters[id] = list()
            else:
                member_id = re.findall(r">(.*)\.{3}", line.strip())

                if not line.strip().endswith("*"):
                    clusters[counter].append(member_id[0])
                # If it is the cluster centroid, append at the beggining of the list
                else:
                    clusters[counter].insert(0, member_id[0])
                    # clusters[id].append(member_id[0])

    return clusters


def cluster_attributes(taxafile):
    """
    Parse taxonomy file. Taxonomy files must be a tsv, with two
    columns: Feature ID and Taxon. Taxon must be seven-lineage taxa
    with lineage prefixes (e.g. k__Bacteria; p__Phylum; ...).
    
    Arguments
    ----------
        taxafile: str
            Path of the taxonomy file
    
    Returns
    -------
        full_taxonomy, species, genus, family: dict
            4 dictionaries (full taxonomy, species, genus and family)
            with ids as keys
    
    """

    taxonomy = {}
    species = {}
    genus = {}
    family = {}
    counter = -1

    with open(taxafile, 'r') as tfh:
        for line in tfh:
            if line.startswith("Feature"):
                pass
            else:
                data = line.strip().split('\t')
                id = data[0]
                taxa = data[1]

                # For each id, save its full taxonomy, family, genus and specie annotation
                taxonomy[id] = taxa
                for i, labels in zip([4,5,6],(family,genus,species)):
                    # Get the family, genus or specie annotation without the prefix
                    val = taxa.split('; ')[i][3:]
                    # If there is no annoation, add a number (enhancing visualization)
                    if val == "" or val == None:
                        counter = counter + 1
                        val = counter
                    labels[id] = val

    return taxonomy, species, genus, family


def create_graph(clusters):
    """
    Create graph from clusters information.
    
    Arguments
    ----------
        clusters: list of lists
            List containing the cluster number and the list of all 
            the components of the clusters. The cluster centroid is
            the first element of the component list.
    
    Returns
    -------
    NetworkX graph
        Graph with all the information of clusters
    
    See Also
    --------
    read_cluster()
    
    """
    # Create graph
    g = nx.Graph()
    
    # Add members
    for nodes in clusters:
        centroid = nodes[0]
        edges = [(centroid, item) for item in nodes[1:]]
        g.add_nodes_from(nodes)
        g.add_edges_from(edges)

        # Save centroid information. The centroid is the first element of node list
        centroids = { nodes[i] : (True if i ==0 else False) for i in range(len(nodes))}
        nx.set_node_attributes(g, centroids, 'centroids')

    return g

def subset_comp(graph, component):
    """
    Generate subgraph from a graph component
    
    Arguments
    ---------
        graph: nx.graph
        component: int
            Component number of the graph
            
    Returns
    -------
    nx.subgraph

    See Also
    --------
    plot_component()
    """
    return graph.subgraph(list(nx.connected_components(graph))[component])



def plot_component(graph, component, labels):
    """
    Plot graph component with rename nodes according to labels
    
    Arguments
    ---------
        graph: nx.graph
        component: int
            component number of the graph
        labels: list
            List with the lineages (e.g. "family", "genus", "species") to be included in the resulting dataframe. Each component will be labeled acoording to the lineage.
            
    Returns
    -------
    None

    See Also
    --------
    get_graph_information()

    """
    # Subset the graph
    GC = subset_comp(graph, component)

    for lab, subplots in zip(labels, [311, 312, 313]):
        # Labels
        newlabels = dict(GC.nodes.data(lab))
        # Rename
        S = nx.relabel_nodes(GC, mapping=newlabels, copy=True)
        # Plot
        plt.subplot(subplots)

        nx.draw(S,
            alpha=0.9, node_color="#009ACD", node_size= 600,
            with_labels=True, font_size=10, font_color="#030303",
            edge_color="#009ACD")

    plt.show()


def get_graph_information(g):
    """
    Obtain nx.graph metrics:
        - Number (#) of components
        - # connected components
        - # components with 1 member
        - # components with > 3 members
        - Histogram showing the distribution of component's length
    
    Arguments
    ---------
        clusters_df: pandas dataframe
           Df where columns are lineage species in taxon and
           rows are clusters
        clade: str
            taxonomic level: family, genus or species
        name: str
            pattern to find. Species taxonomy are separated with '_'.
            
    Returns
    -------
    None

    See Also
    --------
    read_cluster()
    create_graph()
    
    """
    clusters = list(nx.connected_components(g))
    print(f"Number of clusters = {nx.number_connected_components(g)}")
    print(f"Number of cluster with 1 seq = {sum([True for cl in clusters if len(cl) == 1])}")
    print(f"Number of cluster with > 3 seqs = {sum([True for cl in clusters if len(cl) > 3])}\n")

    # Subset graph by the first component
    n_component = 0
    gc = subset_comp(g, n_component )
    print(f"Length of the {n_component} component = {len(gc.nodes)}\n")

    # Histogram showing distribution of component's length
    sizes = [len(g.subgraph(comp).nodes) for comp in clusters ]
    print(f"Sizes of the clusters:\n{set(sizes)}")
    plt.hist(sizes)
    plt.yscale('log')
    plt.show()

def get_relations_by_clade(clusters_df, clade, name):
    """
    Find the clusters where 'name' is present

    Arguments
    ---------
        clusters_df: pandas dataframe
           Df where columns are lineage species in taxon and
           rows are clusters
        clade: str
            taxonomic level: family, genus or species
        name: str
            pattern to find. Species taxonomy are separated with '_'.
            
    Returns
    -------
    dataframe subset:
        df containing the cluster where name is present


    See Also
    --------
    results_to_df2()
    """

    s = clusters_df[clusters_df[clade].str.contains(name)]
    print(f"Number of clusters for {name} {clade}: {len(s.index)}")

    return s

def get_taxa_from_cluster(clu_sub, only_centroid = True):

    """
    Obtain the associated taxa of a cluster object

    Arguments
    ---------
        clu_sub: nx subgraph object
            One of the clusters of cd-hit result
        only_centroid: boolean
            Return only the centroid for each cluster instead of all the ids of the clusters

    Returns
    -------
    tuple
        Two objects:
        1. centroid-id (str) or cluster-ids (list)
        2. dictionary with the clusterized taxa for the input cluster



    """

    taxonomy_list = list(nx.get_node_attributes(clu_sub, "taxonomy").values())

    # Get sets for each clade
    tax = {}
    for item in taxonomy_list:
        tax_split = item.split(';')
        tax_split = [re.sub(r"^[kpcofgsd]__", "", ele.strip().strip('"')) for ele in tax_split]
        tax_split = ["Unknown" if ele == "" else ele for ele in tax_split ]

        for i in [0,1,2,3,4,5,6]:
            tax.setdefault(i, set()).add(tax_split[i])

    db_clu = {}
    # Clusterize taxonomy
    for key, s in tax.items():
        # Order alphabetically
        s = sorted(list(s))

        db_clu.setdefault('taxa', [])

        # If clade contains less than 10 different nomenclatures
        # and the clade is not species
        if len(s) < 10 and key != 6:
            db_clu['taxa'].append('-'.join(s))

        # If clade contains less than 10 different nomenclatures
        # and the clade is species
        elif len(s) < 10 and key == 6:
            names = {}

            # Separate genus and species from nomenclature and
            # Join species from same genus
            for sp in s:
                try:
                    genus = sp.split('_')[0]
                    val = sp.split('_')[1]
                except IndexError:
                    print(sp)
                    break
                names.setdefault(genus, []).append(val)

            taxa_sp = []
            for genus, val in names.items():
                taxa_sp.append(f"{genus}_{'-'.join(val)}")

            db_clu['taxa'].append(':'.join(taxa_sp))
        else:
            db_clu['taxa'].append('Unknown')

    if only_centroid:        
        # Get centroid and return centroid id and clustered taxa
        item = [key for key, val in nx.get_node_attributes(
            clu_sub, 'centroids').items() if val == True]
        centroid_id = item[0]
        return (centroid_id, db_clu)
    else:
        # return all ids (id list) and clustered taxa
        all_ids = list(clu_sub.nodes)
        return (all_ids, db_clu)


def graph_relations_for_comp(comp_num, comp, labels):
    """
    Rename a cluster component according to the labels.
    For each taxon, choose the labels, rename the subgraph and save the edge relations
    
    Arguments
    ---------
        comp_num: int
            Number of the cluster component
        comp: nx subgraph object
            Subgraph (one component) of the main cluster graph.
        labels: list
            List with the lineages (e.g. "family", "genus", "species") to be included in the resulting dataframe. Each cluster will be labeled acoording to the lineage.

    Returns
    -------
        comp_num: int
        relations: list
            List containing the subgraphs rename according to the labels.
            The lenght of the list equals the len of the label list.
    
    """
    # output
    relations = []

    GC = comp

    # For each taxon, choose the labels, rename the subgraph and save the edge relations
    for tax in labels:
        # Labels
        newlabels = dict(GC.nodes.data(tax))
        # Rename
        try:
            S = nx.relabel_nodes(GC, mapping=newlabels, copy=True)
            # Edges relation
            cluster = set()
            for node in set(S.nodes):
                # Unknown clades has been previously converted to numbers
                if isinstance(node, int):
                    cluster.add("Unknown")
                else:
                    cluster.add(node)
            cluster_str = ", ".join(cluster)
            relations.append(cluster_str)

        except ValueError as e:
            pass

    return (comp_num, relations)


def results_to_df2(g, taxon, key, outdir, force=False, threads=1):
    """
    Summarize the information of a NetworkX graph object into
    a dataframe.
    
    Arguments
    ---------
        g: nx subgraph object
            One of the clusters of cd-hit result
        taxon: list
            List with the lineages (e.g. "family", "genus", "species") to be included in the resulting dataframe. Each cluster will be labeled acoording to the lineage.
        outdir: str
            path of the output directory
        force: bool
        threads: int

    Returns
    -------
    dataframe
        Dataframe where columns are lineage species in taxon and
        rows are clusters
    
    """
    relations_out = join(outdir, f"{key}_relations.csv")

    if not exists(relations_out) or force:

        # args = each component of the graph 
        comps = [g.subgraph(comp) for comp in list(nx.connected_components(g))]
        args = [(comps.index(comp), comp, taxon) for comp in comps]

        with mp.Pool(threads) as pool:
            try:
                results = pool.starmap(graph_relations_for_comp, args)
                
                relations = {item[0]:item[1] for item in results}
                clusters_df = pd.DataFrame.from_dict(relations, orient="index", columns = ["family", "genus", "species"])
                clusters_df.to_csv(relations_out)
                return clusters_df

            except KeyboardInterrupt:
                pool.terminate()

    else:
        print(f"{relations_out} exists. Loading csv into Df {time.ctime()}")
        clusters_df = pd.read_csv(relations_out, index_col=0)
        return clusters_df
