from ete3 import *
from collections import defaultdict
import warnings

# Avoid annoying warning when importing seaborn
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import seaborn as sns


def parse_newick(infile, mafft=False):
    """
    Parse a newick string/file with ETE3

    :param infile: infile or string
    :param mafft: if True, replace "_R_" prefixes from the sequence labels
    :return: ETE3 tree object
    """
    if mafft:
        tree = Tree(infile)
        for l in tree:
            l.name = l.name.replace('_R_', '')
        return tree
    else:
        return Tree(infile)


def get_leaves_by_prefix(tree, prefix):
    """

    :param tree:
    :param prefix:
    :return:
    """
    nodes = []
    for l in tree.get_leaves():
        if l.name.startswith(prefix):
            nodes.append(l)
    return nodes


def set_node_style(tree, style, leaves=False, condition=None):
    """
    Sets the node style in a tree. If leaves is set as True it only consider
    the leaves. If condition is provided, it only set the style on the nodes
    that meet the conditions
    :param tree: ETE3 Tree Object
    :param style: ETE3 NodeStyle Object
    :param leaves: Boolean
    :param condition: Tuple of three items (operator, node_feature, value)
                      operator: operator function. E.g. lt, gt, eq, le, ge, ne
                      node_feature: node feature to compare
                      value: value to compare with
            Example: (gt, 'counts', 1) will just set the node style in those
            nodes where the counts feature exists and is greater than 1
    :return:
    """
    # Set style for nodes
    for n in tree.traverse():
        # If leaves is True, skip those nodes that are not leaves
        if leaves:
            if n.is_leaf() is False:
                continue

        if condition is None:
            n.img_style = style()
        else:
            operator, node_feature, value = condition
            if node_feature not in n.features:
                continue
            if operator(getattr(n, node_feature), value):
                n.img_style = style()


def defaultdict_string():
    return ""


def defaultdict_array():
    return []


def get_clade_names(taxon_taxid, rank):
    """
    For a number of taxa, get the names correspodnign to the taxid
    at a specified taxonomic level from the ncbi taxonomuy

    :param taxon_taxid: a dictionary containing sequence name to taxonID
    :param rank: desired taxonomic rank, e.g. class, family etc.
    :return: an array containing two dictionaries, clade to taxa and taxa to clade
    """
    taxon_clade = defaultdict(defaultdict_string)
    clade_taxa = defaultdict(defaultdict_array)
    ncbi = ncbi_taxonomy.NCBITaxa()
    for k, v in taxon_taxid.items():
        lin = ncbi.get_lineage(v)
        clade = get_tax_name_ranked(lin, rank)
        taxon_clade[k] = clade
    for k, v in taxon_clade.items():
        clade_taxa[v].append(k)
    return [clade_taxa, taxon_clade]


def get_tax_name_ranked(lineage, rank):
    """
    given a lineage array and a tax level, return the name of the taxon at that level

    :param lineage: ete3 object returned from get_lineage
    :param rank: desired taxonomic rank
    :return: name corresponding to taxid at the specified level
    """
    ncbi = ncbi_taxonomy.NCBITaxa()
    for l in lineage:
        if ncbi.get_rank([l])[l] == rank:
            return ncbi.get_taxid_translator([l])[l]
