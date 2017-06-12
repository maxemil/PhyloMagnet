#! ${params.python3}

from ETE3_Utils import *
from ETE3_styles import *
from misc_utils import *
import ete3
import copy

def collapse_node(node):
    if len(node2labels[node]) == 1 and not collapse_node(node.up):#or not node.clade
        print('new leaf')
        return True
    else:
        return False

def my_initiate_clades(tree, taxon_clade):
    for l in tree.traverse():
        l.add_feature(pr_name='clade', pr_value=taxon_clade[l.name])

tree = parse_newick("SRR3987482-85.trim.treefile")
set_node_style(tree, node_style_basic)
taxon_clade, clade_taxon, taxon_prefix, prefix_taxon = read_prefix_map(tree, "tax.syn")
clade_taxon_mod, taxon_clade_mod = get_clade_names(taxon_clade, 'class')
my_initiate_clades(tree, taxon_clade_mod)

leaves = set(tree.iter_leaves())
while leaves:
  node = leaves.pop()
  mono_clade = get_mono_clade(node)
  ancestor = get_ancestor(list(mono_clade))
  if len(mono_clade) > 1 and node.clade:
    ancestor.name = node.clade
    for l in mono_clade:
      leaves.discard(l)
    for c in ancestor.get_children():
        ancestor.remove_child(c)


tree.ladderize(direction=1)
ts = tree_style_basic(layout_node_color, "test")
tree.render("${tree.baseName}.pdf", w=1500, units="px", tree_style=ts)
