#! ${params.python3}
import sys
sys.path.append('/local/two/Software/python_lib/')
from ETE3_Utils import *
from ETE3_styles import *
from misc_utils import *
import ete3

def my_initiate_clades(tree, taxon_clade):
  for l in tree.traverse():
    l.add_feature(pr_name='clade', pr_value=taxon_clade[l.name])

tree = parse_newick("$tree")
set_node_style(tree, node_style_basic)
taxon_clade, clade_taxon, taxon_prefix, prefix_taxon = read_prefix_map(tree, "$tax_map_concat")
clade_taxon_mod, taxon_clade_mod = get_clade_names(taxon_clade, "$params.taxonomy_level_trees")
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
ts = tree_style_basic(layout_node_color, "${tree.baseName}")
tree.render("${tree.baseName}.pdf", tree_style=ts)#w=1500, units="px",
