#! ${params.python3}
import sys
sys.path.append('/local/two/Software/python_lib/')
from ETE3_Utils import *
from ETE3_styles import *
from misc_utils import *
import ete3
import pickle

def my_initiate_clades(tree, taxon_clade):
  for l in tree.traverse():
    l.add_feature(pr_name='clade', pr_value=taxon_clade[l.name])

tree = parse_newick("$tree")
set_node_style(tree, node_style_basic)
with open("tax.map", 'rb') as pickled_map:
     taxon_clade_mod = pickle.load(pickled_map)
my_initiate_clades(tree, taxon_clade_mod)

outgroup = tree.get_midpoint_outgroup()
tree.set_outgroup(outgroup)

contigs = get_leaves_by_prefix(tree, 'Contig')
lineage_present = False
for leaf in contigs:
    subtree = leaf.up
    try:
        if subtree.check_monophyly(values=["$params.lineage"], target_attr='clade', ignore_missing=True):
            print('I found something! I think Contig %s in tree %s belongs to %s' % (leaf.name, "$tree", "$params.lineage"))
            lineage_present = True
    except ete3.coretype.tree.TreeError:
        print("clade not found in subtree", file=sys.stderr)
try:
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
except AttributeError:
    print("did not find the attribute clade, because i looked at an internal node", file=sys.stderr)

tree.ladderize(direction=1)
ts = tree_style_basic(layout_node_color, "${tree.baseName}")
tree.render("${tree.baseName}.pdf", tree_style=ts)

with open('decision.txt', 'w') as decision:
    decision.write("\\t".join(["$tree","$params.lineage", str(lineage_present)]) + "\\n")
