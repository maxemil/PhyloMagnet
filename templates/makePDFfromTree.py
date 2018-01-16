#! ${params.python3}
import sys
from ETE3_Utils import *
from ETE3_styles import *
import ete3
import pickle
from xvfbwrapper import Xvfb
import os

vdisplay = Xvfb()
vdisplay.start()

def my_initiate_clades(tree, taxon_clade):
  for l in tree.traverse():
    l.add_feature(pr_name='clade', pr_value=taxon_clade[l.name])

cog_base = "${tree.baseName.minus(~/^.+-/).minus(~/.trim/)}"
for line in open("class.map"):
    line = line.strip().split()
    if line[0] == cog_base:
        os.link("%s_taxid.map" % line[1],  "tax.map")

tree = parse_newick("$tree")
set_node_style(tree, node_style_basic)
with open("tax.map", 'rb') as pickled_map:
     taxon_clade_mod = pickle.load(pickled_map)
my_initiate_clades(tree, taxon_clade_mod)

outgroup = tree.get_midpoint_outgroup()
tree.set_outgroup(outgroup)

contigs = get_leaves_by_prefix(tree, 'Contig')
lineage_present = False
candidate_contigs = set()
for leaf in contigs:
    subtree = leaf.up
    try:
        mono_clade = subtree.check_monophyly(values=["$params.lineage"], target_attr='clade', ignore_missing=True)
        if mono_clade[0]:
            print('I found something! I think Contig %s in tree %s belongs to %s' % (leaf.name, "$tree", "$params.lineage"))
            lineage_present = True
            candidate_contigs.add(leaf.name)
        elif (len(subtree.check_monophyly(values=["$params.lineage"], target_attr='clade', ignore_missing=True)[2]) / len(list(subtree.traverse()))) < 0.15:
            print('I found something! I think Contig %s in tree %s belongs to %s, but its only >85perc monophyletic' % (leaf.name, "$tree", "$params.lineage"))
            lineage_present = True
            candidate_contigs.add(leaf.name)
    except ete3.coretype.tree.TreeError:
        print("clade not found in subtree", file=sys.stderr)
# try:
#     leaves = set(tree.iter_leaves())
#     while leaves:
#       node = leaves.pop()
#       mono_clade = get_mono_clade(node)
#       ancestor = get_ancestor(list(mono_clade))
#       if len(mono_clade) > 1 and node.clade:
#         ancestor.name = node.clade
#         for l in mono_clade:
#           leaves.discard(l)
#         for c in ancestor.get_children():
#           ancestor.remove_child(c)
# except AttributeError:
#     print("did not find the attribute clade, because i looked at an internal node", file=sys.stderr)

tree.ladderize(direction=1)
ts = tree_style_basic(layout_node_color, "${tree.baseName}")
tree.render("${tree.baseName}.pdf", tree_style=ts)

with open('decision.txt', 'w') as decision:
    decision.write("\\t".join(["$tree","$params.lineage", str(lineage_present), ";".join(candidate_contigs)]) + "\\n")

# explicitly kill the xvfb display, because vdisplay.stop() is not working as expected
os.system('(sleep 5 && kill -9 %d) &' % vdisplay.proc.pid)
