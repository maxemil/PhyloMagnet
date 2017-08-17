#! ${params.python3}
import pandas as pd
import re
from numpy import mean

df = pd.read_csv("$tree_decisions", sep='\\t',
                names=['tree', 'lineage', 'decision', 'contigs'])
df['sample'] = df.tree.apply(lambda x: re.sub(pattern='-[0-9]*.trim.treefile',
                                            repl='', string=x))
df['decision'] = df['decision'].astype(int)
sample_decisions = df.pivot_table(columns='lineage', index='sample',
                                values='decision', aggfunc=mean)
sample_decisions.sort_values(list(sample_decisions.columns)).to_csv(
                                            'sample_decisions.txt', sep='\\t')

# df = df.fillna('')
# df.contigs = df.contigs.astype(str)
# df.pivot_table(columns=['sample','lineage'],values='contigs', aggfunc=";".join)
