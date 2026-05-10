import sys
import os

# Command line argument is expected to be the gene name
gene = sys.argv[1]

dir = '/home/sehwahong/GuidedMission/01.result'
os.chdir(dir)
print(os.getcwd())

import pandas as pd
import numpy as np

pileup = pd.read_csv(f'CLIP-{gene}-gene.pileup', sep='\t', names=['chrom', 'pos', 'ref', 'count', 'basereads', 'quals'])
pileup.tail()

pileup['matches'] = pileup['basereads'].str.replace(r'\^.', '', regex=True)
pileup['matches'] = pileup['matches'].str.replace(r'[<>$*#]', '', regex=True)

pileup[['chrom', 'pos', 'matches']]

### step 1. base count 

pileup['A'] = pileup['matches'].str.count('A')
pileup['G'] = pileup['matches'].str.count('G')
pileup['C'] = pileup['matches'].str.count('C')
pileup['T'] = pileup['matches'].str.count('T')

pileup['base_count'] = pileup[['A', 'G', 'C', 'T']].sum(axis=1)

pileup[['chrom', 'pos', 'A', 'G', 'C', 'T', 'base_count']]

# check: target position = 106056094
#pileup[pileup['pos'] == 106056094].iloc[0][["chrom", "pos", "count", "matches", "A", "C", "G", "T", "base_count"]]

# remove base_count == 0 
pileup = pileup[pileup['base_count'] > 0]

### step 2. Shannon entropy #####

pileup['entropy'] = - pileup[['A', 'G', 'C', 'T']].apply(lambda x: ((x[x > 0] / x.sum()) * np.log2(x[x > 0] / x.sum())).sum(), axis=1)

pileup[['chrom', 'pos', 'entropy']]
#pileup[pileup['pos'] == 106056094].iloc[0][["chrom", "pos", "count", "matches", "A", "C", "G", "T", "base_count", "entropy"]]

pileup.sort_values("entropy", ascending=False)[
    ["chrom", "pos", "A", "C", "G", "T", "base_count", "entropy"]
].head(20)

# save to csv

pileup[['chrom', 'pos', "matches", 'A', 'C', 'G', 'T', 'base_count', 'entropy']].to_csv(f'CLIP-{gene}-gene_entrophy.csv', sep='\t', index=False, header=True)



### step 3. bedgraph format ###

bedgraph = pileup[['chrom', 'pos', 'pos', 'entropy']]
bedgraph.columns = ['chrom', 'start', 'end', 'value']
bedgraph['start'] = bedgraph['start'] - 1

bedgraph.to_csv(f'CLIP-{gene}-gene.bedgraph', sep='\t', index=False, header=False)



