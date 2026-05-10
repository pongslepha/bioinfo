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


def remove_indels_from_pileup_string(s):
    s = str(s)
    out = []
    i = 0

    while i < len(s):
        ch = s[i]

        # insertion/deletion annotation: +2gc, -1n, -12ACGT...
        if ch in ['+', '-']:
            i += 1

            # 숫자 부분 읽기
            num = ''
            while i < len(s) and s[i].isdigit():
                num += s[i]
                i += 1

            # 숫자만큼 뒤의 sequence 제거
            if num:
                i += int(num)

            continue

        # indel annotation이 아니면 그대로 남김
        out.append(ch)
        i += 1

    return ''.join(out)


pileup['matches'] = pileup['matches'].apply(remove_indels_from_pileup_string)
pileup['matches'] = pileup['matches'].str.upper()

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

print(pileup.sort_values("entropy", ascending=False)[
    ["chrom", "pos", "A", "C", "G", "T", "base_count", "entropy"]
].head(10))


# save to csv

pileup[['chrom', 'pos', "matches", 'A', 'C', 'G', 'T', 'base_count', 'entropy']].to_csv(f'CLIP-{gene}-gene_entrophy.csv', sep='\t', index=False, header=True)



### step 3. bedgraph format ###

bedgraph = pileup[['chrom', 'pos', 'pos', 'entropy']]
bedgraph.columns = ['chrom', 'start', 'end', 'value']
bedgraph['start'] = bedgraph['start'] - 1

# -0.0 제거
bedgraph.loc[bedgraph['value'].abs() < 1e-12, 'value'] = 0

bedgraph_only_nonzer = bedgraph[bedgraph['value'] != 0]
#bedgraph_only_nonzer.to_csv(f'CLIP-{gene}-gene.bedgraph', sep='\t', index=False, header=False)

print(bedgraph_only_nonzer.sort_values("value", ascending=False).head(10))


# UCSC custom track용 bedGraph 저장
with open(f'CLIP-{gene}-gene.ucsc.bedgraph', 'w') as f:
    f.write(
        'track type=bedGraph '
        f'name="CLIP-{gene} entropy" '
        f'description="Shannon entropy of CLIP base composition at Mir{gene}" '
        'visibility=full '
        'autoScale=off '
        'alwaysZero=on '
        'viewLimits=0:2 '
        'graphType=bar '
        'windowingFunction=maximum '
        'maxHeightPixels=60:30:10\n'
    )

    bedgraph_only_nonzer.to_csv(
        f,
        sep='\t',
        header=False,
        index=False,
        float_format='%.6f'
    )



