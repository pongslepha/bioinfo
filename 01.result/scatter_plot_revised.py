# cd /home/sehwahong/GuidedMission/01.result
# conda activate lab
# python scatter_plot_revised.py

#pip install scipy
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from scipy.stats import pearsonr

# read read-counts.txt
cnts = pd.read_csv('read-counts.txt', sep='\t', comment='#', index_col=0)
print(cnts.columns)

# column 이름 정리: featureCounts가 경로를 붙였을 가능성 대비
cnts.columns = [col.split('/')[-1] for col in cnts.columns]
cnts.head()

# count column만 선택
count_cols = [
    'CLIP-35L33G.bam',
    'RNA-control.bam',
    'RNA-siLin28a.bam',
    'RNA-siLuc.bam',
    'RPF-siLin28a.bam',
    'RPF-siLuc.bam'
]

cnts[count_cols].dtypes

# 숫자형으로 변환
for c in count_cols:
    cnts[c] = pd.to_numeric(cnts[c], errors='coerce')

# library size normalization: CPM
cpm = cnts[count_cols].div(cnts[count_cols].sum(axis=0), axis=1) * 1_000_000
print(cpm.head())




# pseudocount로 0 division 방지
pc = 1

plotdf = pd.DataFrame(index=cnts.index)
plotdf['clip_enrichment_log2'] = np.log2(
    (cpm['CLIP-35L33G.bam'] + pc) / (cpm['RNA-control.bam'] + pc)
)

plotdf['rden_change_log2'] = np.log2(
    ((cpm['RPF-siLin28a.bam'] + pc) / (cpm['RNA-siLin28a.bam'] + pc)) /
    ((cpm['RPF-siLuc.bam'] + pc) / (cpm['RNA-siLuc.bam'] + pc))
)



# filtering 위해 RPKM 계산

# featureCounts length column 사용
# RPKM = reads / (gene_length_kb * total_mapped_reads_million)
rna_siluc_total_million = cnts['RNA-siLuc.bam'].sum() / 1_000_000
gene_length_kb = cnts['Length'] / 1000

cnts['RNA_siLuc_RPKM'] = cnts['RNA-siLuc.bam'] / (gene_length_kb * rna_siluc_total_million)

plotdf['RNA_siLuc_RPKM'] = cnts['RNA_siLuc_RPKM']



# 너무 low count인 gene 제거
plotdf['RNA_control_raw'] = cnts['RNA-control.bam']
plotdf['RPF_siLuc_raw'] = cnts['RPF-siLuc.bam']
plotdf = plotdf[
    (plotdf['RNA_control_raw'] >= 30) &
    (plotdf['RPF_siLuc_raw'] >= 80)
].replace([np.inf, -np.inf], np.nan).dropna()

r, p = pearsonr(plotdf['clip_enrichment_log2'], plotdf['rden_change_log2'])

fig, ax = plt.subplots(figsize=(5, 5))
ax.scatter(
    plotdf['clip_enrichment_log2'],
    plotdf['rden_change_log2'],
    s=3,
    alpha=0.3,
    color='black' 
)

# grid
ax.axhline(0, linewidth=1, color='gray')
ax.axvline(0, linewidth=1, color='gray')

# minor ticks + grid
ax.minorticks_on()
ax.grid(which='major', linestyle='-', linewidth=0.5, alpha=0.5)
ax.grid(which='minor', linestyle='--', linewidth=0.3, alpha=0.3)

# 축 범위: 논문 Figure 4D 느낌
ax.set_xlim(-6, 4)
ax.set_ylim(-2, 2)

# label
ax.set_xlabel('LIN28A CLIP enrichment (log$_2$)')
ax.set_ylabel('Ribosome density change\nupon $Lin28a$ knockdown (log$_2$)')
ax.set_title(f'CLIP and ribosome footprinting upon\n$Lin28a$ knockdown')

ax.text(1, -1.5, f'$r$ = {r:.3f}')

plt.subplots_adjust(left=0.18)
plt.tight_layout()
plt.savefig('scatter_plot_revised.png', dpi=300, bbox_inches='tight')
# plt.show()



# ============================================================
# Protein localization 정보를 이용해서 색깔 있는 plot 다시 그리기
# Figure 5B / S6A-like plot
# ============================================================

# mouselocalization 파일 읽기
mouselocal = pd.read_csv('/home/sehwahong/GuidedMission/00.data/mouselocalization-20210507.txt', sep='\t')

print(mouselocal.head())
print(mouselocal.columns)

# plotdf는 현재 index가 Geneid인 상태이므로 column으로 꺼내기
plotdf_local = plotdf.reset_index()

# reset_index 후 첫 column 이름 확인
print(plotdf_local.columns)

plotdf_local['Geneid_clean'] = plotdf_local['Geneid'].astype(str).str.replace(r'\.\d+$', '', regex=True)
mouselocal['gene_id_clean'] = mouselocal['gene_id'].astype(str).str.replace(r'\.\d+$', '', regex=True)

# mouselocal type 확인
print(mouselocal['type'].value_counts())

# merge
plotdf_local = plotdf_local.merge(
    mouselocal[['gene_id_clean', 'type']],
    left_on='Geneid_clean',
    right_on='gene_id_clean',
    how='left'
)

# merge 결과 확인
print(plotdf_local['type'].value_counts(dropna=False))


# filter out
plotdf_local = plotdf_local[
    (plotdf_local['RNA_siLuc_RPKM'] > 15) & 
    (plotdf_local['type'].isin(['nucleus', 'integral membrane', 'cytoplasm']))
].replace([np.inf, -np.inf], np.nan).dropna()

print(plotdf_local['type'].value_counts())


# 색 지정
color_dict = {
    'nucleus': '#377EB8',             # Set1 blue
    'integral membrane': '#E41A1C',   # Set1 red
    'cytoplasm': '#4DAF4A'            # Set1 green
}

label_dict = {
    'nucleus': 'Nucleus',
    'integral membrane': 'Integral membrane',
    'cytoplasm': 'Cytoplasm'
}


# 새 figure 생성
fig, ax = plt.subplots(figsize=(5.5, 5))

# plotting
for t in ['nucleus', 'integral membrane', 'cytoplasm']:
    sub = plotdf_local[plotdf_local['type'] == t]

    ax.scatter(
        sub['clip_enrichment_log2'],
        sub['rden_change_log2'],
        s=6,
        alpha=0.9,
        color=color_dict[t],
        label=f'{label_dict[t]}'
    )

# 기준선
ax.axhline(0, linewidth=1, color='gray')
ax.axvline(0, linewidth=1, color='gray')

# grid
ax.minorticks_on()
ax.grid(which='major', linestyle='-', linewidth=0.5, alpha=0.5)
ax.grid(which='minor', linestyle='--', linewidth=0.3, alpha=0.3)

# 축 범위: 논문 Figure 5B 느낌
ax.set_xlim(-6, 4)
ax.set_ylim(-2, 3)

# label
ax.set_xlabel('LIN28A CLIP enrichment (log$_2$)')
ax.set_ylabel('Ribosome density change\nupon $Lin28a$ knockdown (log$_2$)')
ax.set_title('Protein localization')

# legend
ax.legend(
    frameon=True,
    fancybox=False,
    edgecolor='black',
    fontsize=10,
    loc='upper left'
)

plt.subplots_adjust(left=0.20)
plt.tight_layout()
plt.savefig('scatter_plot_localization.png', dpi=300, bbox_inches='tight')
# plt.show()




