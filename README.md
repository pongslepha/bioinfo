# bioinfo

### Guided mission 1
script: 
    
    xx.script/scatter_plot_revised.py
    
result:
    
    01.result/scatter_plot_revised.png
    01.result/scatter_plot_localization.png

### Guided mission 2
script: 
    
    xx.script/plot_start_density.R
    
result:
    
    01.result/RPF_density_near_start_codon.png

### Guided mission 3
script: 
    
    xx.script/CIMS.py
    
result:

    01.result/CLIP-let7g-gene.ucsc.bedgraph
    01.result/CLIP-let7d-gene.ucsc.bedgraph
    01.result/CLIP-let7f-1-gene.ucsc.bedgraph
    01.result/hgt_genome_let7g.pdf
    01.result/hgt_genome_let7d.pdf
    01.result/hgt_genome_let7f-1.pdf
    

### Own analysis 1
script:

    # python
    xx.script/01.download_geo.py
    
    # R
    xx.script/01.download_geo.R
    
    
### Own analysis 2
foler:

    # agent setting
    agent/
    

    
### Own analysis 3
foler:

    # agent setting for pathway analysis
    agent/

script:

    xx.script/03.pathway_analysis.R
    xx.script/04.pathway_analysis_GSE134033.R
    xx.script/05.pathway_analysis_GSE39872.R

result:

    # GSE37114 (01.result/own_analysis/GSE37114_analysis/)
    gene_level_CLIP_RPF_localization_table.csv
    localization_group_summary_and_test.csv
    localization_group_bubble_plot.png
    GO_{CC,BP,MF}_enrichment_stats.csv
    GO_{CC,BP,MF}_Figure5A_like.png
    GO_{CC,BP,MF}_dotplot.png

    # GSE134033 (01.result/own_analysis/GSE134033_analysis/)
    gene_level_CLIP_RNA_localization_table.csv
    localization_group_summary_and_test.csv
    localization_group_bubble_plot.png
    GO_{CC,BP,MF}_enrichment_stats.csv
    GO_{CC,BP,MF}_Figure5A_like.png
    GO_{CC,BP,MF}_dotplot.png

    # GSE39872 (01.result/own_analysis/GSE39872_analysis/)
    gene_level_CLIP_RNA_table.csv
    GO_{CC,BP,MF}_enrichment_stats.csv
    GO_{CC,BP,MF}_Figure5A_like.png
    GO_{CC,BP,MF}_dotplot.png

note:

    Same framework as 03 (CLIP enrichment x-axis vs functional change y-axis -> GO Wilcoxon enrichment). Run with the `lab` conda env R.
    - GSE37114: mouse ESC LIN28A (03.pathway_analysis.R). y-axis = ribosome density change on Lin28a knockdown (CLIP + RPF). Localization by ENSMUSG.
    - GSE134033: mouse testes LIN28A. y-axis = Lin28a-KO mRNA fold change (no ribosome profiling). Localization joined by ENSMUSG.
    - GSE39872: human LIN28 (H9 hESC, hg18). org.Hs.eg.db; CLIP cluster BEDs overlapped on hg18 knownGene coords for per-gene CLIP signal
      (00.data/hg18_ucsc/); y-axis = LIN28-KD mRNA fold change; localization dropped (mouse-only file).
    

