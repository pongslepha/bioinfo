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

script:

    xx.script/03.pathway_analysis.R
    xx.script/04.pathway_analysis_GSE134033.R
    xx.script/05.pathway_analysis_GSE39872.R
    xx.script/06.filter_GO_dotplot_top15.R
    xx.script/07.comparison_between_dataset.R

result:

    # GSE37114 (01.result/own_analysis/GSE37114_analysis/)
    gene_level_CLIP_RPF_localization_table.csv
    localization_group_summary_and_test.csv
    localization_group_bubble_plot.png
    GO_{CC,BP,MF,ALL}_enrichment_stats.csv
    GO_{CC,BP,MF,ALL}_Figure5A_like.png
    GO_{CC,BP,MF,ALL}_dotplot.png
    GO_{CC,BP,MF,ALL}_dotplot_top15.csv
    GO_combined_dotplot_top15.csv

    # GSE134033 (01.result/own_analysis/GSE134033_analysis/)
    gene_level_CLIP_RNA_localization_table.csv
    localization_group_summary_and_test.csv
    localization_group_bubble_plot.png
    GO_{CC,BP,MF,ALL}_enrichment_stats.csv
    GO_{CC,BP,MF,ALL}_Figure5A_like.png
    GO_{CC,BP,MF,ALL}_dotplot.png
    GO_{CC,BP,MF,ALL}_dotplot_top15.csv
    GO_combined_dotplot_top15.csv

    # GSE39872 (01.result/own_analysis/GSE39872_analysis/)
    gene_level_CLIP_RNA_table.csv
    GO_{CC,BP,MF,ALL}_enrichment_stats.csv
    GO_{CC,BP,MF,ALL}_Figure5A_like.png
    GO_{CC,BP,MF,ALL}_dotplot.png
    GO_{CC,BP,MF,ALL}_dotplot_top15.csv
    GO_combined_dotplot_top15.csv

    # comparison (01.result/own_analysis/comparison_between_datasets/)
    GO_CLIP_enrichment_comparison_top5.png
    GO_CLIP_enrichment_comparison_top5.csv

