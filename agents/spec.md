# Bioinformatics Agent Specification

## Overview
The Bioinformatics Agent is a specialized workflow orchestrator for omics data analysis, supporting GEO dataset download, gene expression analysis, CLIP-seq peak analysis, and related bioinformatics tasks. It intelligently switches between R, bash, and Python to optimize for task type and output format requirements.

---

## Supported Workflows

### 1. GEO Dataset Download & Metadata Extraction
**Input**: GEO accession (e.g., GSE52778)  
**Process**:
- Query NCBI GEO API via R/GEOquery
- Fetch expression matrices and sample metadata
- Filter by inclusion criteria (tissue, treatment, platform)

**Output**:
- Metadata: `{GSE_ID}_metadata.csv` (samples × attributes)
- Expression: `{GSE_ID}_expression.csv` (genes × samples) or `.RData` if >1GB
- Log: `{GSE_ID}_download.log` (date, source, row/col counts)

**Tools**: R (GEOquery), bash (logging)

---

### 2. Gene Expression Analysis Pipeline
**Input**: Expression matrix (RNA-seq counts, microarray intensities), metadata file  
**Process**:
1. **QC Step**: Library size distribution, per-gene mean/variance, PCA
2. **Normalization**: Batch correction if needed; DESeq2/TMM normalization
3. **Differential Expression**: Fit GLM model → logFC, p-value, adjusted p-value per gene
4. **Visualization**: Volcano plot, MA plot, heatmap
5. **Pathway Analysis**: Run g:Profiler via `gprofiler2` on DEG gene sets and compare human/mouse pathways

**Output** (in order of priority):
- DEG table: `{prefix}_DEG_results.csv` (Gene, logFC, log10(p), baseMean, padj)
- Pathway table: `{prefix}_gprofiler_results.csv` (term, source, p-value, padj, genes)
- Pathway plot: `{prefix}_gprofiler_dotplot.png` or `{prefix}_gprofiler_barplot.png`
- QC plots: `{prefix}_qc_distribution.png`, `{prefix}_pca.png`
- Expression matrix (normalized): `{prefix}_normalized.csv` or `.RData`
- Full report: `{prefix}_analysis_log.txt` (steps, thresholds, gene counts, pathway summary)

**Tools**: R (DESeq2, limma, ggplot2, gprofiler2), bash (pipeline orchestration)

---

### 3. CLIP-seq Peak Analysis
**Input**: Peak files (.bed), alignment files (.bam), reference genome  
**Process**:
1. **Peak Merging**: Combine replicates; consolidate overlapping peaks
2. **Annotation**: Map peaks to genes, exons, UTRs, promoters via GenomicRanges
3. **Enrichment**: Motif discovery; pathway/gene-set enrichment using g:Profiler
4. **Integration**: Correlate with expression data (if available)

**Output**:
- Merged peaks: `{protein}_merged_peaks.bed` (chrom, start, end, peak_id, score)
- Annotations: `{protein}_peak_annotations.csv` (peak, gene, feature_type, strand)
- Pathway table: `{protein}_gprofiler_peak_results.csv`
- Motifs: `{protein}_motif_enrichment.txt` (motif, p-value, n_peaks)
- Visualizations: `{protein}_peak_distribution.png`, `{protein}_upset_plot.png`

**Tools**: R (GenomicRanges, rtracklayer, gprofiler2), bash (bedtools, samtools)

---

## Pathway Enrichment & Comparative Analysis
**Task**: Assess whether LIN28-associated gene regulation is ESC-specific or conserved across tissues and species.

**Input**: DEG sets from human GSE39872 and mouse GSE134033, LIN28/CLIP peak-associated gene lists  
**Process**:
1. Run g:Profiler on human and mouse gene sets separately
2. Map mouse gene IDs to human orthologs when needed for comparison
3. Identify conserved pathways, GO terms, and regulatory signatures
4. Compare hESC results to mouse testis and other tissues for shared versus ESC-specific functions

**Output**:
- Cross-species summary: `LIN28_hESC_mouse_testis_pathway_comparison.csv`
- Ortholog mapping: `mouse2human_orthologs.csv`
- Enrichment summary: `LIN28_cross_tissue_pathway_summary.txt`
- Comparison plots: `LIN28_pathway_conservation.png`

**Tools**: R (`gprofiler2`, `biomaRt`, `clusterProfiler`), bash, g:Profiler API

---

## Data Format Rules

### Output Priority
1. **Small tabular data** (<100MB): Export as `.csv` or `.txt`
   - Use `write.csv()` (R) or `write.table()` (R)
   - Include headers and row names where applicable

2. **Plots & Visualizations**: Always `.png` (raster) or `.pdf` (vector)
   - Use `png()`, `pdf()` (R) or `matplotlib.savefig()` (Python)
   - 300 DPI minimum for publication quality

3. **Large intermediate data** (>500MB): Use `.RData` or `.R` script
   - Save as `save(object, file='...RData')`
   - Or create standalone `.R` script to regenerate on demand

4. **Sequence/genomic data**: `.bed`, `.gff`, `.fasta` as standard

### File Naming Convention
```
{dataset_id or protein_name}_{analysis_type}_{date}.{extension}

Examples:
- GSE52778_DEG_results_2026-05-28.csv
- TP53_clip_merged_peaks_2026-05-28.bed
- sample_qc_report_2026-05-28.txt
```

---

## Tool Usage Matrix

| Task | R | Bash | Python | Tool/Package |
|------|---|------|--------|--------------|
| GEO download | ✓ | ✓ | | GEOquery, wget |
| QC & normalization | ✓ | | | DESeq2, limma, edgeR |
| DEG analysis | ✓ | | | DESeq2, limma |
| Visualization | ✓ | | ✓ | ggplot2, matplotlib |
| Peak annotation | ✓ | ✓ | | GenomicRanges, bedtools |
| File conversion | ✓ | ✓ | | samtools, bedtools |
| Batch processing | | ✓ | | grep, awk, parallel |

---

## Error Handling & Validation

### Pre-execution Checks
- Verify GEO accession format before download
- Check file existence and readability before processing
- Validate matrix dimensions (genes > 0, samples > 0)
- Confirm genome reference matches data

### During Execution
- Log all tool versions and parameters
- Check for memory constraints; abort if predicted >available RAM
- Monitor for infinite loops or stalled processes (timeout: 1 hour per task)

### Post-execution Checks
- Verify output files exist and are non-empty
- Validate CSV headers and row counts match expected
- Spot-check plot generation (file size > 50KB typical)
- Report summary statistics: n genes, n DEGs, p-value range

---

## Default Parameters

### Differential Expression
- Fold-change threshold: |log2FC| > 1.0
- Adjusted p-value cutoff: padj < 0.05
- Normalization method: DESeq2 (default) or TMM (edgeR)

### CLIP-seq
- Peak overlap threshold: 50% reciprocal for merging
- Annotation distance: ±2kb for promoter assignment
- Motif E-value: < 1e-5

### Visualization
- PCA: top 5000 variable genes
- Heatmap: top 50 DEGs or all if <50 total
- Plot dimensions: 8×6 inches (300 DPI)

---

## Example Invocation Patterns

1. **Download and analyze GEO dataset**
   ```
   "Download GSE52778, perform DEG analysis, export results as CSV"
   ```

2. **CLIP-seq peak merging and annotation**
   ```
   "Merge TP53 CLIP-seq peaks from 3 replicates, annotate to genes, create distribution plot"
   ```

3. **Custom pipeline**
   ```
   "Normalize expression data, filter low-abundance genes, run DESeq2 with batch correction, output volcano plot"
   ```

---

## Known Limitations

- GEO queries are limited to publicly available datasets
- CLIP-seq analysis requires reference genome annotation (GTF/GFF)
- Large files (>5GB) may require external storage or cloud processing
- Real-time interactive plots (Shiny, plotly) are rendered as static PNG/PDF only
