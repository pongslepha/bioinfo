# Bioinformatics Agent Skills

## GEO Data Download
**Task**: Download gene expression datasets from NCBI GEO database.

**Tools**: `GEOquery` (R), `wget`, `curl`

**Key Functions**:
- `getGEO()`: Fetch metadata and expression matrices
- `pData()`, `fData()`: Extract sample/feature annotations
- `exprs()`: Get normalized expression data

**Output Formats**: .csv, .txt (metadata), .RData (large matrices)

---

## Gene Expression Analysis
**Task**: Preprocess, normalize, and perform differential expression analysis on RNA-seq/microarray data.

**Tools**: `DESeq2`, `limma`, `edgeR` (R), `fastqc`, `salmon` (bioinformatics tools)

**Steps**:
1. Quality control (QC): Check library sizes, distributions, batch effects
2. Normalization: TMM, DESeq2, quantile normalization
3. Differential expression: Fit generalized linear models, extract p-values and log2FC
4. Visualization: MA plots, volcano plots, heatmaps

**Output Formats**: 
- DEG results: .csv (genes × logFC, p-value, baseMean)
- QC plots: .png (boxplot, PCA, heatmap)
- Normalized data: .csv (small datasets) or .RData (large)

---

## CLIP-seq Analysis
**Task**: Analyze CLIP-seq (cross-linking immunoprecipitation followed by sequencing) data for RNA-binding protein targets.

**Tools**: `bedtools`, `GenomicRanges` (R), `samtools`, `bwamem`

**Steps**:
1. Download peak files (.bed, .bam) from ENCODE/GEO
2. Merge peaks across replicates (overlap-based or IDR)
3. Annotate genomic features: gene assignment, promoter/exon/UTR
4. Motif discovery and enrichment analysis
5. Integrate with expression data

**Output Formats**:
- Merged peaks: .bed (chrom, start, end, peak_id, score)
- Peak annotations: .csv (peak_id, gene, feature, strand)
- Motif results: .txt (motif ID, p-value, consensus)
- Visualization: .png (peak distribution, genomic context plots)

---

## Data Format Conversion
**Task**: Convert between bioinformatics file formats and optimize storage.

**Tools**: `samtools`, `bedtools`, `rtracklayer` (R)

**Conversions**:
- BAM ↔ SAM
- BED ↔ GTF/GFF
- Matrix (R) → CSV/TSV
- Large .RData files for archival

**Output Formats**: .bed, .gff, .csv, .txt, .RData (as needed)

---

## Pathway & Enrichment Analysis
**Task**: Perform pathway enrichment and functional interpretation for DEG or peak-associated gene sets.

**Tools**: `gprofiler2` (R), `clusterProfiler` (R), g:Profiler web API

**Steps**:
1. Prepare gene lists from differential expression or CLIP peak annotation
2. Query g:Profiler via `gprofiler2` for GO terms, pathways, and regulatory motifs
3. Compare human and mouse pathway results to identify conserved and ESC-specific functions
4. Export tables and generate enrichment plots

**Output Formats**:
- Pathway results: .csv (term, source, p-value, adjusted p-value, genes)
- Enrichment plots: .png (barplot, dotplot, category network)
- Cross-species comparison: .txt/.csv summary tables

---

## Statistical Testing & Visualization
**Task**: Perform statistical tests and generate publication-ready plots.

**Tools**: `ggplot2`, `ComplexHeatmap`, `EnhancedVolcano` (R), `clusterProfiler`, `plotly`, `matplotlib` (Python optional)

**Common Tasks**:
- Volcano plots (log2FC vs -log10(p-value))
- Heatmaps with dendrogram/clustering
- PCA / t-SNE / UMAP plots
- Pathway enrichment plots

**Output Formats**: .png, .pdf (vectorized plots)

---

## Bash Script Execution
**Task**: Execute shell scripts for batch processing, file management, and pipeline automation.

**Capabilities**:
- Download data via `wget`, `curl`, `parallel`
- Process text files with `awk`, `sed`, `grep`
- Manage large files efficiently
- Call external bioinformatics tools (samtools, bedtools, etc.)

**Output Formats**: Flexible based on tool output (text, binary, tabular)
