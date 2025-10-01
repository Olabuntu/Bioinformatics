# RNA-seq Pipeline - Complete Process Breakdown

## 🔍 **Updated Pipeline Flow Analysis (HISAT2 + DGE-Ready Output)**

### **🧬 Quick Biology Refresher:**
- **RNA** comes from DNA and has direction (5' → 3')
- **Genes** can be on either DNA strand (forward or reverse)
- **Stranded libraries** preserve which strand the RNA came from
- **Unstranded libraries** lose this information during preparation
- **Most RNA-seq** uses unstranded libraries (simpler and cheaper)

---

## **PART 1: SETUP & CONFIGURATION (Lines 164-380)**

### **1A. Variable Initialization (Lines 182-216)**
```bash
# INPUT FILES (REQUIRED)
SAMPLE_MANIFEST=""              # Your samples.csv file
REFERENCE_FASTA=""              # Genome file (for alignment mode)
TRANSCRIPTS_FASTA=""            # Transcriptome file (for pseudoalignment mode)
ANNOTATION_FILE=""              # Gene annotation (GTF/GFF3)
WORKFLOW_MODE=""                # ALIGNMENT or PSEUDOALIGNMENT
OUTPUT_DIR=""                   # Where results go

# RESOURCE SETTINGS
THREADS=8                       # CPU cores to use

# PROCESSING OPTIONS
TRIM_READS="true"               # Set to "false" to skip trimming
MIN_READ_LENGTH=50              # Shortest read to keep after trimming
QUALITY_THRESHOLD=20            # Minimum quality score (Q20 = 99% accuracy)
STRANDEDNESS="FR"               # HISAT2 strandedness: F, R, FR, RF
```

**What this does:**
- Sets up all the variables the pipeline needs
- Defines quality thresholds and resource limits
- **STRANDEDNESS**: Controls how reads are counted (F=forward, R=reverse, FR=unstranded, RF=unstranded)

### **1B. Command Line Parsing (Lines 219-330)**
```bash
while [[ $# -gt 0 ]]; do
    case $1 in
        -m|--manifest)
            SAMPLE_MANIFEST="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -w|--workflow)
            WORKFLOW_MODE="$2"
            shift 2
            ;;
        --strandedness)
            STRANDEDNESS="$2"
            shift 2
            ;;
        --no-trim)
            TRIM_READS="false"
            shift
            ;;
```

**What this does:**
- Reads the command line arguments you provide
- Sets the variables based on your flags
- **New**: `--strandedness` parameter for HISAT2 alignment
- **New**: `--no-trim` to skip trimming step

### **1C. Input Validation (Lines 331-380)**
```bash
if [[ -z "$SAMPLE_MANIFEST" ]]; then
    echo "ERROR: Sample manifest required"
    exit 1
fi
```

**What this does:**
- Checks that all required files exist
- Validates file formats and paths
- Ensures you have the right inputs for your chosen workflow

---

## **PART 2: DIRECTORY SETUP (Lines 381-420)**

### **2A. Create Output Structure**
```bash
mkdir -p "$OUTPUT_DIR"/{qc/{raw_fastqc,posttrim_fastqc},trimmed_fastq,indices/{HISAT2,salmon},alignments/{sorted,metrics},quantification,counts_matrix,logs/samples,multiqc,results}
```

**Directory Structure Created:**
```
results/
├── qc/                          ← Quality control reports
│   ├── raw_fastqc/             ← Pre-trimming QC
│   └── posttrim_fastqc/        ← Post-trimming QC
├── trimmed_fastq/              ← Trimmed reads
├── indices/                    ← Reference indices
│   ├── HISAT2/                ← HISAT2 genome index
│   └── salmon/                 ← Salmon transcriptome index
├── alignments/                 ← Alignment results
│   ├── sorted/                 ← Sorted BAM files
│   └── metrics/                ← Alignment statistics
├── quantification/             ← Salmon results (pseudoalignment mode)
├── counts_matrix/              ← Final count matrices
├── logs/                       ← All log files
├── multiqc/                    ← MultiQC reports
└── results/                    ← Final summaries
```

---

## **PART 3: QUALITY CONTROL (Lines 421-520)**

### **3A. Raw Read QC (Lines 421-480)**
```bash
run_raw_qc() {
    echo "Running FastQC on raw reads..."
    for sample in samples; do
        fastqc --outdir "$fastqc_out" --threads "$THREADS" --quiet "$r1_path" "$r2_path"
    done
    multiqc "$fastqc_out" -o "$OUTPUT_DIR/multiqc" -n "raw_fastqc_report"
}
```

**What this does:**
- Runs FastQC on your original FASTQ files
- Generates quality reports for each sample
- Creates a combined MultiQC report
- **Output**: `multiqc/raw_fastqc_report.html`

### **3B. Post-Trim QC (Lines 481-520)**
```bash
run_post_trim_qc() {
    echo "Running FastQC on trimmed reads..."
    for sample in samples; do
        fastqc --outdir "$fastqc_out" --threads "$THREADS" --quiet "$r1_trimmed" "$r2_trimmed"
    done
    multiqc "$fastqc_out" -o "$OUTPUT_DIR/multiqc" -n "posttrim_fastqc_report"
}
```

**What this does:**
- Runs FastQC on trimmed reads
- Compares quality before and after trimming
- **Output**: `multiqc/posttrim_fastqc_report.html`

---

## **PART 4: READ TRIMMING (Lines 521-600)**

### **4A. FastP Trimming (Lines 521-580)**
```bash
run_trimming() {
    echo "Running read trimming with fastp..."
    for sample in samples; do
        fastp \
            --in1 "$r1_path" \
            --in2 "$r2_path" \
            --out1 "$r1_trimmed" \
            --out2 "$r2_trimmed" \
            --length_required "$MIN_READ_LENGTH" \
            --qualified_quality_phred "$QUALITY_THRESHOLD" \
            --html "$trimmed_dir/${sample_id}_fastp.html" \
            --json "$trimmed_dir/${sample_id}_fastp.json" \
            --detect_adapter_for_pe
    done
}
```

**What this does:**
- Trims low-quality bases and adapters
- Removes reads shorter than 50bp
- Keeps only reads with Q20+ quality
- **Auto-detects adapters** for paired-end reads
- **Output**: Trimmed FASTQ files + trimming reports

**Parameters Explained:**
- `--length_required 50`: Keep only reads ≥50bp
- `--qualified_quality_phred 20`: Q20 = 99% accuracy
- `--detect_adapter_for_pe`: Auto-detect adapters

---

## **PART 5: REFERENCE INDEXING (Lines 601-680)**

### **5A. HISAT2 Index Building (Lines 601-650)**
```bash
build_reference_index() {
    if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
        echo "Building HISAT2 genome index..."
        
        # Basic index building
        hisat2-build -p "$THREADS" "$REFERENCE_FASTA" "$hisat2_index_prefix"
        
        # Add annotation file if available (improves alignment accuracy)
        if [[ -f "$ANNOTATION_FILE" ]]; then
            echo "Adding annotation file for improved alignment..."
            # Extract splice sites and exons from GTF/GFF3
            hisat2-build --ss splice_sites.txt --exon exons.txt \
                -p "$THREADS" "$REFERENCE_FASTA" "$hisat2_index_prefix"
        fi
    fi
}
```

**What this does:**
- Builds HISAT2 index from your reference genome
- **Uses annotation file** to extract splice sites and exons
- **Improves alignment accuracy** for RNA-seq data
- **One-time process** - index is reused for all samples
- **Output**: `indices/HISAT2/genome.*.ht2` files + splice site info

### **5B. Salmon Index Building (Lines 651-680)**
```bash
if [[ "$WORKFLOW_MODE" == "PSEUDOALIGNMENT" ]]; then
    echo "Building Salmon transcriptome index..."
    salmon index --transcripts "$TRANSCRIPTS_FASTA" --index "$salmon_index" --threads "$THREADS"
fi
```

**What this does:**
- Builds Salmon index from transcriptome
- **Output**: `indices/salmon/` directory

---

## **PART 6: ALIGNMENT (Lines 681-780)**

### **6A. HISAT2 Alignment (Lines 681-750)**
```bash
run_alignment() {
    echo "Running alignment-based mapping with HISAT2..."
    for sample in samples; do
        hisat2 \
            -x "$hisat2_index" \
            -1 "$r1_input" \
            -2 "$r2_input" \
            -S "$sam_output" \
            --threads "$THREADS" \
            --dta \
            --rna-strandness "$STRANDEDNESS" \
            --max-intronlen 1000000 \
            --min-intronlen 20 \
            --max-seeds 20 \
            --secondary \
            --summary-file "$metrics_dir/${sample_id}_hisat2_summary.txt"
    done
}
```

**What this does:**
- Aligns reads to reference genome using HISAT2
- **Converts SAM to sorted BAM** automatically
- **Indexes BAM files** for downstream analysis
- **Generates alignment statistics**

**Parameters Explained:**
- `--rna-strandness "$STRANDEDNESS"`: Controls strand counting
- `--dta`: Downstream transcript assembly (for StringTie)
- `--max-intronlen 1000000`: Maximum intron length
- `--min-intronlen 20`: Minimum intron length
- `--secondary`: Report secondary alignments

**STRANDEDNESS Options:**

**What does "stranded" vs "unstranded" mean?**

**Unstranded Libraries (`FR` or `RF`):**
- The library preparation method **doesn't preserve strand information**
- You can't tell which strand the RNA originally came from
- **Count reads from both strands** (forward and reverse)
- Most common type of RNA-seq library
- **Use `FR` (default) for most experiments**

**Stranded Libraries (`F` or `R`):**
- The library preparation method **preserves strand information**
- You can tell which strand the RNA originally came from
- **Count only reads from the correct strand**
- More informative but less common
- **Use `F` or `R` only if you know your library is stranded**

**HISAT2 Options:**
- `F`: Forward-stranded (count only forward strand reads)
- `R`: Reverse-stranded (count only reverse strand reads)
- `FR`: Unstranded (count both strands) - **DEFAULT**
- `RF`: Unstranded (count both strands)

---

## **PART 7: EXPRESSION QUANTIFICATION (Lines 781-900)**

### **7A. FeatureCounts Gene Counting (Lines 781-850)**
```bash
run_expression_quantification() {
    echo "Using featureCounts for gene-level quantification..."
    
    # Determine strandedness parameter for featureCounts
    local strand_param=0
    case "$STRANDEDNESS" in
        "F") strand_param=1 ;;
        "R") strand_param=2 ;;
        *) strand_param=0 ;;
    esac
    
    featureCounts \
        -p \
        -T "$THREADS" \
        -t exon \
        -g "$gene_attr" \
        -s "$strand_param" \
        -a "$ANNOTATION_FILE" \
        -o "$counts_dir/raw_gene_counts.txt" \
        "${bam_files[@]}"
}
```

**What this does:**
- Counts reads per gene using featureCounts
- **Handles GFF3 and GTF files** automatically
- **Respects strandedness** settings
- **Output**: `counts_matrix/raw_gene_counts.txt`

### **7B. DGE-Ready Output Generation (Lines 851-900)**
```bash
# Create DGE-ready CSV format
echo "Creating DGE-ready count matrix..."
awk 'NR==1 {print "GeneID," substr($0, index($0, $2))} NR>1 {print $1"," substr($0, index($0, $2))}' "$counts_dir/gene_counts_matrix.tsv" | tr '\t' ',' > "$counts_dir/count_matrix.csv"

# Create sample metadata for DGE analysis
echo "Creating sample metadata..."
echo "sample_id,condition,replicate" > "$counts_dir/sample_metadata.csv"
local replicate=1
for sample in "${sample_names[@]}"; do
    echo "${sample},condition1,${replicate}" >> "$counts_dir/sample_metadata.csv"
    replicate=$((replicate + 1))
done
```

**What this does:**
- **Creates CSV format** for DESeq2/edgeR
- **Generates sample metadata** for DGE analysis
- **Outputs multiple formats** for different tools

---

## **PART 8: FINAL REPORTS (Lines 901-1000)**

### **8A. MultiQC Final Report (Lines 901-950)**
```bash
generate_final_reports() {
    echo "Running final MultiQC report..."
    multiqc "$OUTPUT_DIR" -o "$OUTPUT_DIR/multiqc" -n "multiqc_report_final"
}
```

**What this does:**
- Combines all QC reports into one comprehensive report
- **Output**: `multiqc/multiqc_report_final.html`

### **8B. Output Summary (Lines 951-1000)**
```bash
echo "Key outputs:"
echo "  • Gene counts matrix (TSV): $OUTPUT_DIR/counts_matrix/gene_counts_matrix.tsv"
echo "  • DGE-ready count matrix (CSV): $OUTPUT_DIR/counts_matrix/count_matrix.csv"
echo "  • Sample metadata: $OUTPUT_DIR/counts_matrix/sample_metadata.csv"
echo "  • Final QC report: $OUTPUT_DIR/multiqc/multiqc_report_final.html"
```

---

## **📊 OUTPUT FILES EXPLAINED**

### **Main Outputs:**
1. **`count_matrix.csv`** - DGE-ready count matrix (CSV format)
   - Rows: Genes
   - Columns: Samples
   - Values: Read counts
   - **Ready for DESeq2, edgeR, limma**

2. **`sample_metadata.csv`** - Sample information
   - `sample_id,condition,replicate`
   - **Required for DGE analysis**

3. **`gene_counts_matrix.tsv`** - Original TSV format
   - Same data as CSV, different format
   - **Compatible with most tools**

### **Quality Reports:**
- **`multiqc_report_final.html`** - Comprehensive QC report
- **`*_fastp.html`** - Individual trimming reports
- **`*_hisat2_summary.txt`** - Alignment statistics

### **Intermediate Files:**
- **`*.sorted.bam`** - Sorted, indexed BAM files
- **`*.ht2`** - HISAT2 genome index (reusable)

---

## **🚀 USAGE EXAMPLES**

### **Basic Alignment Workflow:**
```bash
./rna_seq_pipeline.sh \
    -m samples.csv \
    -o results/ \
    -w ALIGNMENT \
    -g genome.fa \
    -a genes.gtf
```

### **With Custom Strandedness:**
```bash
./rna_seq_pipeline.sh \
    -m samples.csv \
    -o results/ \
    -w ALIGNMENT \
    -g genome.fa \
    -a genes.gtf \
    --strandedness F
```

### **Skip Trimming:**
```bash
./rna_seq_pipeline.sh \
    -m samples.csv \
    -o results/ \
    -w ALIGNMENT \
    -g genome.fa \
    -a genes.gtf \
    --no-trim
```

---

## **🔧 PARAMETER REFERENCE**

### **Required Parameters:**
- `-m, --manifest FILE` - Sample manifest (CSV)
- `-o, --output DIR` - Output directory
- `-w, --workflow MODE` - ALIGNMENT or PSEUDOALIGNMENT
- `-g, --genome FILE` - Reference genome (alignment mode)
- `-a, --annotation FILE` - Gene annotation (GTF/GFF3)

### **Optional Parameters:**
- `--threads N` - Number of threads (default: 8)
- `--strandedness STRAND` - F, R, FR, RF (default: FR)
- `--no-trim` - Skip trimming step

### **Internal Parameters:**
- `MIN_READ_LENGTH=50` - Minimum read length after trimming
- `QUALITY_THRESHOLD=20` - Minimum quality score (Q20)
- `STRANDEDNESS="FR"` - Default strandedness setting

---

## **📈 DOWNSTREAM ANALYSIS**

### **R with DESeq2:**
```r
library(DESeq2)
counts <- read.csv('results/counts_matrix/count_matrix.csv', row.names=1)
metadata <- read.csv('results/counts_matrix/sample_metadata.csv')
dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~condition)
dds <- DESeq(dds)
results <- results(dds)
```

### **Python with DESeq2:**
```python
import pandas as pd
counts = pd.read_csv('results/counts_matrix/count_matrix.csv', index_col=0)
metadata = pd.read_csv('results/counts_matrix/sample_metadata.csv')
# Use rpy2 or subprocess to call DESeq2
```

---

## **✅ PIPELINE FEATURES**

### **What's New:**
- ✅ **HISAT2 alignment** (reliable, fast)
- ✅ **DGE-ready output** (CSV + metadata)
- ✅ **Strandedness control** (F, R, FR, RF)
- ✅ **Automatic format detection** (GFF3/GTF)
- ✅ **Comprehensive QC** (MultiQC reports)
- ✅ **Error handling** (detailed logging)

### **What's Removed:**
- ❌ STAR alignment (unreliable)
- ❌ Unused parameters (MEMORY, ADAPTERS_FILE, etc.)
- ❌ Complex configuration options
- ❌ Outdated documentation

---

## **🎯 SUMMARY**

This updated pipeline provides a **streamlined, reliable RNA-seq analysis workflow** that:

1. **Processes your samples** through QC → Trimming → Alignment → Counting
2. **Generates DGE-ready outputs** for differential expression analysis
3. **Handles strandedness** correctly for accurate gene counting
4. **Provides comprehensive QC** reports for quality assessment
5. **Works with both GFF3 and GTF** annotation files
6. **Produces publication-ready** results

**The pipeline is now production-ready and optimized for modern RNA-seq analysis!** 🚀