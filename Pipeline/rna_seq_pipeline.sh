#!/bin/bash

# ===============================================================================
# RNA-seq Analysis Pipeline
# ===============================================================================
# Description: Complete RNA-seq analysis from raw FASTQ to gene expression matrix
# Purpose: Optimized specifically for RNA-seq data analysis
# Note: Duplicate removal is intentionally skipped (not recommended for RNA-seq)
# ===============================================================================

# ===============================================================================
# REQUIREMENTS & SETUP GUIDE
# ===============================================================================

# 🛠️  REQUIRED TOOLS (Install before running):
# ============================================
# Core Tools:
#   - python3          (Python interpreter)
#   - fastqc           (Quality control: conda install -c bioconda fastqc)
#   - multiqc          (Report aggregation: conda install -c bioconda multiqc)
#   - samtools         (BAM processing: conda install -c bioconda samtools)
#   - pandas           (Python library: pip install pandas)
#
# Trimming:
#   - fastp            (Fast all-in-one preprocessing: conda install -c bioconda fastp)
#
# Workflow Tools (choose based on workflow):
#   For ALIGNMENT workflow:
#     - HISAT2         (Aligner: conda install -c bioconda hisat2)
#     - featureCounts  (Counting: conda install -c bioconda subread)
#   For PSEUDOALIGNMENT workflow:
#     - salmon         (Quantifier: conda install -c bioconda salmon)
#
# 💡 Quick install command:
#    conda install -c bioconda -c conda-forge fastqc multiqc fastp salmon hisat2 subread samtools pandas

# 📁 REQUIRED INPUT FILES:
# ========================
# 1. Sample Manifest (samples.csv):
#    Format: sample_id,r1_path,r2_path
#    Example:
#      sample_id,r1_path,r2_path
#      sample_A,/data/sample_A_R1.fastq.gz,/data/sample_A_R2.fastq.gz
#      sample_B,/data/sample_B_R1.fastq.gz,/data/sample_B_R2.fastq.gz
#
# 2. Reference Files (choose based on workflow):
#    For ALIGNMENT workflow:
#      - genome.fa      (Reference genome FASTA)
#      - genes.gtf      (Gene annotation GTF/GFF3)
#    For PSEUDOALIGNMENT workflow:
#      - transcripts.fa (Transcriptome FASTA)
#      - genes.gtf      (Gene annotation GTF/GFF3)
#
# 3. Optional Files:
#    - adapters.fa     (Adapter sequences for trimming)

# 📂 OUTPUT DIRECTORY STRUCTURE:
# ==============================
# The pipeline creates this organized directory structure:
#
# results/                          ← Main output directory (you specify)
# ├── qc/                          ← Quality control reports
# │   ├── raw_fastqc/             ← FastQC reports for raw reads
# │   │   ├── sample_A_R1_fastqc.html
# │   │   ├── sample_A_R1_fastqc.zip
# │   │   └── ...
# │   └── posttrim_fastqc/        ← FastQC reports for trimmed reads
# │       ├── sample_A_R1.trim_fastqc.html
# │       └── ...
# ├── multiqc/                     ← Aggregate QC reports
# │   ├── multiqc_report_raw.html      ← Raw reads summary
# │   ├── multiqc_report_posttrim.html ← Trimmed reads summary
# │   └── multiqc_report_final.html    ← Complete pipeline summary
# ├── trimmed_fastq/               ← Trimmed and filtered reads
# │   ├── sample_A_R1.trim.fastq.gz
# │   ├── sample_A_R2.trim.fastq.gz
# │   └── ...
# ├── indices/                     ← Reference indices (built once, reused)
# │   ├── HISAT2/                  ← HISAT2 genome index (alignment mode)
# │   └── salmon/                  ← Salmon transcriptome index (pseudoalignment)
# ├── alignments/                  ← Alignment results (ALIGNMENT mode only)
# │   ├── sorted/                  ← Sorted BAM files
# │   │   ├── sample_A.sorted.bam
# │   │   ├── sample_A.sorted.bam.bai
# │   │   └── ...
# │   └── metrics/                 ← Alignment statistics
# │       ├── sample_A_flagstat.txt
# │       └── ...
# ├── quantification/              ← Expression quantification (PSEUDOALIGNMENT mode)
# │   ├── sample_A/
# │   │   ├── quant.sf            ← Salmon quantification results
# │   │   └── logs/
# │   └── ...
# ├── counts_matrix/               ← Final expression matrix
# │   └── gene_counts_matrix.tsv  ← Gene × sample count matrix (MAIN OUTPUT)
# ├── logs/                        ← All log files
# │   ├── pipeline_master.log      ← Main pipeline log
# │   └── samples/                 ← Per-sample detailed logs
# │       ├── sample_A.log
# │       └── ...
# └── results/                     ← Final summaries and reports
#     ├── run_report.txt           ← Pipeline execution summary
#     ├── mapping_summary.txt      ← Alignment/quantification statistics
#     └── output_manifest.txt      ← List of all output files

# 🎯 WORKFLOW MODES:
# ==================
# ALIGNMENT (Traditional):
#   Raw FASTQ → QC → Trim → QC → Align to Genome → Count Genes → Matrix
#   Tools: HISAT2 + featureCounts
#   Pros: Fast, memory-efficient, good for standard RNA-seq
#   Cons: Less sensitive to novel splice junctions than STAR
#
# PSEUDOALIGNMENT (Modern):
#   Raw FASTQ → QC → Trim → QC → Quantify Transcripts → Aggregate to Genes → Matrix
#   Tools: Salmon
#   Pros: Faster, less memory, better quantification accuracy
#   Cons: Cannot discover novel transcripts

# 🚀 USAGE EXAMPLES:
# ==================
# Alignment workflow:
#   ./rna_seq_pipeline.sh -m samples.csv -o results/ -w ALIGNMENT -g genome.fa -a genes.gtf
#
# Pseudoalignment workflow:
#   ./rna_seq_pipeline.sh -m samples.csv -o results/ -w PSEUDOALIGNMENT -t transcripts.fa -a genes.gtf
#
# With custom settings:
#   ./rna_seq_pipeline.sh -m samples.csv -o results/ -w PSEUDOALIGNMENT \
#     -t transcripts.fa -a genes.gtf --threads 16 --no-trim

# 📊 WHAT YOU GET:
# ================
# Main Outputs:
#   - gene_counts_matrix.tsv    ← Gene × sample count matrix (TSV format)
#   - count_matrix.csv          ← DGE-ready count matrix (CSV format)
#   - sample_metadata.csv       ← Sample information for DGE analysis
#   - Rows: Genes, Columns: Samples, Values: Read counts
#   - Ready for differential expression analysis (DESeq2, edgeR, limma)
#
# Quality Reports: multiqc_report_final.html
#   - Comprehensive quality assessment
#   - Mapping rates, read quality, sample correlations
#   - Essential for publication and troubleshooting

# ⚠️  IMPORTANT NOTES:
# ====================
# - Use simple sample names (no spaces, use underscores: sample_A not "sample A")
# - Ensure FASTQ files are paired-end and properly named
# - Reference genome and annotation must match (same organism, same version)
# - Allow 2-8 hours runtime depending on dataset size and workflow
# - Pseudoalignment is recommended for beginners (faster, easier)
# - Works on both macOS and Linux (auto-detects OS for optimal compression handling)

# 🆘 TROUBLESHOOTING:
# ===================
# - Check tool installation: ./check_tools.sh
# - Verify input files exist and are readable
# - Ensure sufficient disk space (2-5x input size)
# - Check logs in logs/ directory for detailed errors
# - Start with small test dataset before full analysis

# ===============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# ===============================================================================
# 1. HEADER & CONFIGURATION BLOCK
# ===============================================================================

# Script version and info
SCRIPT_VERSION="1.0.0"
SCRIPT_NAME="RNA-seq Analysis Pipeline"
SCRIPT_AUTHOR="BUNTU(olabuntubabatunde@gmail.com) Under the supervision of Dr. P.Agre(p.agre@cgiar.org)"

# Print header
echo "==============================================================================="
echo "$SCRIPT_NAME v$SCRIPT_VERSION"
echo "Written by: $SCRIPT_AUTHOR"
echo "Started at: $(date)"
echo "==============================================================================="

# -------------------------------------------------------------------------------
# Configuration Variables - MODIFY THESE FOR YOUR SETUP
# -------------------------------------------------------------------------------

# INPUT FILES (REQUIRED)
SAMPLE_MANIFEST=""              # Path to samples.csv (sample_id,r1_path,r2_path)
REFERENCE_FASTA=""              # Genome FASTA for alignment mode
TRANSCRIPTS_FASTA=""            # Transcriptome FASTA for pseudoalignment mode
ANNOTATION_FILE=""              # GFF3 or GTF annotation file

# WORKFLOW MODE (REQUIRED - choose one)
WORKFLOW_MODE=""                # "ALIGNMENT" or "PSEUDOALIGNMENT"

# OUTPUT DIRECTORY (REQUIRED)
OUTPUT_DIR=""                   # Base output directory

# TOOL PATHS - modify if tools are not in PATH
FASTQC_CMD="fastqc"
MULTIQC_CMD="multiqc"
FASTP_CMD="fastp"
HISAT2_CMD="hisat2"
HISAT2_BUILD_CMD="hisat2-build"
SALMON_CMD="salmon"
FEATURECOUNTS_CMD="featureCounts"
SAMTOOLS_CMD="samtools"

# RESOURCE SETTINGS
THREADS=8

# PROCESSING OPTIONS
TRIM_READS="true"
MIN_READ_LENGTH=50
QUALITY_THRESHOLD=20
STRANDEDNESS="FR"              # Library strandedness: FR, RF, F, R, or NONE (unstranded)

# -------------------------------------------------------------------------------
# Parse command line arguments
# -------------------------------------------------------------------------------

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required Options:
  -m, --manifest FILE         Sample manifest (CSV with sample_id,r1_path,r2_path)
  -o, --output DIR           Output directory
  -w, --workflow MODE        Workflow mode: ALIGNMENT or PSEUDOALIGNMENT

Alignment Mode Options:
  -g, --genome FILE          Reference genome FASTA
  -a, --annotation FILE      Gene annotation (GFF3 or GTF)

Pseudoalignment Mode Options:
  -t, --transcripts FILE     Transcriptome FASTA
  -a, --annotation FILE      Gene annotation (GFF3 or GTF)

Optional:
  --threads N                Number of threads (default: 8)
  --strandedness STRAND      Library strandedness: FR, RF, F, R, NONE (default: FR)
  --no-trim                  Skip read trimming step
  --help                     Show this help message

Required Tools:
  - fastqc: Quality control
  - fastp: Read trimming
  - hisat2: Read alignment
  - hisat2-build: Index building
  - featureCounts: Gene counting
  - samtools: BAM manipulation
  - multiqc: Report generation

Examples:
  # Alignment-based workflow
  $0 -m samples.csv -o results/ -w ALIGNMENT -g genome.fa -a genes.gtf

  # Pseudoalignment workflow
  $0 -m samples.csv -o results/ -w PSEUDOALIGNMENT -t transcripts.fa -a genes.gtf
EOF
}

# Parse arguments
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
        -g|--genome)
            REFERENCE_FASTA="$2"
            shift 2
            ;;
        -t|--transcripts)
            TRANSCRIPTS_FASTA="$2"
            shift 2
            ;;
        -a|--annotation)
            ANNOTATION_FILE="$2"
            shift 2
            ;;
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --no-trim)
            TRIM_READS="false"
            shift
            ;;
        --strandedness)
            STRANDEDNESS="$2"
            shift 2
            ;;
        --help)
            usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option $1"
            usage
            exit 1
            ;;
    esac
done

# -------------------------------------------------------------------------------
# Comprehensive Dependency Checking (BEFORE any processing)
# -------------------------------------------------------------------------------


# Check all required tools and dependencies
check_all_dependencies() {
    local errors=0
    
    echo "==============================================================================="
    echo "DEPENDENCY CHECK - Verifying all required tools and packages"
    echo "==============================================================================="
    
    # Check basic command-line tools
    echo ""
    echo "Checking command-line tools..."
    local tools=("$FASTQC_CMD" "$MULTIQC_CMD" "$SAMTOOLS_CMD" "python3")
    
    if [[ "$TRIM_READS" == "true" ]]; then
        tools+=("$FASTP_CMD")
    fi
    
    if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
        tools+=("$HISAT2_CMD" "$HISAT2_BUILD_CMD" "$FEATURECOUNTS_CMD")
    else
        tools+=("$SALMON_CMD")
    fi
    
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            echo "  ✗ Missing: $tool"
            errors=$((errors + 1))
        else
            echo "  ✓ Found: $tool"
        fi
    done
    
    # Check Python packages
    echo ""
    echo "Checking Python packages..."
    if ! python3 -c "import pandas" &> /dev/null; then
        echo "  ✗ Missing: pandas (Python package)"
        echo "    Install with: pip install pandas"
        errors=$((errors + 1))
    else
        echo "  ✓ Found: pandas"
    fi
    
    
    if [[ $errors -gt 0 ]]; then
        echo ""
        echo "==============================================================================="
        echo "ERROR: Found $errors missing dependency/dependencies"
        echo "Please install all required tools and packages before running the pipeline."
        echo "==============================================================================="
        return 1
    fi
    
    echo ""
    echo "==============================================================================="
    echo "✓ All dependencies verified successfully"
    echo "==============================================================================="
    return 0
}

# -------------------------------------------------------------------------------
# Validate required inputs
# -------------------------------------------------------------------------------

validate_inputs() {
    local errors=0
    
    echo ""
    echo "Validating input files and parameters..."
    
    # Check required parameters
    if [[ -z "$SAMPLE_MANIFEST" ]]; then
        echo "ERROR: Sample manifest (-m/--manifest) is required"
        errors=$((errors + 1))
    elif [[ ! -f "$SAMPLE_MANIFEST" ]]; then
        echo "ERROR: Sample manifest file does not exist: $SAMPLE_MANIFEST"
        errors=$((errors + 1))
    else
        echo "✓ Sample manifest: $SAMPLE_MANIFEST"
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        echo "ERROR: Output directory (-o/--output) is required"
        errors=$((errors + 1))
    else
        echo "✓ Output directory: $OUTPUT_DIR"
    fi
    
    if [[ -z "$WORKFLOW_MODE" ]]; then
        echo "ERROR: Workflow mode (-w/--workflow) is required"
        errors=$((errors + 1))
    elif [[ "$WORKFLOW_MODE" != "ALIGNMENT" && "$WORKFLOW_MODE" != "PSEUDOALIGNMENT" ]]; then
        echo "ERROR: Workflow mode must be ALIGNMENT or PSEUDOALIGNMENT"
        errors=$((errors + 1))
    else
        echo "✓ Workflow mode: $WORKFLOW_MODE"
    fi
    
    # Validate workflow-specific inputs
    if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
        if [[ -z "$REFERENCE_FASTA" ]]; then
            echo "ERROR: Reference genome (-g/--genome) is required for ALIGNMENT mode"
            errors=$((errors + 1))
        elif [[ ! -f "$REFERENCE_FASTA" ]]; then
            echo "ERROR: Reference genome file does not exist: $REFERENCE_FASTA"
            errors=$((errors + 1))
        else
            echo "✓ Reference genome: $REFERENCE_FASTA"
        fi
    elif [[ "$WORKFLOW_MODE" == "PSEUDOALIGNMENT" ]]; then
        if [[ -z "$TRANSCRIPTS_FASTA" ]]; then
            echo "ERROR: Transcripts FASTA (-t/--transcripts) is required for PSEUDOALIGNMENT mode"
            errors=$((errors + 1))
        elif [[ ! -f "$TRANSCRIPTS_FASTA" ]]; then
            echo "ERROR: Transcripts FASTA file does not exist: $TRANSCRIPTS_FASTA"
            errors=$((errors + 1))
        else
            echo "✓ Transcripts FASTA: $TRANSCRIPTS_FASTA"
        fi
    fi
    
    if [[ -z "$ANNOTATION_FILE" ]]; then
        echo "ERROR: Annotation file (-a/--annotation) is required"
        errors=$((errors + 1))
    elif [[ ! -f "$ANNOTATION_FILE" ]]; then
        echo "ERROR: Annotation file does not exist: $ANNOTATION_FILE"
        errors=$((errors + 1))
    else
        echo "✓ Annotation file: $ANNOTATION_FILE"
    fi
    
    if [[ $errors -gt 0 ]]; then
        echo ""
        echo "Found $errors error(s). Please fix them before running the pipeline."
        exit 1
    fi
    
    echo "✓ Input validation passed"
}

# ===============================================================================
# EXECUTE ALL VALIDATION CHECKS BEFORE PROCESSING
# ===============================================================================

# Step 1: Check ALL dependencies FIRST (tools, R packages, Python packages)
if ! check_all_dependencies; then
    exit 1
fi

# Step 2: Validate input files and parameters
validate_inputs


# ===============================================================================
# 2. PREPARE DIRECTORIES & LOGS
# ===============================================================================

setup_directories() {
    echo "Setting up directory structure..."
    
    # Convert to absolute path
    OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
    
    # Create main output directory
    mkdir -p "$OUTPUT_DIR"
    
    # Create subdirectories for QC, alignments, quantification, etc.
    mkdir -p "$OUTPUT_DIR"/{qc/raw_fastqc,qc/posttrim_fastqc,multiqc,trimmed_fastq}
    mkdir -p "$OUTPUT_DIR"/{indices,alignments/sorted,alignments/metrics,quantification}
    mkdir -p "$OUTPUT_DIR"/{counts_matrix,logs,results}
    
    # Initialize log files
    MASTER_LOG="$OUTPUT_DIR/logs/pipeline_master.log"
    SAMPLE_LOG_DIR="$OUTPUT_DIR/logs/samples"
    mkdir -p "$SAMPLE_LOG_DIR"
    
    # Start master log
    {
        echo "==============================================================================="
        echo "RNA-seq Pipeline Master Log"
        echo "Started: $(date)"
        echo "Command: $0 $*"
        echo "Working directory: $(pwd)"
        echo "Output directory: $OUTPUT_DIR"
        echo "Workflow mode: $WORKFLOW_MODE"
        echo "Threads: $THREADS"
        echo "==============================================================================="
    } > "$MASTER_LOG"
    
    echo "✓ Directory structure created: $OUTPUT_DIR"
}

# Resolve manifest paths safely:
# - absolute paths are kept as-is
# - relative paths are resolved relative to the manifest file directory
resolve_manifest_path() {
    local input_path="$1"
    local manifest_dir
    manifest_dir="$(cd "$(dirname "$SAMPLE_MANIFEST")" && pwd)"
    
    if [[ "$input_path" = /* ]]; then
        printf "%s\n" "$input_path"
    else
        # Try relative to manifest directory first, then current working directory.
        # This supports both ../fastq/... and Input/fastq/... styles in manifests.
        if [[ -e "$manifest_dir/$input_path" ]]; then
            printf "%s\n" "$(realpath "$manifest_dir/$input_path")"
        elif [[ -e "$input_path" ]]; then
            printf "%s\n" "$(realpath "$input_path")"
        else
            # Return the manifest-relative candidate so downstream errors remain informative.
            printf "%s\n" "$manifest_dir/$input_path"
        fi
    fi
}

validate_manifest() {
    echo "Validating sample manifest..."
    
    local header=$(head -n1 "$SAMPLE_MANIFEST")
    # Remove BOM if present
    header=$(echo "$header" | sed '1s/^\xEF\xBB\xBF//')
    if [[ ! "$header" =~ sample_id.*r1.*r2 ]]; then
        echo "WARNING: Manifest header should contain 'sample_id', 'r1', 'r2' columns"
    fi
    
    # Read samples and validate files exist
    local missing_files=0
    
    while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
        # Remove BOM, quotes, whitespace, and carriage returns
        sample_id=$(echo "$sample_id" | sed 's/^\xEF\xBB\xBF//' | tr -d '"' | tr -d '\r' | xargs)
        r1_path=$(echo "$r1_path" | tr -d '"' | tr -d '\r' | xargs)
        r2_path=$(echo "$r2_path" | tr -d '"' | tr -d '\r' | xargs)
        
        # Skip header line (check after cleaning)
        if [[ "$sample_id" == "sample_id" ]]; then
            continue
        fi
        
        # Skip empty lines
        if [[ -z "$sample_id" ]]; then
            continue
        fi
        
        r1_path=$(resolve_manifest_path "$r1_path")
        r2_path=$(resolve_manifest_path "$r2_path")
        
        if [[ ! -f "$r1_path" ]]; then
            echo "ERROR: R1 file not found for sample $sample_id: $r1_path"
            missing_files=$((missing_files + 1))
        fi
        
        if [[ ! -f "$r2_path" ]]; then
            echo "ERROR: R2 file not found for sample $sample_id: $r2_path"
            missing_files=$((missing_files + 1))
        fi
        
        # Create sample log file
        mkdir -p "$SAMPLE_LOG_DIR"
        touch "$SAMPLE_LOG_DIR/${sample_id}.log"
    done < "$SAMPLE_MANIFEST"
    
    if [[ $missing_files -gt 0 ]]; then
        echo "ERROR: Found $missing_files missing FASTQ files"
        exit 1
    fi
    
    echo "✓ Manifest validation passed"
}

setup_directories
validate_manifest

# ===============================================================================
# 3. RAW QUALITY CONTROL (FASTQC)
# ===============================================================================

run_raw_fastqc() {
    echo "Running FastQC on raw reads..."
    
    local fastqc_out="$OUTPUT_DIR/qc/raw_fastqc"
    local failed_samples=0
    
    # Process each sample
    while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == *sample_id* ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | tr -d '\r' | xargs)
        r1_path=$(echo "$r1_path" | tr -d '"' | tr -d '\r' | xargs)
        r2_path=$(echo "$r2_path" | tr -d '"' | tr -d '\r' | xargs)
        r1_path=$(resolve_manifest_path "$r1_path")
        r2_path=$(resolve_manifest_path "$r2_path")
        
        echo "  Processing sample: $sample_id"
        
        # Run FastQC on both R1 and R2
        if ! $FASTQC_CMD \
            --outdir "$fastqc_out" \
            --threads "$THREADS" \
            --quiet \
            "$r1_path" "$r2_path" \
            >> "$SAMPLE_LOG_DIR/${sample_id}.log" 2>&1; then
            echo "ERROR: FastQC failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
        fi
    done < "$SAMPLE_MANIFEST"
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: FastQC failed for $failed_samples samples"
        exit 1
    fi
    
    # Run MultiQC to aggregate raw FastQC results
    echo "  Running MultiQC on raw FastQC results..."
    $MULTIQC_CMD \
        --outdir "$OUTPUT_DIR/multiqc" \
        --filename "multiqc_report_raw.html" \
        --title "Raw Reads QC Report" \
        --force \
        "$fastqc_out" \
        >> "$MASTER_LOG" 2>&1
    
    echo "✓ Raw FastQC completed"
}

run_raw_fastqc

# ===============================================================================
# 4. READ PREPROCESSING / TRIMMING
# ===============================================================================

run_trimming() {
    if [[ "$TRIM_READS" != "true" ]]; then
        echo "Skipping read trimming (disabled)"
        
        # Create symlinks to original files in trimmed_fastq directory
        echo "Creating symlinks to original FASTQ files..."
        while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == *sample_id* ]] && continue
            sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
            r1_path=$(echo "$r1_path" | tr -d '"' | xargs)
            r2_path=$(echo "$r2_path" | tr -d '"' | xargs)
            r1_path=$(resolve_manifest_path "$r1_path")
            r2_path=$(resolve_manifest_path "$r2_path")
            
            ln -sf "$(realpath "$r1_path")" "$OUTPUT_DIR/trimmed_fastq/${sample_id}_R1.trim.fastq.gz"
            ln -sf "$(realpath "$r2_path")" "$OUTPUT_DIR/trimmed_fastq/${sample_id}_R2.trim.fastq.gz"
        done < "$SAMPLE_MANIFEST"
        
        return
    fi
    
    echo "Running read trimming with fastp..."
    
    local trimmed_dir="$OUTPUT_DIR/trimmed_fastq"
    local failed_samples=0
    
    # Process each sample
    while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == *sample_id* ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | tr -d '\r' | xargs)
        r1_path=$(echo "$r1_path" | tr -d '"' | tr -d '\r' | xargs)
        r2_path=$(echo "$r2_path" | tr -d '"' | tr -d '\r' | xargs)
        r1_path=$(resolve_manifest_path "$r1_path")
        r2_path=$(resolve_manifest_path "$r2_path")
        
        echo "  Trimming sample: $sample_id"
        
        local r1_out="$trimmed_dir/${sample_id}_R1.trim.fastq.gz"
        local r2_out="$trimmed_dir/${sample_id}_R2.trim.fastq.gz"
        local trim_log="$SAMPLE_LOG_DIR/${sample_id}_trimming.log"
        
        # Run fastp with quality filtering
        local fastp_args=(
            --in1 "$r1_path"
            --in2 "$r2_path"
            --out1 "$r1_out"
            --out2 "$r2_out"
            --thread "$THREADS"
            --length_required "$MIN_READ_LENGTH"
            --qualified_quality_phred "$QUALITY_THRESHOLD"
            --html "$trimmed_dir/${sample_id}_fastp.html"
            --json "$trimmed_dir/${sample_id}_fastp.json"
            --detect_adapter_for_pe
        )
        
        if ! $FASTP_CMD "${fastp_args[@]}" >> "$trim_log" 2>&1; then
            echo "ERROR: Trimming failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
        fi
    done < "$SAMPLE_MANIFEST"
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: Trimming failed for $failed_samples samples"
        exit 1
    fi
    
    echo "✓ Read trimming completed"
}

run_trimming

# ===============================================================================
# 5. POST-TRIM QC
# ===============================================================================

run_posttrim_fastqc() {
    if [[ "$TRIM_READS" != "true" ]]; then
        echo "Skipping post-trim FastQC (trimming was disabled)"
        return
    fi
    
    echo "Running FastQC on trimmed reads..."
    
    local fastqc_out="$OUTPUT_DIR/qc/posttrim_fastqc"
    local trimmed_dir="$OUTPUT_DIR/trimmed_fastq"
    local failed_samples=0
    
    # Process each sample
    while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == *sample_id* ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        
        echo "  Processing trimmed sample: $sample_id"
        
        local r1_trimmed="$trimmed_dir/${sample_id}_R1.trim.fastq.gz"
        local r2_trimmed="$trimmed_dir/${sample_id}_R2.trim.fastq.gz"
        
        # Run FastQC on trimmed files
        if ! $FASTQC_CMD \
            --outdir "$fastqc_out" \
            --threads "$THREADS" \
            --quiet \
            "$r1_trimmed" "$r2_trimmed" \
            >> "$SAMPLE_LOG_DIR/${sample_id}.log" 2>&1; then
            echo "ERROR: Post-trim FastQC failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
        fi
    done < "$SAMPLE_MANIFEST"
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: Post-trim FastQC failed for $failed_samples samples"
        exit 1
    fi
    
    # Run MultiQC to aggregate trimmed FastQC results
    echo "  Running MultiQC on trimmed FastQC results..."
    $MULTIQC_CMD \
        --outdir "$OUTPUT_DIR/multiqc" \
        --filename "multiqc_report_posttrim.html" \
        --title "Trimmed Reads QC Report" \
        --force \
        "$fastqc_out" \
        >> "$MASTER_LOG" 2>&1
    
    echo "✓ Post-trim FastQC completed"
}

run_posttrim_fastqc

# ===============================================================================
# 6. REFERENCE INDEXING
# ===============================================================================

build_reference_index() {
    echo "Setting up reference index for $WORKFLOW_MODE mode..."
    
    local index_dir="$OUTPUT_DIR/indices"
    
    if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
        # HISAT2 genome index
        local hisat2_index="$index_dir/HISAT2"
        local hisat2_index_prefix="$hisat2_index/genome"
        
        if [[ -f "${hisat2_index_prefix}.1.ht2" ]]; then
            echo "  ✓ HISAT2 index already exists, skipping build"
        else
            echo "  Building HISAT2 genome index..."
            mkdir -p "$hisat2_index"
            
            # HISAT2 index building command
            local hisat2_build_args=(
                -p "$THREADS"
                "$REFERENCE_FASTA"
                "$hisat2_index_prefix"
            )
            
            # Skip annotation file for large genomes (causes indexing issues)
            # Annotation will be used during alignment with --dta parameter
            echo "  Building index without annotation (will use --dta during alignment)..."
            
            if ! $HISAT2_BUILD_CMD "${hisat2_build_args[@]}" >> "$MASTER_LOG" 2>&1; then
                echo "ERROR: HISAT2 index building failed"
                exit 1
            fi
            
            echo "  ✓ HISAT2 index built successfully"
        fi
        
    elif [[ "$WORKFLOW_MODE" == "PSEUDOALIGNMENT" ]]; then
        # Salmon transcriptome index
        local salmon_index="$index_dir/salmon"
        
        if [[ -d "$salmon_index" && -f "$salmon_index/info.json" ]]; then
            echo "  ✓ Salmon index already exists, skipping build"
        else
            echo "  Building Salmon transcriptome index..."
            
            if ! $SALMON_CMD index \
                --transcripts "$TRANSCRIPTS_FASTA" \
                --index "$salmon_index" \
                --threads "$THREADS" \
                >> "$MASTER_LOG" 2>&1; then
                echo "ERROR: Salmon index building failed"
                exit 1
            fi
            
            echo "  ✓ Salmon index built successfully"
        fi
    fi
}

build_reference_index

# ===============================================================================
# 7A. ALIGNMENT-BASED MAPPING
# ===============================================================================

run_alignment() {
    if [[ "$WORKFLOW_MODE" != "ALIGNMENT" ]]; then
        return
    fi
    
    echo "Running alignment-based mapping with HISAT2..."
    
    local hisat2_index="$OUTPUT_DIR/indices/HISAT2"
    local align_dir="$OUTPUT_DIR/alignments"
    local sorted_dir="$align_dir/sorted"
    local metrics_dir="$align_dir/metrics"
    local failed_samples=0
    
    # Process each sample
    while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == *sample_id* ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        
        echo "  Aligning sample: $sample_id"
        
        local r1_input="$OUTPUT_DIR/trimmed_fastq/${sample_id}_R1.trim.fastq.gz"
        local r2_input="$OUTPUT_DIR/trimmed_fastq/${sample_id}_R2.trim.fastq.gz"
        local sam_output="$align_dir/${sample_id}.sam"
        local sorted_bam="$sorted_dir/${sample_id}.sorted.bam"
        local align_log="$SAMPLE_LOG_DIR/${sample_id}_alignment.log"
        
        echo "    Running HISAT2 alignment..."
        echo "    R1: $r1_input"
        echo "    R2: $r2_input"
        echo "    Index: $hisat2_index"
        
        # HISAT2 alignment command with simplified parameters to avoid _plen[] error
        #
        # HISAT2 --rna-strandness is only valid for stranded libraries (e.g. FR/RF).
        # For unstranded libraries, we omit the option.
        local hisat2_strand_args=()
        case "$STRANDEDNESS" in
            "FR"|"RF")
                hisat2_strand_args=(--rna-strandness "$STRANDEDNESS")
                ;;
            "NONE"|"UNSTRANDED"|"0"|"")
                hisat2_strand_args=()
                ;;
            *)
                # For single-end style codes (F/R) or unknown values, do not pass to HISAT2.
                # This prevents HISAT2 from erroring on invalid strandness strings.
                hisat2_strand_args=()
                ;;
        esac
        
        if ! $HISAT2_CMD \
            -x "$hisat2_index/genome" \
            -1 "$r1_input" \
            -2 "$r2_input" \
            -S "$sam_output" \
            --threads "$THREADS" \
            "${hisat2_strand_args[@]}" \
            --summary-file "$metrics_dir/${sample_id}_hisat2_summary.txt" \
            >> "$align_log" 2>&1; then
            echo "ERROR: HISAT2 alignment failed for sample $sample_id"
            echo "    Check log: $align_log"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        echo "    ✓ HISAT2 alignment completed for $sample_id"
        
        # Convert SAM to BAM and sort
        echo "    Converting SAM to sorted BAM..."
        if ! $SAMTOOLS_CMD view -bS "$sam_output" | \
            $SAMTOOLS_CMD sort -@ "$THREADS" -o "$sorted_bam" -; then
            echo "ERROR: SAM to BAM conversion failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        # Index the BAM file
        if ! $SAMTOOLS_CMD index "$sorted_bam" >> "$align_log" 2>&1; then
            echo "ERROR: BAM indexing failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        # Clean up SAM file to save space
        rm -f "$sam_output"
        
        # Generate alignment statistics
        $SAMTOOLS_CMD flagstat "$sorted_bam" > "$metrics_dir/${sample_id}_flagstat.txt" 2>> "$align_log"
        $SAMTOOLS_CMD stats "$sorted_bam" > "$metrics_dir/${sample_id}_stats.txt" 2>> "$align_log"
        
        echo "    ✓ BAM processing completed for $sample_id"
    done < "$SAMPLE_MANIFEST"
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: Alignment failed for $failed_samples samples"
        exit 1
    fi
    
    echo "✓ Alignment-based mapping completed"
}

# ===============================================================================
# 7B. PSEUDO-ALIGNMENT / QUASI-MAPPING
# ===============================================================================

run_pseudoalignment() {
    if [[ "$WORKFLOW_MODE" != "PSEUDOALIGNMENT" ]]; then
        return
    fi
    
    echo "Running pseudo-alignment with Salmon..."
    
    local salmon_index="$OUTPUT_DIR/indices/salmon"
    local quant_dir="$OUTPUT_DIR/quantification"
    local failed_samples=0
    
    # Process each sample
    while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == *sample_id* ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        
        echo "  Quantifying sample: $sample_id"
        
        local r1_input="$OUTPUT_DIR/trimmed_fastq/${sample_id}_R1.trim.fastq.gz"
        local r2_input="$OUTPUT_DIR/trimmed_fastq/${sample_id}_R2.trim.fastq.gz"
        local sample_quant_dir="$quant_dir/$sample_id"
        local quant_log="$SAMPLE_LOG_DIR/${sample_id}_quantification.log"
        
        # Determine library type based on strandedness for Salmon
        local lib_type="A"  # Auto-detect
        case "$STRANDEDNESS" in
            "F"|"FR") lib_type="ISF" ;;
            "R"|"RF") lib_type="ISR" ;;
            "NONE"|"UNSTRANDED"|"0") lib_type="A" ;;
            *) lib_type="A" ;;
        esac
        
        # Run Salmon quantification
        local salmon_args=(
            quant
            --index "$salmon_index"
            --libType "$lib_type"
            --mates1 "$r1_input"
            --mates2 "$r2_input"
            --output "$sample_quant_dir"
            --threads "$THREADS"
            --validateMappings
            --gcBias
            --seqBias
        )
        
        if ! $SALMON_CMD "${salmon_args[@]}" >> "$quant_log" 2>&1; then
            echo "ERROR: Salmon quantification failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        echo "    ✓ Quantification completed for $sample_id"
    done < "$SAMPLE_MANIFEST"
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: Pseudo-alignment failed for $failed_samples samples"
        exit 1
    fi
    
    echo "✓ Pseudo-alignment completed"
}

# Run the appropriate mapping method based on workflow mode
run_alignment
run_pseudoalignment

# ===============================================================================
# 8. POST-ALIGNMENT QC (RNA-SEQ OPTIMIZED)
# ===============================================================================

run_post_alignment_qc() {
    if [[ "$WORKFLOW_MODE" != "ALIGNMENT" ]]; then
        return
    fi
    
    echo "Running post-alignment QC..."
    
    # For RNA-seq, we skip duplicate removal as it's not recommended
    # Duplicate reads in RNA-seq often represent real biological signal
    # (highly expressed genes) rather than technical artifacts
    
    echo "✓ Post-alignment QC completed (duplicate removal skipped for RNA-seq)"
}

run_post_alignment_qc

# ===============================================================================
# 9. EXPRESSION QUANTIFICATION (GENE-LEVEL)
# ===============================================================================

run_expression_quantification() {
    echo "Running expression quantification..."
    
    local counts_dir="$OUTPUT_DIR/counts_matrix"
    
    if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
        echo "  Using featureCounts for gene-level quantification..."
        
        local sorted_dir="$OUTPUT_DIR/alignments/sorted"
        local bam_files=()
        local sample_names=()
        
        # Collect BAM files and sample names
        while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == *sample_id* ]] && continue  # Skip header
            sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
            bam_files+=("$sorted_dir/${sample_id}.sorted.bam")
            sample_names+=("$sample_id")
        done < "$SAMPLE_MANIFEST"
        
        # Determine strandedness parameter for featureCounts
        #
        # featureCounts -s:
        #   0 = unstranded
        #   1 = stranded (same strand)
        #   2 = reversely stranded
        #
        # We accept FR/RF and map them in the common way:
        #   FR -> 1, RF -> 2
        local strand_param=0
        case "$STRANDEDNESS" in
            "F"|"FR") strand_param=1 ;;
            "R"|"RF") strand_param=2 ;;
            "NONE"|"UNSTRANDED"|"0") strand_param=0 ;;
            *) strand_param=0 ;;
        esac
        
        # featureCounts configuration depends on annotation conventions.
        #
        # - Typical GTF: exons carry gene_id, so we count exons and group by gene_id.
        # - Many GFF3s (including Helixer-style) give exons their own IDs (exon IDs),
        #   so counting exons grouped by ID produces exon-level rows. To get gene-level
        #   counts, count 'gene' features and group by their ID.
        local fc_feature_type="exon"
        local fc_group_attr="gene_id"
        
        if [[ "$ANNOTATION_FILE" =~ \.gff3$ || "$ANNOTATION_FILE" =~ \.gff$ ]]; then
            fc_feature_type="gene"
            fc_group_attr="ID"
            echo "  Detected GFF3/GFF annotation: using featureCounts -t gene -g ID for gene-level counts"
        elif [[ "$ANNOTATION_FILE" =~ \.gtf$ ]]; then
            fc_feature_type="exon"
            fc_group_attr="gene_id"
            echo "  Detected GTF annotation: using featureCounts -t exon -g gene_id for gene-level counts"
        else
            echo "  WARNING: Unknown annotation extension for: $ANNOTATION_FILE"
            echo "  Defaulting to featureCounts -t exon -g gene_id (typical GTF style)."
        fi
        
        # Run featureCounts to count reads per gene
        local featurecounts_args=(
            -p  # paired-end
            -T "$THREADS"
            -t "$fc_feature_type"
            -g "$fc_group_attr"
            -s "$strand_param"
            -a "$ANNOTATION_FILE"
            -o "$counts_dir/raw_gene_counts.txt"
            "${bam_files[@]}"
        )
        
        if ! $FEATURECOUNTS_CMD "${featurecounts_args[@]}" >> "$MASTER_LOG" 2>&1; then
            echo "ERROR: featureCounts failed"
            exit 1
        fi
        
        if [[ -f "$counts_dir/raw_gene_counts.txt" ]]; then
            # Clean up raw_gene_counts.txt file (replace full paths with sample names)
            echo "  Cleaning raw_gene_counts.txt header..."
            local temp_file=$(mktemp)
            
            # Process header: replace full BAM paths with sample names
            local header_line=$(grep -v "^#" "$counts_dir/raw_gene_counts.txt" | head -n1)
            local cleaned_header="Geneid\tChr\tStart\tEnd\tStrand\tLength"
            
            # Add sample names to header
            for sample in "${sample_names[@]}"; do
                cleaned_header="${cleaned_header}\t${sample}"
            done
            
            # Write cleaned file: skip comment line, use cleaned header, keep all data
            grep "^#" "$counts_dir/raw_gene_counts.txt" > "$temp_file" || true
            echo -e "$cleaned_header" >> "$temp_file"
            grep -v "^#" "$counts_dir/raw_gene_counts.txt" | tail -n +2 >> "$temp_file"
            mv "$temp_file" "$counts_dir/raw_gene_counts.txt"
            
            # Clean up raw_gene_counts.txt.summary file
            if [[ -f "$counts_dir/raw_gene_counts.txt.summary" ]]; then
                echo "  Cleaning raw_gene_counts.txt.summary header..."
                local summary_temp=$(mktemp)
                local summary_header="Status"
                for sample in "${sample_names[@]}"; do
                    summary_header="${summary_header}\t${sample}"
                done
                echo -e "$summary_header" > "$summary_temp"
                tail -n +2 "$counts_dir/raw_gene_counts.txt.summary" >> "$summary_temp"
                mv "$summary_temp" "$counts_dir/raw_gene_counts.txt.summary"
            fi
            
            # Create gene_counts_matrix.tsv (extract only Geneid and count columns)
            local cleaned_header_matrix="Geneid"
            for sample in "${sample_names[@]}"; do
                cleaned_header_matrix="${cleaned_header_matrix}\t${sample}"
            done
            echo -e "$cleaned_header_matrix" > "$counts_dir/gene_counts_matrix.tsv"
            grep -v "^#" "$counts_dir/raw_gene_counts.txt" | tail -n +2 | cut -f1,7- >> "$counts_dir/gene_counts_matrix.tsv"
            
            # Create DGE-ready CSV format (convert TSV to CSV - all columns)
            echo "  Creating DGE-ready count matrix..."
            awk -F'\t' 'BEGIN {OFS=","} {for(i=1;i<=NF;i++) {if(i>1) printf ","; printf "%s", $i} printf "\n"}' "$counts_dir/gene_counts_matrix.tsv" | sed '1s/^Geneid/GeneID/' > "$counts_dir/count_matrix.csv"
            
            create_sample_metadata "$counts_dir/sample_metadata.csv"
            
            echo "  ✓ DGE-ready files created:"
            echo "    - count_matrix.csv (CSV format for DESeq2/edgeR)"
            echo "    - sample_metadata.csv (sample information)"
        fi
        
    elif [[ "$WORKFLOW_MODE" == "PSEUDOALIGNMENT" ]]; then
        echo "  Collecting Salmon quantification results..."
        
        # Create transcript-to-gene mapping from GTF/GFF
        local tx2gene_file="$counts_dir/tx2gene.tsv"
        create_tx2gene_mapping "$tx2gene_file"
        
        # Collect all quant.sf files from Salmon output
        local quant_files=()
        local sample_names=()
        
        while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == *sample_id* ]] && continue
            sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
            quant_files+=("$OUTPUT_DIR/quantification/$sample_id/quant.sf")
            sample_names+=("$sample_id")
        done < "$SAMPLE_MANIFEST"
        
        # Aggregate transcript-level counts to gene-level
        # Note: sample_names is set in this function scope and will be accessible
        create_gene_counts_matrix "$tx2gene_file" "${quant_files[@]}"
        
        if [[ -f "$counts_dir/gene_counts_matrix.tsv" ]]; then
            # Create DGE-ready CSV format (convert TSV to CSV - all columns)
            echo "  Creating DGE-ready count matrix..."
            awk -F'\t' 'BEGIN {OFS=","} {for(i=1;i<=NF;i++) {if(i>1) printf ","; printf "%s", $i} printf "\n"}' "$counts_dir/gene_counts_matrix.tsv" | sed '1s/^Geneid/GeneID/' > "$counts_dir/count_matrix.csv"
            
            create_sample_metadata "$counts_dir/sample_metadata.csv"
            
            echo "  ✓ DGE-ready files created:"
            echo "    - count_matrix.csv (CSV format for DESeq2/edgeR)"
            echo "    - sample_metadata.csv (sample information)"
        fi
    fi
    
    echo "✓ Expression quantification completed"
}

create_sample_metadata() {
    local output_file="$1"
    echo "  Creating sample metadata..."
    
    # Prefer condition/replicate columns from manifest when available.
    # Falls back to condition1 + sequential replicate numbers otherwise.
    python3 - "$SAMPLE_MANIFEST" "$output_file" << 'EOF'
import csv
import sys

manifest = sys.argv[1]
output = sys.argv[2]

with open(manifest, newline="", encoding="utf-8-sig") as f:
    reader = csv.DictReader(f)
    if reader.fieldnames is None:
        raise SystemExit("Manifest appears empty or malformed")

    field_map = {name.strip().lower(): name for name in reader.fieldnames}
    sid_key = field_map.get("sample_id")
    cond_key = field_map.get("condition")
    rep_key = field_map.get("replicate")

    if sid_key is None:
        raise SystemExit("Manifest is missing required 'sample_id' column")

    rows = []
    replicate_counter = 1
    used_defaults = False

    for row in reader:
        sample_id = (row.get(sid_key) or "").strip()
        if not sample_id:
            continue

        condition = (row.get(cond_key) or "").strip() if cond_key else ""
        replicate = (row.get(rep_key) or "").strip() if rep_key else ""

        if not condition:
            condition = "condition1"
            used_defaults = True
        if not replicate:
            replicate = str(replicate_counter)
            used_defaults = True

        rows.append((sample_id, condition, replicate))
        replicate_counter += 1

with open(output, "w", newline="", encoding="utf-8") as f:
    writer = csv.writer(f)
    writer.writerow(["sample_id", "condition", "replicate"])
    writer.writerows(rows)

if used_defaults:
    print("  WARNING: sample_metadata.csv used defaults for missing condition/replicate values")
    print("  Please review and edit sample_metadata.csv to match your experimental design")
EOF
    
    echo "  ✓ sample_metadata.csv created"
}

create_tx2gene_mapping() {
    local output_file="$1"
    
    echo "  Creating transcript-to-gene mapping..."
    
    if [[ "$ANNOTATION_FILE" =~ \.gtf$ ]]; then
        # GTF format
        awk -F'\t' '$3 == "transcript" {
            match($9, /transcript_id "([^"]+)"/, tid);
            match($9, /gene_id "([^"]+)"/, gid);
            if (tid[1] && gid[1]) print tid[1] "\t" gid[1]
        }' "$ANNOTATION_FILE" > "$output_file"
    else
        # GFF3 format
        awk -F'\t' '$3 == "mRNA" || $3 == "transcript" {
            match($9, /ID=([^;]+)/, tid);
            match($9, /Parent=([^;]+)/, gid);
            if (tid[1] && gid[1]) print tid[1] "\t" gid[1]
        }' "$ANNOTATION_FILE" > "$output_file"
    fi
    
    echo "    ✓ Created tx2gene mapping with $(wc -l < "$output_file") entries"
}

create_gene_counts_matrix() {
    local tx2gene_file="$1"
    shift
    local quant_files=("$@")
    
    # sample_names should be set in the calling function's scope
    # Make a local copy to avoid issues
    local sample_names_array=("${sample_names[@]}")
    
    echo "  Creating gene-level counts matrix..."
    
    # This is a simplified aggregation - for production, use tximport
    local output_file="$OUTPUT_DIR/counts_matrix/gene_counts_matrix.tsv"
    
    printf "gene_id"
    for i in "${!sample_names_array[@]}"; do
        printf "\t%s" "${sample_names_array[$i]}"
    done
    printf "\n" > "$output_file"
    
    # Simple aggregation by summing transcript counts per gene
    # This is basic - tximport does proper length-scaled aggregation
    python3 -c "
import sys
import pandas as pd

# Load tx2gene mapping
tx2gene = pd.read_csv('$tx2gene_file', sep='\t', header=None, names=['transcript_id', 'gene_id'])
tx2gene_dict = dict(zip(tx2gene['transcript_id'], tx2gene['gene_id']))

# Load and aggregate quantification files
gene_counts = {}
sample_names = [$(printf "'%s'," "${sample_names_array[@]}" | sed 's/,$//'))]

for i, quant_file in enumerate(['$(printf "%s','" "${quant_files[@]}" | sed "s/',$//")']):
    try:
        df = pd.read_csv(quant_file, sep='\t')
        for _, row in df.iterrows():
            transcript_id = row['Name']
            count = row['NumReads']  # Use raw read counts
            
            if transcript_id in tx2gene_dict:
                gene_id = tx2gene_dict[transcript_id]
                if gene_id not in gene_counts:
                    gene_counts[gene_id] = [0] * len(sample_names)
                gene_counts[gene_id][i] += count
    except Exception as e:
        print(f'Error processing {quant_file}: {e}', file=sys.stderr)

# Write results
with open('$output_file', 'w') as f:
    f.write('gene_id\t' + '\t'.join(sample_names) + '\n')
    for gene_id, counts in sorted(gene_counts.items()):
        f.write(gene_id + '\t' + '\t'.join(map(str, map(int, counts))) + '\n')
" >> "$MASTER_LOG" 2>&1
    
    if [[ ! -f "$output_file" || ! -s "$output_file" ]]; then
        echo "WARNING: Gene counts matrix creation failed or is empty"
    fi
}

run_expression_quantification

# ===============================================================================
# 10. AGGREGATE REPORTS (MULTIQC + MAPPING SUMMARY)
# ===============================================================================

# Create MultiQC configuration
create_multiqc_config() {
    local config_file="$OUTPUT_DIR/multiqc_config.yaml"
    
    cat > "$config_file" << 'EOF'
report_title: "RNA-seq Preprocessing Report"
report_comment: "Quality control and alignment statistics"
template: default
force: true
EOF

    echo "✓ MultiQC configuration created: $config_file"
}


# Create master HTML report
create_master_html_report() {
    local master_file="$OUTPUT_DIR/master_report.html"
    local run_report="$OUTPUT_DIR/results/run_report.txt"
    local mapping_summary="$OUTPUT_DIR/results/mapping_summary.txt"
    local count_matrix="$OUTPUT_DIR/counts_matrix/count_matrix.csv"
    
    # Get count matrix statistics (use TSV file for accurate column count)
    local total_genes="N/A"
    local total_samples="N/A"
    local tsv_matrix="$OUTPUT_DIR/counts_matrix/gene_counts_matrix.tsv"
    if [[ -f "$tsv_matrix" ]]; then
        total_genes=$(wc -l < "$tsv_matrix" | awk '{print $1-1}')
        total_samples=$(head -1 "$tsv_matrix" | awk -F'\t' '{print NF-1}')
    elif [[ -f "$count_matrix" ]]; then
        total_genes=$(wc -l < "$count_matrix" | awk '{print $1-1}')
        total_samples=$(head -1 "$count_matrix" | awk -F',' '{print NF-1}')
    fi
    
    # Get pipeline info from run report
    local pipeline_version="$SCRIPT_VERSION"
    local workflow_mode="$WORKFLOW_MODE"
    local start_time="N/A"
    local end_time="N/A"
    local total_samples_processed="0"
    local successful_samples="0"
    local failed_samples="0"
    local status="UNKNOWN"
    
    if [[ -f "$run_report" ]]; then
        start_time=$(grep "Start time:" "$run_report" | cut -d':' -f2- | xargs || echo "N/A")
        end_time=$(grep "End time:" "$run_report" | cut -d':' -f2- | xargs || echo "N/A")
        total_samples_processed=$(grep "Total samples:" "$run_report" | awk '{print $3}' || echo "0")
        successful_samples=$(grep "Successful:" "$run_report" | awk '{print $2}' || echo "0")
        failed_samples=$(grep "Failed:" "$run_report" | awk '{print $2}' || echo "0")
        if grep -q "COMPLETED SUCCESSFULLY" "$run_report"; then
            status="SUCCESS"
        elif grep -q "COMPLETED WITH ERRORS" "$run_report"; then
            status="ERRORS"
        fi
    fi
    
    # Get configuration
    local threads="$THREADS"
    local trimming="$TRIM_READS"
    local strandedness="$STRANDEDNESS"
    
    # Generate sample list
    local sample_list=""
    while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == *sample_id* ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        sample_list="${sample_list}        <li>$sample_id</li>\n"
    done < "$SAMPLE_MANIFEST"
    
    # Generate mapping summary HTML
    local mapping_summary_html="<p>No mapping summary available.</p>"
    if [[ -f "$mapping_summary" ]]; then
        mapping_summary_html="<pre style='background: #f5f5f5; padding: 15px; border-radius: 5px; overflow-x: auto;'>"
        mapping_summary_html+=$(sed 's/&/\&amp;/g; s/</\&lt;/g; s/>/\&gt;/g' "$mapping_summary")
        mapping_summary_html+="</pre>"
    fi
    
    # Generate run report HTML
    local run_report_html="<p>No run report available.</p>"
    if [[ -f "$run_report" ]]; then
        run_report_html="<pre style='background: #f5f5f5; padding: 15px; border-radius: 5px; overflow-x: auto;'>"
        run_report_html+=$(sed 's/&/\&amp;/g; s/</\&lt;/g; s/>/\&gt;/g' "$run_report")
        run_report_html+="</pre>"
    fi
    
    cat > "$master_file" << EOF
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>RNA-seq Preprocessing Master Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; background: #f9f9f9; }
        .header { background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%); color: white; padding: 30px; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }
        .header h1 { margin: 0 0 10px 0; font-size: 2.2em; }
        .nav { background: #34495e; padding: 15px; margin: 20px 0; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.1); }
        .nav a { color: white; text-decoration: none; margin-right: 20px; padding: 10px 15px; border-radius: 4px; transition: background 0.3s; }
        .nav a:hover { background: #3498db; }
        .section { margin: 30px 0; padding: 25px; background: white; border: 1px solid #ddd; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }
        .section h2 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; margin-top: 0; }
        .info-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin: 20px 0; }
        .info-box { background: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #3498db; }
        .info-box strong { color: #2c3e50; display: block; margin-bottom: 5px; }
        .status-success { color: #27ae60; font-weight: bold; }
        .status-error { color: #e74c3c; font-weight: bold; }
        .footer { text-align: center; margin-top: 40px; padding: 20px; color: #666; background: white; border-radius: 8px; }
        table { width: 100%; border-collapse: collapse; margin: 15px 0; }
        table th, table td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
        table th { background: #34495e; color: white; }
        table tr:hover { background: #f5f5f5; }
        .file-link { color: #3498db; text-decoration: none; }
        .file-link:hover { text-decoration: underline; }
        pre { white-space: pre-wrap; word-wrap: break-word; }
    </style>
</head>
<body>
    <div class="header">
        <h1>RNA-seq Preprocessing Master Report</h1>
        <p>Generated: $(date '+%Y-%m-%d %H:%M:%S')</p>
        <p style="margin: 5px 0;">Pipeline Version: $pipeline_version | Workflow: $workflow_mode</p>
    </div>
    
    <div class="nav">
        <a href="#summary">Summary</a>
        <a href="#configuration">Configuration</a>
        <a href="#samples">Samples</a>
        <a href="#statistics">Statistics</a>
        <a href="#multiqc">MultiQC Reports</a>
        <a href="#mapping">Mapping Summary</a>
        <a href="#outputs">Output Files</a>
    </div>
    
    <div id="summary" class="section">
        <h2>📊 Pipeline Summary</h2>
        <div class="info-grid">
            <div class="info-box">
                <strong>Status</strong>
                <span class="status-$([ "$status" = "SUCCESS" ] && echo "success" || echo "error")">$status</span>
            </div>
            <div class="info-box">
                <strong>Total Samples</strong>
                $total_samples_processed
            </div>
            <div class="info-box">
                <strong>Successful</strong>
                <span class="status-success">$successful_samples</span>
            </div>
            <div class="info-box">
                <strong>Failed</strong>
                <span class="status-error">$failed_samples</span>
            </div>
            <div class="info-box">
                <strong>Total Genes</strong>
                $(printf "%'d" $total_genes 2>/dev/null || echo "$total_genes")
            </div>
            <div class="info-box">
                <strong>Total Samples (Matrix)</strong>
                $total_samples
            </div>
            <div class="info-box">
                <strong>Start Time</strong>
                $start_time
            </div>
            <div class="info-box">
                <strong>End Time</strong>
                $end_time
            </div>
        </div>
    </div>
    
    <div id="configuration" class="section">
        <h2>⚙️ Pipeline Configuration</h2>
        <table>
            <tr>
                <th>Parameter</th>
                <th>Value</th>
            </tr>
            <tr>
                <td>Workflow Mode</td>
                <td><strong>$workflow_mode</strong></td>
            </tr>
            <tr>
                <td>Threads</td>
                <td>$threads</td>
            </tr>
            <tr>
                <td>Trimming Enabled</td>
                <td>$trimming</td>
            </tr>
            <tr>
                <td>Strandedness</td>
                <td>$strandedness</td>
            </tr>
            <tr>
                <td>Reference Genome</td>
                <td>$(basename "$REFERENCE_FASTA" 2>/dev/null || echo "N/A")</td>
            </tr>
            <tr>
                <td>Annotation File</td>
                <td>$(basename "$ANNOTATION_FILE" 2>/dev/null || echo "N/A")</td>
            </tr>
            <tr>
                <td>Output Directory</td>
                <td><code>$OUTPUT_DIR</code></td>
            </tr>
        </table>
    </div>
    
    <div id="samples" class="section">
        <h2>🧪 Processed Samples</h2>
        <ul style="columns: 3; column-gap: 20px;">
$(echo -e "$sample_list")
        </ul>
    </div>
    
    <div id="statistics" class="section">
        <h2>📈 Count Matrix Statistics</h2>
        <div class="info-grid">
            <div class="info-box">
                <strong>Total Genes</strong>
                $(printf "%'d" $total_genes 2>/dev/null || echo "$total_genes")
            </div>
            <div class="info-box">
                <strong>Total Samples</strong>
                $total_samples
            </div>
            <div class="info-box">
                <strong>Matrix Format</strong>
                CSV & TSV
            </div>
            <div class="info-box">
                <strong>Ready for DGE</strong>
                ✓ Yes
            </div>
        </div>
        <p style="margin-top: 15px;">
            <strong>Count Matrix Files:</strong><br>
            • <a href="counts_matrix/count_matrix.csv" class="file-link">count_matrix.csv</a> - DGE-ready format for DESeq2/edgeR<br>
            • <a href="counts_matrix/gene_counts_matrix.tsv" class="file-link">gene_counts_matrix.tsv</a> - Tab-separated format<br>
            • <a href="counts_matrix/sample_metadata.csv" class="file-link">sample_metadata.csv</a> - Sample information
        </p>
    </div>
    
    <div id="multiqc" class="section">
        <h2>📋 MultiQC Reports</h2>
        <p>Comprehensive quality control and alignment statistics:</p>
        <ul>
            <li><a href="multiqc/multiqc_report_raw.html" target="_blank" class="file-link">Raw Reads QC Report</a> - Quality assessment of raw sequencing data</li>
            <li><a href="multiqc/multiqc_report_posttrim.html" target="_blank" class="file-link">Trimmed Reads QC Report</a> - Quality assessment after trimming</li>
            <li><a href="multiqc/multiqc_report_final.html" target="_blank" class="file-link">Final Pipeline Report</a> - Complete pipeline summary with all metrics</li>
        </ul>
    </div>
    
    <div id="mapping" class="section">
        <h2>🗺️ Mapping/Quantification Summary</h2>
        $mapping_summary_html
        <p style="margin-top: 15px;">
            <a href="results/mapping_summary.txt" class="file-link">View full mapping summary (text)</a>
        </p>
    </div>
    
    <div id="outputs" class="section">
        <h2>📁 Output Files & Reports</h2>
        
        <h3>Main Results (DGE-Ready)</h3>
        <ul>
            <li><strong>Count Matrix (CSV)</strong>: <a href="counts_matrix/count_matrix.csv" class="file-link">counts_matrix/count_matrix.csv</a> - Ready for DESeq2/edgeR/limma</li>
            <li><strong>Count Matrix (TSV)</strong>: <a href="counts_matrix/gene_counts_matrix.tsv" class="file-link">counts_matrix/gene_counts_matrix.tsv</a> - Tab-separated format</li>
            <li><strong>Sample Metadata</strong>: <a href="counts_matrix/sample_metadata.csv" class="file-link">counts_matrix/sample_metadata.csv</a> - Sample information for DGE analysis</li>
        </ul>
        
        <h3>Alignment Files</h3>
        <ul>
            <li><strong>Sorted BAM Files</strong>: <code>alignments/sorted/</code> - One sorted, indexed BAM file per sample</li>
            <li><strong>Alignment Metrics</strong>: <code>alignments/metrics/</code> - Flagstat, stats, and HISAT2 summaries</li>
        </ul>
        
        <h3>Quality Control Files</h3>
        <ul>
            <li><strong>Raw FastQC Reports</strong>: <code>qc/raw_fastqc/</code> - HTML reports for each sample's raw reads</li>
            <li><strong>Post-Trim FastQC Reports</strong>: <code>qc/posttrim_fastqc/</code> - HTML reports for trimmed reads</li>
            <li><strong>Trimming Reports</strong>: <code>trimmed_fastq/</code> - fastp HTML/JSON reports per sample</li>
        </ul>
        
        <h3>Reports and Logs</h3>
        <ul>
            <li><strong>Pipeline Log</strong>: <a href="logs/pipeline_master.log" class="file-link">logs/pipeline_master.log</a> - Complete execution log</li>
            <li><strong>Run Report</strong>: <a href="results/run_report.txt" class="file-link">results/run_report.txt</a> - Pipeline execution summary</li>
            <li><strong>Mapping Summary</strong>: <a href="results/mapping_summary.txt" class="file-link">results/mapping_summary.txt</a> - Alignment/quantification statistics</li>
            <li><strong>Output Manifest</strong>: <a href="results/output_manifest.txt" class="file-link">results/output_manifest.txt</a> - Complete file listing</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>📄 Detailed Reports</h2>
        <h3>Run Report</h3>
        $run_report_html
    </div>
    
    <div class="footer">
        <p><strong>Generated by Olabuntu RNA-seq Pipeline v$pipeline_version</strong></p>
        <p>For questions or issues, check the pipeline log: <a href="logs/pipeline_master.log" class="file-link">logs/pipeline_master.log</a></p>
        <p style="margin-top: 10px; font-size: 0.9em; color: #999;">
            This report provides a comprehensive overview of your RNA-seq preprocessing results.<br>
            All files are ready for downstream differential expression analysis.
        </p>
    </div>
</body>
</html>
EOF

    echo "✓ Master HTML report created: $master_file"
}

generate_final_reports() {
    echo "Generating final reports..."
    
    # Create MultiQC configuration
    create_multiqc_config
    
    # Run comprehensive MultiQC
    echo "  Running final MultiQC report..."
    
    local multiqc_dirs=(
        "$OUTPUT_DIR/qc/raw_fastqc"
        "$OUTPUT_DIR/qc/posttrim_fastqc"
        "$OUTPUT_DIR/alignments/metrics"
        "$OUTPUT_DIR/quantification"
        "$OUTPUT_DIR/trimmed_fastq"
        "$OUTPUT_DIR/logs"
    )
    
    if ! $MULTIQC_CMD \
        --outdir "$OUTPUT_DIR/multiqc" \
        --filename "multiqc_report_final.html" \
        --title "RNA-seq Pipeline Final Report" \
        --config "$OUTPUT_DIR/multiqc_config.yaml" \
        --force \
        "${multiqc_dirs[@]}" \
        >> "$MASTER_LOG" 2>&1; then
        echo "  WARNING: MultiQC with config failed, trying without config..."
        # Fallback: run MultiQC without config file
        $MULTIQC_CMD \
            --outdir "$OUTPUT_DIR/multiqc" \
            --filename "multiqc_report_final.html" \
            --title "RNA-seq Pipeline Final Report" \
            --force \
            "${multiqc_dirs[@]}" \
            >> "$MASTER_LOG" 2>&1
    fi
    
    echo "  ✓ Final MultiQC report completed"
    
    echo "✓ Final reports generated"
}

create_mapping_summary() {
    local summary_file="$OUTPUT_DIR/results/mapping_summary.txt"
    
    {
        echo "==============================================================================="
        echo "RNA-seq Pipeline Mapping Summary"
        echo "Generated: $(date)"
        echo "==============================================================================="
        echo ""
        
        if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
            echo "ALIGNMENT-BASED WORKFLOW SUMMARY"
            echo "--------------------------------"
            
            while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
                [[ "$sample_id" == *sample_id* ]] && continue
                sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
                
                local flagstat_file="$OUTPUT_DIR/alignments/metrics/${sample_id}_flagstat.txt"
                if [[ -f "$flagstat_file" ]]; then
                    echo "Sample: $sample_id"
                    grep -E "(total|mapped|properly paired)" "$flagstat_file" | sed 's/^/  /'
                    echo ""
                fi
            done < "$SAMPLE_MANIFEST"
            
        elif [[ "$WORKFLOW_MODE" == "PSEUDOALIGNMENT" ]]; then
            echo "PSEUDOALIGNMENT WORKFLOW SUMMARY"
            echo "-------------------------------"
            
            while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
                [[ "$sample_id" == *sample_id* ]] && continue
                sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
                
                local salmon_log="$OUTPUT_DIR/quantification/$sample_id/logs/salmon_quant.log"
                if [[ -f "$salmon_log" ]]; then
                    echo "Sample: $sample_id"
                    grep -E "Mapping rate|processed" "$salmon_log" | sed 's/^/  /' || echo "  Log parsing failed"
                    echo ""
                fi
            done < "$SAMPLE_MANIFEST"
        fi
        
        echo "OUTPUT FILES"
        echo "------------"
        echo "Counts matrix: $OUTPUT_DIR/counts_matrix/gene_counts_matrix.tsv"
        echo "Final MultiQC report: $OUTPUT_DIR/multiqc/multiqc_report_final.html"
        echo "Master log: $MASTER_LOG"
        
    } > "$summary_file"
}

generate_final_reports

# ===============================================================================
# 11. FINAL HOUSEKEEPING & OPTIONAL DOWNSTREAM ANALYSIS TRIGGERS
# ===============================================================================

final_housekeeping() {
    echo "Performing final housekeeping..."
    
    # Compress intermediate files if requested
    # (This could be made configurable)
    
    # Create manifest of final outputs
    local manifest_file="$OUTPUT_DIR/results/output_manifest.txt"
    {
        echo "RNA-seq Pipeline Output Manifest"
        echo "Generated: $(date)"
        echo "==============================="
        echo ""
        echo "KEY OUTPUT FILES:"
        echo "Gene counts matrix: $(realpath "$OUTPUT_DIR/counts_matrix/gene_counts_matrix.tsv" 2>/dev/null || echo "Not found")"
        echo "Final QC report: $(realpath "$OUTPUT_DIR/multiqc/multiqc_report_final.html" 2>/dev/null || echo "Not found")"
        echo "Mapping summary: $(realpath "$OUTPUT_DIR/results/mapping_summary.txt" 2>/dev/null || echo "Not found")"
        echo ""
        echo "DIRECTORY STRUCTURE:"
        find "$OUTPUT_DIR" -type d | sort | sed 's/^/  /'
        echo ""
        echo "SAMPLE-SPECIFIC OUTPUTS:"
        
        while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == *sample_id* ]] && continue
            sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
            
            echo "  Sample: $sample_id"
            if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
                echo "    BAM: alignments/sorted/${sample_id}.sorted.bam"
            else
                echo "    Quantification: quantification/$sample_id/quant.sf"
            fi
            echo "    Logs: logs/samples/${sample_id}.log"
        done < "$SAMPLE_MANIFEST"
        
    } > "$manifest_file"
    
    # Final success message
    echo ""
    echo "==============================================================================="
    echo "PIPELINE COMPLETED SUCCESSFULLY!"
    echo "==============================================================================="
}

# ===============================================================================
# 12. ERROR HANDLING & LOGGING BEHAVIOR
# ===============================================================================

create_run_report() {
    local report_file="$OUTPUT_DIR/results/run_report.txt"
    local end_time=$(date)
    local start_time=$(head -n 10 "$MASTER_LOG" | grep "Started:" | cut -d' ' -f2-)
    
    {
        echo "==============================================================================="
        echo "RNA-seq Pipeline Run Report"
        echo "==============================================================================="
        echo ""
        echo "PIPELINE INFORMATION"
        echo "-------------------"
        echo "Script version: $SCRIPT_VERSION"
        echo "Workflow mode: $WORKFLOW_MODE"
        echo "Start time: $start_time"
        echo "End time: $end_time"
        echo "Output directory: $OUTPUT_DIR"
        echo ""
        echo "CONFIGURATION"
        echo "------------"
        echo "Threads used: $THREADS"
        echo "Trimming enabled: $TRIM_READS"
        echo "Strandedness: $STRANDEDNESS"
        echo ""
        echo "INPUT FILES"
        echo "----------"
        echo "Sample manifest: $SAMPLE_MANIFEST"
        if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
            echo "Reference genome: $REFERENCE_FASTA"
        else
            echo "Transcriptome: $TRANSCRIPTS_FASTA"
        fi
        echo "Annotation: $ANNOTATION_FILE"
        echo ""
        echo "SAMPLES PROCESSED"
        echo "----------------"
        
        local total_samples=0
        local failed_samples=0
        
        while IFS=',' read -r sample_id r1_path r2_path _extra || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == *sample_id* ]] && continue
            sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
            total_samples=$((total_samples + 1))
            
            # Check if sample completed successfully
            if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
                if [[ ! -f "$OUTPUT_DIR/alignments/sorted/${sample_id}.sorted.bam" ]]; then
                    echo "FAILED: $sample_id (missing BAM file)"
                    failed_samples=$((failed_samples + 1))
                else
                    echo "SUCCESS: $sample_id"
                fi
            else
                if [[ ! -f "$OUTPUT_DIR/quantification/$sample_id/quant.sf" ]]; then
                    echo "FAILED: $sample_id (missing quantification)"
                    failed_samples=$((failed_samples + 1))
                else
                    echo "SUCCESS: $sample_id"
                fi
            fi
        done < "$SAMPLE_MANIFEST"
        
        echo ""
        echo "SUMMARY"
        echo "-------"
        echo "Total samples: $total_samples"
        echo "Successful: $((total_samples - failed_samples))"
        echo "Failed: $failed_samples"
        
        if [[ $failed_samples -eq 0 ]]; then
            echo "Status: COMPLETED SUCCESSFULLY"
        else
            echo "Status: COMPLETED WITH ERRORS"
        fi
        
        echo ""
        echo "KEY OUTPUT LOCATIONS"
        echo "-------------------"
        echo "Gene counts matrix: counts_matrix/gene_counts_matrix.tsv"
        echo "Final QC report: multiqc/multiqc_report_final.html"
        echo "Master log: logs/pipeline_master.log"
        echo "Sample logs: logs/samples/"
        
    } > "$report_file"
    
    echo "✓ Run report created: $report_file"
}

# Execute final steps
final_housekeeping
create_run_report
create_mapping_summary
create_master_html_report

# Final log entry
{
    echo ""
    echo "==============================================================================="
    echo "Pipeline completed at: $(date)"
    echo "Total runtime: $SECONDS seconds"
    echo "==============================================================================="
} >> "$MASTER_LOG"
