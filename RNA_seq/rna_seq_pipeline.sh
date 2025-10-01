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
#   - R                (Statistical analysis: conda install -c conda-forge r-base)
#   - pandoc           (R Markdown rendering: brew install pandoc or conda install -c conda-forge pandoc)
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
#    conda install -c bioconda -c conda-forge fastqc multiqc fastp salmon hisat2 subread samtools pandas r-base pandoc

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
# ├── raw_fastq/                   ← Links to original FASTQ files
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
FASTP_CMD="fastp"               # Fast all-in-one preprocessing
HISAT2_CMD="hisat2"             # Primary aligner
HISAT2_BUILD_CMD="hisat2-build" # HISAT2 index building
SALMON_CMD="salmon"             # For pseudoalignment mode
FEATURECOUNTS_CMD="featureCounts"
SAMTOOLS_CMD="samtools"
R_CMD="R"                       # R for report generation

# RESOURCE SETTINGS
THREADS=8                       # Number of CPU threads to use

# PROCESSING OPTIONS
TRIM_READS="true"               # Set to "false" to skip trimming
MIN_READ_LENGTH=50              # Minimum read length after trimming
QUALITY_THRESHOLD=20            # Quality threshold for trimming
STRANDEDNESS="FR"               # HISAT2 strandedness: F, R, FR, RF

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
  --strandedness STRAND      HISAT2 strandedness: F, R, FR, RF (default: FR)
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
  - R: Statistical analysis and reporting
  - pandoc: R Markdown rendering

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
# Validate required inputs
# -------------------------------------------------------------------------------

validate_inputs() {
    local errors=0
    
    echo "Validating inputs..."
    
    # Check required parameters
    if [[ -z "$SAMPLE_MANIFEST" ]]; then
        echo "ERROR: Sample manifest (-m/--manifest) is required"
        errors=$((errors + 1))
    elif [[ ! -f "$SAMPLE_MANIFEST" ]]; then
        echo "ERROR: Sample manifest file does not exist: $SAMPLE_MANIFEST"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$OUTPUT_DIR" ]]; then
        echo "ERROR: Output directory (-o/--output) is required"
        errors=$((errors + 1))
    fi
    
    if [[ -z "$WORKFLOW_MODE" ]]; then
        echo "ERROR: Workflow mode (-w/--workflow) is required"
        errors=$((errors + 1))
    elif [[ "$WORKFLOW_MODE" != "ALIGNMENT" && "$WORKFLOW_MODE" != "PSEUDOALIGNMENT" ]]; then
        echo "ERROR: Workflow mode must be ALIGNMENT or PSEUDOALIGNMENT"
        errors=$((errors + 1))
    fi
    
    # Validate workflow-specific inputs
    if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
        if [[ -z "$REFERENCE_FASTA" ]]; then
            echo "ERROR: Reference genome (-g/--genome) is required for ALIGNMENT mode"
            errors=$((errors + 1))
        elif [[ ! -f "$REFERENCE_FASTA" ]]; then
            echo "ERROR: Reference genome file does not exist: $REFERENCE_FASTA"
            errors=$((errors + 1))
        fi
    elif [[ "$WORKFLOW_MODE" == "PSEUDOALIGNMENT" ]]; then
        if [[ -z "$TRANSCRIPTS_FASTA" ]]; then
            echo "ERROR: Transcripts FASTA (-t/--transcripts) is required for PSEUDOALIGNMENT mode"
            errors=$((errors + 1))
        elif [[ ! -f "$TRANSCRIPTS_FASTA" ]]; then
            echo "ERROR: Transcripts FASTA file does not exist: $TRANSCRIPTS_FASTA"
            errors=$((errors + 1))
        fi
    fi
    
    if [[ -z "$ANNOTATION_FILE" ]]; then
        echo "ERROR: Annotation file (-a/--annotation) is required"
        errors=$((errors + 1))
    elif [[ ! -f "$ANNOTATION_FILE" ]]; then
        echo "ERROR: Annotation file does not exist: $ANNOTATION_FILE"
        errors=$((errors + 1))
    fi
    
    # Check tool availability
    local tools=("$FASTQC_CMD" "$MULTIQC_CMD" "$SAMTOOLS_CMD" "$R_CMD" "pandoc")
    
    if [[ "$TRIM_READS" == "true" ]]; then
        tools+=("$FASTP_CMD")
    fi
    
    if [[ "$WORKFLOW_MODE" == "ALIGNMENT" ]]; then
        tools+=("$HISAT2_CMD" "$HISAT2_BUILD_CMD" "$FEATURECOUNTS_CMD")
    else
        tools+=("$SALMON_CMD")
    fi
    
    echo "Checking required tools..."
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            echo "ERROR: Required tool not found in PATH: $tool"
            errors=$((errors + 1))
        else
            echo "✓ Found: $tool"
        fi
    done
    
    if [[ $errors -gt 0 ]]; then
        echo "Found $errors error(s). Please fix them before running the pipeline."
        exit 1
    fi
    
    echo "✓ Input validation passed"
}

# Check R packages
check_r_packages() {
    echo "Checking R packages..."
    
    local required_packages=("rmarkdown" "knitr" "DT" "plotly")
    local missing_packages=()
    
    for package in "${required_packages[@]}"; do
        if ! Rscript -e "library($package)" &> /dev/null; then
            missing_packages+=("$package")
        fi
    done
    
    if [ ${#missing_packages[@]} -ne 0 ]; then
        echo "ERROR: Missing R packages: ${missing_packages[*]}"
        echo "Please install with: Rscript -e \"install.packages(c('${missing_packages[*]}'))\""
        echo "Or: conda install -c conda-forge r-rmarkdown r-knitr r-dt r-plotly"
        return 1
    fi
    
    echo "✓ All required R packages found"
    return 0
}

validate_inputs

# Check R packages (required for R Markdown)
if ! check_r_packages; then
    echo "ERROR: R packages are required for complete pipeline functionality."
    echo "Please install the missing packages and run the pipeline again."
    exit 1
fi

SKIP_R_REPORT=false

# ===============================================================================
# 2. PREPARE DIRECTORIES & LOGS
# ===============================================================================

setup_directories() {
    echo "Setting up directory structure..."
    
    # Convert to absolute path
    OUTPUT_DIR=$(realpath "$OUTPUT_DIR")
    
    # Create main output directory
    mkdir -p "$OUTPUT_DIR"
    
    # Create subdirectories
    mkdir -p "$OUTPUT_DIR"/{raw_fastq,qc/raw_fastqc,qc/posttrim_fastqc,multiqc,trimmed_fastq}
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

validate_manifest() {
    echo "Validating sample manifest..."
    
    # Check if manifest has header and proper format
    local header=$(head -n1 "$SAMPLE_MANIFEST")
    if [[ ! "$header" =~ sample_id.*r1.*r2 ]]; then
        echo "WARNING: Manifest header should contain 'sample_id', 'r1', 'r2' columns"
    fi
    
    # Read samples and validate files exist
    local sample_count=0
    local missing_files=0
    
    # Skip header and process each sample
    tail -n +2 "$SAMPLE_MANIFEST" | while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
        sample_count=$((sample_count + 1))
        
        # Remove any quotes and whitespace
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        r1_path=$(echo "$r1_path" | tr -d '"' | xargs)
        r2_path=$(echo "$r2_path" | tr -d '"' | xargs)
        
        # Check if files exist
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
    done
    
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
    tail -n +2 "$SAMPLE_MANIFEST" | while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
        # Clean sample data
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        r1_path=$(echo "$r1_path" | tr -d '"' | xargs)
        r2_path=$(echo "$r2_path" | tr -d '"' | xargs)
        
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
    done
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: FastQC failed for $failed_samples samples"
        exit 1
    fi
    
    # Run MultiQC on raw FastQC results
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
        tail -n +2 "$SAMPLE_MANIFEST" | while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
            sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
            r1_path=$(echo "$r1_path" | tr -d '"' | xargs)
            r2_path=$(echo "$r2_path" | tr -d '"' | xargs)
            
            ln -sf "$(realpath "$r1_path")" "$OUTPUT_DIR/trimmed_fastq/${sample_id}_R1.trim.fastq.gz"
            ln -sf "$(realpath "$r2_path")" "$OUTPUT_DIR/trimmed_fastq/${sample_id}_R2.trim.fastq.gz"
        done
        
        return
    fi
    
    echo "Running read trimming with fastp..."
    
    local trimmed_dir="$OUTPUT_DIR/trimmed_fastq"
    local failed_samples=0
    
    # Process each sample
    tail -n +2 "$SAMPLE_MANIFEST" | while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        r1_path=$(echo "$r1_path" | tr -d '"' | xargs)
        r2_path=$(echo "$r2_path" | tr -d '"' | xargs)
        
        echo "  Trimming sample: $sample_id"
        
        local r1_out="$trimmed_dir/${sample_id}_R1.trim.fastq.gz"
        local r2_out="$trimmed_dir/${sample_id}_R2.trim.fastq.gz"
        local trim_log="$SAMPLE_LOG_DIR/${sample_id}_trimming.log"
        
        # Run fastp
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
    done
    
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
    tail -n +2 "$SAMPLE_MANIFEST" | while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
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
    done
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: Post-trim FastQC failed for $failed_samples samples"
        exit 1
    fi
    
    # Run MultiQC on trimmed FastQC results
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
    tail -n +2 "$SAMPLE_MANIFEST" | while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
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
        
        # HISAT2 alignment command
        if ! $HISAT2_CMD \
            -x "$hisat2_index/genome" \
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
    done
    
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
    tail -n +2 "$SAMPLE_MANIFEST" | while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        
        echo "  Quantifying sample: $sample_id"
        
        local r1_input="$OUTPUT_DIR/trimmed_fastq/${sample_id}_R1.trim.fastq.gz"
        local r2_input="$OUTPUT_DIR/trimmed_fastq/${sample_id}_R2.trim.fastq.gz"
        local sample_quant_dir="$quant_dir/$sample_id"
        local quant_log="$SAMPLE_LOG_DIR/${sample_id}_quantification.log"
        
        # Determine library type based on strandedness
        local lib_type="A"  # Auto-detect
        case "$STRANDEDNESS" in
            "F") lib_type="ISF" ;;
            "R") lib_type="ISR" ;;
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
    done
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: Pseudo-alignment failed for $failed_samples samples"
        exit 1
    fi
    
    echo "✓ Pseudo-alignment completed"
}

# Run the appropriate mapping method
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
        while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == "sample_id" ]] && continue  # Skip header
            sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
            bam_files+=("$sorted_dir/${sample_id}.sorted.bam")
            sample_names+=("$sample_id")
        done < "$SAMPLE_MANIFEST"
        
        # Determine strandedness parameter for featureCounts
        local strand_param=0
        case "$STRANDEDNESS" in
            "F") strand_param=1 ;;
            "R") strand_param=2 ;;
            *) strand_param=0 ;;
        esac
        
        # Run featureCounts
        # Determine the correct attribute based on file format
        local gene_attr="gene_id"
        if [[ "$ANNOTATION_FILE" =~ \.gff3$ ]]; then
            gene_attr="ID"  # GFF3 uses ID attribute
        fi
        
        local featurecounts_args=(
            -p  # paired-end
            -T "$THREADS"
            -t exon
            -g "$gene_attr"
            -s "$strand_param"
            -a "$ANNOTATION_FILE"
            -o "$counts_dir/raw_gene_counts.txt"
            "${bam_files[@]}"
        )
        
        if ! $FEATURECOUNTS_CMD "${featurecounts_args[@]}" >> "$MASTER_LOG" 2>&1; then
            echo "ERROR: featureCounts failed"
            exit 1
        fi
        
        # Clean up the counts file header (remove path prefixes)
        if [[ -f "$counts_dir/raw_gene_counts.txt" ]]; then
            # Create a cleaner version with just sample names
            head -n1 "$counts_dir/raw_gene_counts.txt" | sed 's|[^[:space:]]*\/||g' > "$counts_dir/gene_counts_matrix.tsv"
            tail -n +2 "$counts_dir/raw_gene_counts.txt" | cut -f1,7- >> "$counts_dir/gene_counts_matrix.tsv"
            
            # Create DGE-ready CSV format
            echo "  Creating DGE-ready count matrix..."
            awk 'NR==1 {print "GeneID," substr($0, index($0, $2))} NR>1 {print $1"," substr($0, index($0, $2))}' "$counts_dir/gene_counts_matrix.tsv" | tr '\t' ',' > "$counts_dir/count_matrix.csv"
            
            # Create sample metadata for DGE analysis
            echo "  Creating sample metadata..."
            echo "sample_id,condition,replicate" > "$counts_dir/sample_metadata.csv"
            local replicate=1
            for sample in "${sample_names[@]}"; do
                echo "${sample},condition1,${replicate}" >> "$counts_dir/sample_metadata.csv"
                replicate=$((replicate + 1))
            done
            
            echo "  ✓ DGE-ready files created:"
            echo "    - count_matrix.csv (CSV format for DESeq2/edgeR)"
            echo "    - sample_metadata.csv (sample information)"
        fi
        
    elif [[ "$WORKFLOW_MODE" == "PSEUDOALIGNMENT" ]]; then
        echo "  Collecting Salmon quantification results..."
        
        # Create transcript-to-gene mapping from GTF/GFF
        local tx2gene_file="$counts_dir/tx2gene.tsv"
        create_tx2gene_mapping "$tx2gene_file"
        
        # Collect all quant.sf files
        local quant_files=()
        local sample_names=()
        
        while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == "sample_id" ]] && continue  # Skip header
            sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
            quant_files+=("$OUTPUT_DIR/quantification/$sample_id/quant.sf")
            sample_names+=("$sample_id")
        done < "$SAMPLE_MANIFEST"
        
        # Create a simple gene-level counts matrix (basic aggregation)
        # Note: For production use, consider using tximport in R
        create_gene_counts_matrix "$tx2gene_file" "${quant_files[@]}"
    fi
    
    echo "✓ Expression quantification completed"
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
    
    echo "  Creating gene-level counts matrix..."
    
    # This is a simplified aggregation - for production, use tximport
    local output_file="$OUTPUT_DIR/counts_matrix/gene_counts_matrix.tsv"
    
    # Create header
    printf "gene_id"
    for i in "${!sample_names[@]}"; do
        printf "\t%s" "${sample_names[$i]}"
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
sample_names = [$(printf "'%s'," "${sample_names[@]}" | sed 's/,$//'))]

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

# Create R Markdown template
create_rmarkdown_template() {
    local template_file="$OUTPUT_DIR/report.Rmd"
    
    cat > "$template_file" << 'EOF'
---
title: "RNA-seq Preprocessing Summary Report"
author: "RNA-seq Pipeline"
date: "`r Sys.Date()`"
output: 
  html_document:
    toc: true
    toc_float: true
    theme: flatly
    code_folding: show
    df_print: paged
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
library(DT)
library(plotly)
```

# Sample Information

```{r sample-info}
# Read sample metadata
sample_metadata <- read.csv("counts_matrix/sample_metadata.csv")
DT::datatable(sample_metadata, caption = "Sample Metadata")
```

# Quality Control Summary

## Raw Reads Quality
- **MultiQC Report**: [View Raw QC Report](multiqc/multiqc_report_raw.html)

## Trimmed Reads Quality  
- **MultiQC Report**: [View Trimmed QC Report](multiqc/multiqc_report_posttrim.html)

# Alignment Statistics

## HISAT2 Alignment Summary
- **MultiQC Report**: [View Alignment Report](multiqc/multiqc_report_final.html)

# Gene Expression Summary

```{r count-summary}
# Read count matrix
count_matrix <- read.csv("counts_matrix/count_matrix.csv", row.names = 1)

# Summary statistics
cat("Total genes:", nrow(count_matrix), "\n")
cat("Total samples:", ncol(count_matrix), "\n")
cat("Mean counts per gene:", round(mean(rowSums(count_matrix)), 2), "\n")
cat("Mean counts per sample:", round(mean(colSums(count_matrix)), 2), "\n")

# Display count matrix summary
summary_stats <- data.frame(
  Gene = rownames(count_matrix),
  Total_Counts = rowSums(count_matrix),
  Mean_Counts = rowMeans(count_matrix),
  Max_Counts = apply(count_matrix, 1, max)
)

DT::datatable(head(summary_stats, 100), caption = "Gene Count Summary (Top 100 genes)")
```

# Pipeline Information

- **Pipeline Version**: 1.0.0
- **Workflow Mode**: `r ifelse(exists("WORKFLOW_MODE"), WORKFLOW_MODE, "Unknown")`
- **Strandedness**: `r ifelse(exists("STRANDEDNESS"), STRANDEDNESS, "Unknown")`
- **Threads Used**: `r ifelse(exists("THREADS"), THREADS, "Unknown")`

# Next Steps

1. **Quality Assessment**: Review MultiQC reports for data quality
2. **Differential Expression**: Import count matrix into R for DGE analysis
3. **Statistical Analysis**: Use DESeq2, edgeR, or limma for analysis

EOF

    echo "✓ R Markdown template created: $template_file"
}

# Create master HTML report
create_master_html_report() {
    local master_file="$OUTPUT_DIR/master_report.html"
    
    cat > "$master_file" << 'EOF'
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>RNA-seq Preprocessing Master Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; }
        .header { background: #2c3e50; color: white; padding: 20px; border-radius: 5px; }
        .nav { background: #34495e; padding: 15px; margin: 20px 0; border-radius: 5px; }
        .nav a { color: white; text-decoration: none; margin-right: 20px; padding: 10px; }
        .nav a:hover { background: #3498db; border-radius: 3px; }
        .section { margin: 30px 0; padding: 20px; border: 1px solid #ddd; border-radius: 5px; }
        .footer { text-align: center; margin-top: 40px; color: #666; }
    </style>
</head>
<body>
    <div class="header">
        <h1>RNA-seq Preprocessing Master Report</h1>
        <p>Generated: <script>document.write(new Date().toLocaleString())</script></p>
    </div>
    
    <div class="nav">
        <a href="#overview">Overview</a>
        <a href="#multiqc">MultiQC Reports</a>
        <a href="#rmarkdown">R Analysis</a>
        <a href="#outputs">Output Files</a>
    </div>
    
    <div id="overview" class="section">
        <h2>Pipeline Overview</h2>
        <p>This report combines all RNA-seq preprocessing results into a single, comprehensive view.</p>
        <ul>
            <li><strong>Quality Control</strong>: FastQC analysis of raw and trimmed reads</li>
            <li><strong>Alignment</strong>: HISAT2 alignment statistics and mapping rates</li>
            <li><strong>Quantification</strong>: Gene expression count matrix</li>
            <li><strong>Analysis</strong>: R-based summaries and statistics</li>
        </ul>
    </div>
    
    <div id="multiqc" class="section">
        <h2>MultiQC Reports</h2>
        <p>Comprehensive quality control and alignment statistics:</p>
        <ul>
            <li><a href="multiqc/multiqc_report_raw.html" target="_blank">Raw Reads QC Report</a></li>
            <li><a href="multiqc/multiqc_report_posttrim.html" target="_blank">Trimmed Reads QC Report</a></li>
            <li><a href="multiqc/multiqc_report_final.html" target="_blank">Final Pipeline Report</a></li>
        </ul>
    </div>
    
    <div id="rmarkdown" class="section">
        <h2>R Analysis Report</h2>
        <p>Custom analysis and summaries generated with R Markdown:</p>
        <ul>
            <li><a href="custom_report.html" target="_blank">R Markdown Analysis Report</a></li>
        </ul>
    </div>
    
    <div id="outputs" class="section">
        <h2>Key Output Files</h2>
        <ul>
            <li><strong>Count Matrix (TSV)</strong>: <code>counts_matrix/gene_counts_matrix.tsv</code></li>
            <li><strong>Count Matrix (CSV)</strong>: <code>counts_matrix/count_matrix.csv</code></li>
            <li><strong>Sample Metadata</strong>: <code>counts_matrix/sample_metadata.csv</code></li>
            <li><strong>Pipeline Log</strong>: <code>logs/pipeline_master.log</code></li>
        </ul>
    </div>
    
    <div class="footer">
        <p>Generated by RNA-seq Pipeline v1.0.0</p>
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
    
    # Generate R Markdown report
    echo "  Generating R Markdown report..."
    create_rmarkdown_template
    
    # Render R Markdown to HTML
    cd "$OUTPUT_DIR"
    if Rscript -e "rmarkdown::render('report.Rmd', output_file='custom_report.html')" >> "$MASTER_LOG" 2>&1; then
        echo "  ✓ R Markdown report generated"
    else
        echo "  ERROR: R Markdown rendering failed"
        exit 1
    fi
    cd - > /dev/null
    
    # Create master HTML report
    create_master_html_report
    
    # Generate mapping summary
    echo "  Creating mapping summary..."
    create_mapping_summary
    
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
            
            while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
                [[ "$sample_id" == "sample_id" ]] && continue
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
            
            while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
                [[ "$sample_id" == "sample_id" ]] && continue
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
        
        while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == "sample_id" ]] && continue
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
        
        while IFS=',' read -r sample_id r1_path r2_path || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == "sample_id" ]] && continue
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

# Final log entry
{
    echo ""
    echo "==============================================================================="
    echo "Pipeline completed at: $(date)"
    echo "Total runtime: $SECONDS seconds"
    echo "==============================================================================="
} >> "$MASTER_LOG"
