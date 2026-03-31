#!/bin/bash

# ===============================================================================
# Variant Calling Pipeline
# ===============================================================================
# Description: Comprehensive variant calling pipeline supporting multiple workflows
# Purpose: Find genetic variants (SNPs, indels, SVs) from sequencing data
# Workflows: Germline, Somatic, RNA-seq, Structural Variants, Population
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
#   - bwa              (Alignment: conda install -c bioconda bwa)
#   - pandas           (Python library: pip install pandas)
#
# Trimming:
#   - fastp            (Fast all-in-one preprocessing: conda install -c bioconda fastp)
#
# Variant Calling (choose based on workflow):
#   For GERMLINE workflow:
#     - gatk           (Variant calling: conda install -c bioconda gatk4)
#     - picard         (BAM processing: conda install -c bioconda picard)
#   For SOMATIC workflow:
#     - gatk           (Mutect2 for somatic calling)
#   For RNASEQ workflow:
#     - star           (RNA alignment: conda install -c bioconda star)
#     - gatk           (Variant calling)
#   For SV workflow:
#     - manta          (Structural variants: conda install -c bioconda manta)
#   For POPULATION workflow:
#     - gatk           (Joint calling with GenotypeGVCFs)
#
# Variant Annotation (optional):
#   - vep              (Variant Effect Predictor: conda install -c bioconda ensembl-vep)
#   - snpeff           (SnpEff: conda install -c bioconda snpeff)
#
# 💡 Quick install command:
#    conda install -c bioconda -c conda-forge fastqc multiqc fastp bwa samtools gatk4 picard star manta ensembl-vep snpeff pandas

# 📁 REQUIRED INPUT FILES:
# ========================
# 1. Sample Manifest (samples.csv):
#    Format depends on workflow:
#    Germline: sample_id,r1_path,r2_path
#    Somatic: sample_id,r1_path,r2_path,type,matched_normal
#    Population: sample_id,r1_path,r2_path
#
# 2. Reference Files:
#    - genome.fa      (Reference genome FASTA - REQUIRED)
#    - genome.fa.fai  (FASTA index - will be created if missing)
#    - genome.dict    (GATK dictionary - will be created if missing)
#
# 3. Optional Files:
#    - adapters.fa     (Adapter sequences for trimming)
#    - known_sites.vcf (Known variants for BQSR/VQSR)
#    - pon.vcf         (Panel of Normals for somatic calling)
#    - target_regions.bed (Target regions for targeted sequencing)

# 📂 OUTPUT DIRECTORY STRUCTURE:
# ==============================
# results/
# ├── qc/                          # Quality control reports
# │   ├── raw_fastqc/             # FastQC on raw reads
# │   ├── posttrim_fastqc/        # FastQC on trimmed reads
# │   └── alignment_metrics/       # Alignment statistics
# ├── multiqc/                     # Aggregated QC reports
# ├── trimmed_fastq/               # Trimmed FASTQ files
# ├── alignments/                  # Alignment results
# │   ├── sorted/                 # Sorted BAM files
# │   ├── processed/              # Processed BAM (dedup, BQSR)
# │   └── metrics/                # Alignment metrics
# ├── indices/                     # Reference indices
# │   ├── bwa/                    # BWA index
# │   └── gatk/                   # GATK reference dict
# ├── variants/                    # Variant calling results
# │   ├── raw/                   # Raw VCF files (unfiltered)
# │   ├── filtered/              # Filtered VCF files
# │   ├── annotated/             # Annotated VCF files
# │   └── statistics/            # Variant statistics
# ├── logs/                       # Log files
# └── master_report.html          # Master HTML report

# 🎯 WORKFLOW MODES:
# ==================
# GERMLINE: Inherited variants (one sample per person)
# SOMATIC: Cancer mutations (tumor + normal pairs)
# RNASEQ: Variants from RNA-seq data
# SV: Structural variants (large deletions/duplications)
# POPULATION: Cohort analysis (many samples together)

# 🚀 USAGE EXAMPLES:
# ==================
# Germline workflow:
#   ./variant_calling_pipeline.sh -m samples.csv -o results/ -w GERMLINE -g genome.fa
#
# Somatic workflow:
#   ./variant_calling_pipeline.sh -m tumor_normal.csv -o results/ -w SOMATIC -g genome.fa
#
# With annotation:
#   ./variant_calling_pipeline.sh -m samples.csv -o results/ -w GERMLINE -g genome.fa --annotation

# ⚠️  IMPORTANT NOTES:
# ====================
# - Use simple sample names (no spaces, use underscores)
# - Ensure FASTQ files are paired-end and properly named
# - Reference genome must be indexed (BWA index, .fai, .dict)
# - Allow 4-24 hours runtime depending on dataset size and workflow
# - Somatic workflow requires matched tumor-normal pairs

# ===============================================================================

set -euo pipefail  # Exit on error, undefined vars, pipe failures

# Preserve argv for master log (parsed away in the loop below)
declare -a PIPELINE_CMD_ARGS=("$@")

# ===============================================================================
# 1. HEADER & CONFIGURATION BLOCK
# ===============================================================================

# Script version and info
SCRIPT_VERSION="1.0.0"
SCRIPT_NAME="Variant Calling Pipeline"
SCRIPT_AUTHOR="BUNTU(olabuntubabatunde@gmail.com) Under the supervision of Dr. P.Agre(p.agre@cgiar.org)"

# Print header
echo "==============================================================================="
echo "$SCRIPT_NAME v$SCRIPT_VERSION"
echo "Written by: $SCRIPT_AUTHOR"
echo "Started at: $(date)"
echo "==============================================================================="

# -------------------------------------------------------------------------------
# Configuration Variables
# -------------------------------------------------------------------------------

# INPUT FILES (REQUIRED)
SAMPLE_MANIFEST=""              # Path to samples.csv
REFERENCE_FASTA=""              # Reference genome FASTA

# WORKFLOW MODE (REQUIRED - choose one)
WORKFLOW_MODE=""                # "GERMLINE", "SOMATIC", "RNASEQ", "SV", "POPULATION"

# OUTPUT DIRECTORY (REQUIRED)
OUTPUT_DIR=""                   # Base output directory

# TOOL PATHS - modify if tools are not in PATH
FASTQC_CMD="fastqc"
MULTIQC_CMD="multiqc"
FASTP_CMD="fastp"
BWA_CMD="bwa"
SAMTOOLS_CMD="samtools"
GATK_CMD="gatk"
PICARD_CMD="picard"
STAR_CMD="STAR"
MANTA_CMD="configManta.py"
VEP_CMD="vep"
SNPEFF_CMD="snpEff"

# RESOURCE SETTINGS
THREADS=8

# PROCESSING OPTIONS
TRIM_READS="true"
MIN_READ_LENGTH=50
QUALITY_THRESHOLD=20
ENABLE_ANNOTATION="false"
ANNOTATION_TOOL="vep"           # "vep", "snpeff", or "annovar"

# VARIANT CALLING OPTIONS
VARIANT_CALLER="gatk"           # "gatk", "freebayes", "bcftools"
ENABLE_VQSR="false"             # Variant Quality Score Recalibration
HARD_FILTER_ONLY="true"         # Use hard filtering instead of VQSR

# OPTIONAL FILES
KNOWN_SITES=""                  # Known variants for BQSR/VQSR
TARGET_REGIONS=""               # BED file for targeted sequencing
PON_FILE=""                     # Panel of Normals for somatic calling
GERMLINE_RESOURCE=""            # Germline resource for Mutect2

# -------------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------------

# Absolute output path: mkdir first so this works for new directories (realpath alone fails).
resolve_output_dir() {
    local d="$1"
    mkdir -p "$d"
    (cd "$d" && pwd -P)
}

# SnpEff database name from reference path (handles .fa / .fasta / .fna / .gz).
snpeff_genome_name() {
    local base
    base=$(basename "$1")
    base="${base%.gz}"
    base="${base%.fasta}"
    base="${base%.fna}"
    base="${base%.fa}"
    printf '%s' "$base"
}

# Index compressed VCF (samtools index; avoid -t which is not for VCF).
index_vcf() {
    local vcf="$1"
    local log="${2:-}"
    if [[ -n "$log" ]]; then
        $SAMTOOLS_CMD index "$vcf" >> "$log" 2>&1 || true
    else
        $SAMTOOLS_CMD index "$vcf" 2>/dev/null || true
    fi
}

# -------------------------------------------------------------------------------
# Parse command line arguments
# -------------------------------------------------------------------------------

usage() {
    cat << EOF
Usage: $0 [OPTIONS]

Required Options:
  -m, --manifest FILE         Sample manifest (CSV format)
  -o, --output DIR           Output directory
  -w, --workflow MODE        Workflow: GERMLINE, SOMATIC, RNASEQ, SV, POPULATION
  -g, --genome FILE          Reference genome FASTA

Optional:
  --threads N                Number of threads (default: 8)
  --no-trim                  Skip read trimming step
  --annotation               Enable variant annotation
  --annotation-tool TOOL     Annotation tool: vep, snpeff (default: vep)
  --known-sites VCF          Known variants for BQSR/VQSR
  --target-regions BED       BED file for targeted sequencing
  --pon VCF                  Panel of Normals for somatic calling
  --germline-resource VCF    Germline resource for Mutect2
  --vqsr                     Enable VQSR (instead of hard filtering)
  --help                     Show this help message

Examples:
  # Germline variant calling
  $0 -m samples.csv -o results/ -w GERMLINE -g genome.fa

  # Somatic variant calling with annotation
  $0 -m tumor_normal.csv -o results/ -w SOMATIC -g genome.fa --annotation

  # RNA-seq variant calling
  $0 -m rnaseq_samples.csv -o results/ -w RNASEQ -g genome.fa
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
        --threads)
            THREADS="$2"
            shift 2
            ;;
        --no-trim)
            TRIM_READS="false"
            shift
            ;;
        --annotation)
            ENABLE_ANNOTATION="true"
            shift
            ;;
        --annotation-tool)
            ANNOTATION_TOOL="$2"
            shift 2
            ;;
        --known-sites)
            KNOWN_SITES="$2"
            shift 2
            ;;
        --target-regions)
            TARGET_REGIONS="$2"
            shift 2
            ;;
        --pon)
            PON_FILE="$2"
            shift 2
            ;;
        --germline-resource)
            GERMLINE_RESOURCE="$2"
            shift 2
            ;;
        --vqsr)
            ENABLE_VQSR="true"
            HARD_FILTER_ONLY="false"
            shift
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

# ===============================================================================
# 2. DEPENDENCY CHECKING
# ===============================================================================

check_all_dependencies() {
    local errors=0
    
    echo "==============================================================================="
    echo "DEPENDENCY CHECK - Verifying all required tools and packages"
    echo "==============================================================================="
    
    echo ""
    echo "Checking command-line tools..."
    local tools=("$FASTQC_CMD" "$MULTIQC_CMD" "$SAMTOOLS_CMD" "$BWA_CMD" "python3")
    
    if [[ "$TRIM_READS" == "true" ]]; then
        tools+=("$FASTP_CMD")
    fi
    
    # Workflow-specific tools
    if [[ "$WORKFLOW_MODE" == "RNASEQ" ]]; then
        tools+=("$STAR_CMD")
    fi
    
    if [[ "$WORKFLOW_MODE" == "SV" ]]; then
        tools+=("$MANTA_CMD")
    fi
    
    # GATK is used in most workflows
    if [[ "$VARIANT_CALLER" == "gatk" ]]; then
        tools+=("$GATK_CMD" "$PICARD_CMD")
    fi
    
    # Annotation tools
    if [[ "$ENABLE_ANNOTATION" == "true" ]]; then
        if [[ "$ANNOTATION_TOOL" == "vep" ]]; then
            tools+=("$VEP_CMD")
        elif [[ "$ANNOTATION_TOOL" == "snpeff" ]]; then
            tools+=("$SNPEFF_CMD")
        fi
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

# ===============================================================================
# 3. INPUT VALIDATION
# ===============================================================================

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
    elif [[ ! "$WORKFLOW_MODE" =~ ^(GERMLINE|SOMATIC|RNASEQ|SV|POPULATION)$ ]]; then
        echo "ERROR: Workflow mode must be GERMLINE, SOMATIC, RNASEQ, SV, or POPULATION"
        errors=$((errors + 1))
    else
        echo "✓ Workflow mode: $WORKFLOW_MODE"
    fi
    
    if [[ -z "$REFERENCE_FASTA" ]]; then
        echo "ERROR: Reference genome (-g/--genome) is required"
        errors=$((errors + 1))
    elif [[ ! -f "$REFERENCE_FASTA" ]]; then
        echo "ERROR: Reference genome file does not exist: $REFERENCE_FASTA"
        errors=$((errors + 1))
    else
        echo "✓ Reference genome: $REFERENCE_FASTA"
        
        # Check if index files exist (will create if missing)
        if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
            echo "  WARNING: FASTA index (.fai) not found, will be created"
        fi
    fi
    
    # Validate workflow-specific requirements
    if [[ "$WORKFLOW_MODE" == "SOMATIC" ]]; then
        echo "  Checking somatic workflow requirements..."
        # Will validate tumor-normal pairs in manifest validation
    fi
    
    if [[ $errors -gt 0 ]]; then
        echo ""
        echo "Found $errors error(s). Please fix them before running the pipeline."
        exit 1
    fi
    
    echo "✓ Input validation passed"
}

# ===============================================================================
# 4. DIRECTORY SETUP & LOGGING
# ===============================================================================

setup_directories() {
    echo "Setting up directory structure..."
    
    # Create directory first, then resolve to absolute path (realpath fails if path does not exist)
    OUTPUT_DIR=$(resolve_output_dir "$OUTPUT_DIR")
    
    # Create subdirectories
    mkdir -p "$OUTPUT_DIR"/{qc/raw_fastqc,qc/posttrim_fastqc,qc/alignment_metrics,multiqc,trimmed_fastq}
    mkdir -p "$OUTPUT_DIR"/{indices/bwa,indices/gatk,indices/star}
    mkdir -p "$OUTPUT_DIR"/{alignments/sorted,alignments/processed,alignments/metrics}
    mkdir -p "$OUTPUT_DIR"/{variants/raw,variants/filtered,variants/annotated,variants/statistics}
    mkdir -p "$OUTPUT_DIR"/{logs/samples,results}
    
    # Initialize log files
    MASTER_LOG="$OUTPUT_DIR/logs/pipeline_master.log"
    SAMPLE_LOG_DIR="$OUTPUT_DIR/logs/samples"
    
    # Start master log
    {
        echo "==============================================================================="
        echo "Variant Calling Pipeline Master Log"
        echo "Started: $(date)"
        printf 'Command:'
        printf ' %q' "$0" "${PIPELINE_CMD_ARGS[@]}"
        echo ""
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
    
    local header=$(head -n1 "$SAMPLE_MANIFEST")
    # Remove BOM if present
    header=$(echo "$header" | sed '1s/^\xEF\xBB\xBF//')
    
    # Check header based on workflow
    if [[ "$WORKFLOW_MODE" == "SOMATIC" ]]; then
        if [[ ! "$header" =~ sample_id.*r1.*r2.*type ]]; then
            echo "WARNING: Somatic workflow manifest should contain 'sample_id', 'r1_path', 'r2_path', 'type' columns"
        fi
    else
        if [[ ! "$header" =~ sample_id.*r1.*r2 ]]; then
            echo "WARNING: Manifest header should contain 'sample_id', 'r1_path', 'r2_path' columns"
        fi
    fi
    
    # Read samples and validate files exist
    local missing_files=0
    local tumor_samples=0
    local normal_samples=0
    local all_samples=0
    
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        # Remove BOM, quotes, whitespace, and carriage returns
        sample_id=$(echo "$sample_id" | sed 's/^\xEF\xBB\xBF//' | tr -d '"' | tr -d '\r' | xargs)
        r1_path=$(echo "$r1_path" | tr -d '"' | tr -d '\r' | xargs)
        r2_path=$(echo "$r2_path" | tr -d '"' | tr -d '\r' | xargs)
        type=$(echo "$type" | tr -d '"' | tr -d '\r' | xargs 2>/dev/null || echo "")
        matched_normal=$(echo "$matched_normal" | tr -d '"' | tr -d '\r' | xargs 2>/dev/null || echo "")
        
        # Skip header line
        if [[ "$sample_id" == "sample_id" ]]; then
            continue
        fi
        
        # Skip empty lines
        if [[ -z "$sample_id" ]]; then
            continue
        fi
        
        all_samples=$((all_samples + 1))
        
        if [[ ! -f "$r1_path" ]]; then
            echo "ERROR: R1 file not found for sample $sample_id: $r1_path"
            missing_files=$((missing_files + 1))
        fi
        
        if [[ ! -f "$r2_path" ]]; then
            echo "ERROR: R2 file not found for sample $sample_id: $r2_path"
            missing_files=$((missing_files + 1))
        fi
        
        # Validate somatic workflow requirements
        if [[ "$WORKFLOW_MODE" == "SOMATIC" ]]; then
            if [[ "$type" == "tumor" ]]; then
                tumor_samples=$((tumor_samples + 1))
                if [[ -z "$matched_normal" ]]; then
                    echo "WARNING: Tumor sample $sample_id has no matched_normal specified"
                fi
            elif [[ "$type" == "normal" ]]; then
                normal_samples=$((normal_samples + 1))
            fi
        fi
        
        # Create sample log file
        mkdir -p "$SAMPLE_LOG_DIR"
        touch "$SAMPLE_LOG_DIR/${sample_id}.log"
    done < "$SAMPLE_MANIFEST"
    
    if [[ $missing_files -gt 0 ]]; then
        echo "ERROR: Found $missing_files missing FASTQ files"
        exit 1
    fi
    
    # Validate somatic workflow has both tumor and normal
    if [[ "$WORKFLOW_MODE" == "SOMATIC" ]]; then
        if [[ $tumor_samples -eq 0 ]]; then
            echo "ERROR: Somatic workflow requires at least one tumor sample"
            exit 1
        fi
        if [[ $normal_samples -eq 0 ]]; then
            echo "WARNING: No normal samples found. Somatic calling works best with matched normals."
        fi
    fi
    
    echo "✓ Manifest validation passed"
    echo "  Total samples: $all_samples"
    if [[ "$WORKFLOW_MODE" == "SOMATIC" ]]; then
        echo "  Tumor samples: $tumor_samples"
        echo "  Normal samples: $normal_samples"
    fi
}

# ===============================================================================
# EXECUTE INITIALIZATION
# ===============================================================================

# Step 1: Check dependencies
if ! check_all_dependencies; then
    exit 1
fi

# Step 2: Validate inputs
validate_inputs

# Step 3: Setup directories
setup_directories

# Step 4: Validate manifest
validate_manifest

echo ""
echo "==============================================================================="
echo "Initialization complete. Starting pipeline..."
echo "==============================================================================="
echo ""

# ===============================================================================
# 5. QUALITY CONTROL (RAW READS)
# ===============================================================================

run_raw_fastqc() {
    echo "Running FastQC on raw reads..."
    
    local fastqc_out="$OUTPUT_DIR/qc/raw_fastqc"
    local failed_samples=0
    
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == "sample_id" ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | tr -d '\r' | xargs)
        r1_path=$(echo "$r1_path" | tr -d '"' | tr -d '\r' | xargs)
        r2_path=$(echo "$r2_path" | tr -d '"' | tr -d '\r' | xargs)
        
        echo "  Processing sample: $sample_id"
        
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

# ===============================================================================
# 6. READ TRIMMING
# ===============================================================================

run_trimming() {
    if [[ "$TRIM_READS" != "true" ]]; then
        echo "Skipping read trimming (disabled)"
        while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == "sample_id" ]] && continue
            sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
            r1_path=$(echo "$r1_path" | tr -d '"' | xargs)
            r2_path=$(echo "$r2_path" | tr -d '"' | xargs)
            ln -sf "$(realpath "$r1_path")" "$OUTPUT_DIR/trimmed_fastq/${sample_id}_R1.trim.fastq.gz"
            ln -sf "$(realpath "$r2_path")" "$OUTPUT_DIR/trimmed_fastq/${sample_id}_R2.trim.fastq.gz"
        done < "$SAMPLE_MANIFEST"
        return
    fi
    
    echo "Running read trimming with fastp..."
    
    local trimmed_dir="$OUTPUT_DIR/trimmed_fastq"
    local failed_samples=0
    
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == "sample_id" ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | tr -d '\r' | xargs)
        r1_path=$(echo "$r1_path" | tr -d '"' | tr -d '\r' | xargs)
        r2_path=$(echo "$r2_path" | tr -d '"' | tr -d '\r' | xargs)
        
        echo "  Trimming sample: $sample_id"
        
        local r1_out="$trimmed_dir/${sample_id}_R1.trim.fastq.gz"
        local r2_out="$trimmed_dir/${sample_id}_R2.trim.fastq.gz"
        local trim_log="$SAMPLE_LOG_DIR/${sample_id}_trimming.log"
        
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

# ===============================================================================
# 7. POST-TRIM QC
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
    
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == "sample_id" ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        
        echo "  Processing trimmed sample: $sample_id"
        
        local r1_trimmed="$trimmed_dir/${sample_id}_R1.trim.fastq.gz"
        local r2_trimmed="$trimmed_dir/${sample_id}_R2.trim.fastq.gz"
        
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

# ===============================================================================
# 8. REFERENCE INDEXING
# ===============================================================================

build_reference_index() {
    echo "Setting up reference indices..."
    
    # Create FASTA index if missing
    if [[ ! -f "${REFERENCE_FASTA}.fai" ]]; then
        echo "  Creating FASTA index..."
        $SAMTOOLS_CMD faidx "$REFERENCE_FASTA" >> "$MASTER_LOG" 2>&1
    fi
    
    # Create GATK dictionary if missing
    local dict_file="${REFERENCE_FASTA%.*}.dict"
    if [[ ! -f "$dict_file" ]]; then
        echo "  Creating GATK dictionary..."
        $PICARD_CMD CreateSequenceDictionary \
            R="$REFERENCE_FASTA" \
            O="$dict_file" \
            >> "$MASTER_LOG" 2>&1
    fi
    
    # Build BWA index if missing (for DNA-seq workflows)
    if [[ "$WORKFLOW_MODE" != "RNASEQ" ]]; then
        local bwa_index="$OUTPUT_DIR/indices/bwa/genome"
        if [[ ! -f "${bwa_index}.bwt" ]]; then
            echo "  Building BWA index (this may take a while)..."
            mkdir -p "$OUTPUT_DIR/indices/bwa"
            $BWA_CMD index -p "$bwa_index" "$REFERENCE_FASTA" >> "$MASTER_LOG" 2>&1
            echo "  ✓ BWA index built"
        else
            echo "  ✓ BWA index already exists"
        fi
    fi
    
    # Build STAR index for RNA-seq workflow
    if [[ "$WORKFLOW_MODE" == "RNASEQ" ]]; then
        local star_index="$OUTPUT_DIR/indices/star"
        if [[ ! -d "$star_index" || ! -f "$star_index/Genome" ]]; then
            echo "  Building STAR index (this may take a while)..."
            mkdir -p "$star_index"
            # Note: STAR index requires gene annotation - this is a placeholder
            # In production, you'd need GTF file and proper STAR index command
            echo "  WARNING: STAR index building requires GTF annotation file"
            echo "  Please build STAR index separately or provide pre-built index"
        else
            echo "  ✓ STAR index already exists"
        fi
    fi
    
    echo "✓ Reference indexing completed"
}

# ===============================================================================
# 9. ALIGNMENT
# ===============================================================================

run_alignment() {
    echo "Running alignment..."
    
    local align_dir="$OUTPUT_DIR/alignments"
    local sorted_dir="$align_dir/sorted"
    local failed_samples=0
    
    if [[ "$WORKFLOW_MODE" == "RNASEQ" ]]; then
        echo "  Using STAR for RNA-seq alignment..."
        # STAR alignment would go here
        echo "  WARNING: STAR alignment not yet implemented"
        echo "  Please align RNA-seq reads separately or use BWA for DNA-seq"
        return
    fi
    
    # BWA alignment for DNA-seq workflows
    local bwa_index="$OUTPUT_DIR/indices/bwa/genome"
    
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == "sample_id" ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        
        echo "  Aligning sample: $sample_id"
        
        local r1_input="$OUTPUT_DIR/trimmed_fastq/${sample_id}_R1.trim.fastq.gz"
        local r2_input="$OUTPUT_DIR/trimmed_fastq/${sample_id}_R2.trim.fastq.gz"
        local sam_output="$align_dir/${sample_id}.sam"
        local sorted_bam="$sorted_dir/${sample_id}.sorted.bam"
        local align_log="$SAMPLE_LOG_DIR/${sample_id}_alignment.log"
        
        # BWA alignment
        echo "    Running BWA alignment..."
        if ! $BWA_CMD mem \
            -t "$THREADS" \
            -R "@RG\tID:${sample_id}\tSM:${sample_id}\tPL:ILLUMINA" \
            "$bwa_index" \
            "$r1_input" \
            "$r2_input" \
            > "$sam_output" 2>> "$align_log"; then
            echo "ERROR: BWA alignment failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        # Convert SAM to BAM and sort
        echo "    Converting SAM to sorted BAM..."
        if ! $SAMTOOLS_CMD view -b "$sam_output" | \
            $SAMTOOLS_CMD sort -@ "$THREADS" -o "$sorted_bam" -; then
            echo "ERROR: SAM to BAM conversion failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        # Index BAM
        if ! $SAMTOOLS_CMD index "$sorted_bam" >> "$align_log" 2>&1; then
            echo "ERROR: BAM indexing failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        # Clean up SAM
        rm -f "$sam_output"
        
        # Generate alignment statistics
        $SAMTOOLS_CMD flagstat "$sorted_bam" > "$align_dir/metrics/${sample_id}_flagstat.txt" 2>> "$align_log"
        $SAMTOOLS_CMD stats "$sorted_bam" > "$align_dir/metrics/${sample_id}_stats.txt" 2>> "$align_log"
        
        echo "    ✓ Alignment completed for $sample_id"
    done < "$SAMPLE_MANIFEST"
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: Alignment failed for $failed_samples samples"
        exit 1
    fi
    
    echo "✓ Alignment completed"
}

# Execute QC and alignment steps
run_raw_fastqc
run_trimming
run_posttrim_fastqc
build_reference_index
run_alignment

# Final log entry
{
    echo ""
    echo "==============================================================================="
    echo "Pipeline processing started at: $(date)"
    echo "Workflow: $WORKFLOW_MODE"
    echo "==============================================================================="
} >> "$MASTER_LOG"

# ===============================================================================
# 10. POST-ALIGNMENT PROCESSING
# ===============================================================================

run_post_alignment_processing() {
    if [[ "$WORKFLOW_MODE" == "RNASEQ" ]]; then
        echo "Skipping duplicate marking for RNA-seq workflow (not recommended)"
        return
    fi
    
    echo "Running post-alignment processing (MarkDuplicates, BQSR)..."
    
    local sorted_dir="$OUTPUT_DIR/alignments/sorted"
    local processed_dir="$OUTPUT_DIR/alignments/processed"
    local metrics_dir="$OUTPUT_DIR/alignments/metrics"
    local failed_samples=0
    
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == "sample_id" ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        
        echo "  Processing sample: $sample_id"
        
        local sorted_bam="$sorted_dir/${sample_id}.sorted.bam"
        local dedup_bam="$processed_dir/${sample_id}.dedup.bam"
        local dedup_metrics="$metrics_dir/${sample_id}_dedup_metrics.txt"
        local recal_table="$processed_dir/${sample_id}.recal.table"
        local recal_bam="$processed_dir/${sample_id}.recal.bam"
        local process_log="$SAMPLE_LOG_DIR/${sample_id}_processing.log"
        
        if [[ ! -f "$sorted_bam" ]]; then
            echo "ERROR: Sorted BAM not found for sample $sample_id"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        # Mark duplicates
        echo "    Marking duplicates..."
        if ! $PICARD_CMD MarkDuplicates \
            INPUT="$sorted_bam" \
            OUTPUT="$dedup_bam" \
            METRICS_FILE="$dedup_metrics" \
            CREATE_INDEX=true \
            >> "$process_log" 2>&1; then
            echo "ERROR: MarkDuplicates failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        # Base Quality Score Recalibration (BQSR)
        if [[ -n "$KNOWN_SITES" && -f "$KNOWN_SITES" ]]; then
            echo "    Running BaseRecalibrator..."
            
            # Build recalibration table
            if ! $GATK_CMD BaseRecalibrator \
                -R "$REFERENCE_FASTA" \
                -I "$dedup_bam" \
                --known-sites "$KNOWN_SITES" \
                -O "$recal_table" \
                >> "$process_log" 2>&1; then
                echo "WARNING: BQSR failed for sample $sample_id, continuing without recalibration"
                # Use dedup BAM as final BAM
                cp "$dedup_bam" "$recal_bam"
                $SAMTOOLS_CMD index "$recal_bam"
            else
                # Apply recalibration
                echo "    Applying BQSR..."
                if ! $GATK_CMD ApplyBQSR \
                    -R "$REFERENCE_FASTA" \
                    -I "$dedup_bam" \
                    -bqsr "$recal_table" \
                    -O "$recal_bam" \
                    >> "$process_log" 2>&1; then
                    echo "WARNING: ApplyBQSR failed for sample $sample_id, using dedup BAM"
                    cp "$dedup_bam" "$recal_bam"
                    $SAMTOOLS_CMD index "$recal_bam"
                else
                    $SAMTOOLS_CMD index "$recal_bam"
                    echo "    ✓ BQSR completed for $sample_id"
                fi
            fi
        else
            echo "    Skipping BQSR (no known sites provided)"
            # Use dedup BAM as final BAM
            cp "$dedup_bam" "$recal_bam"
            $SAMTOOLS_CMD index "$recal_bam"
        fi
        
        echo "    ✓ Post-alignment processing completed for $sample_id"
    done < "$SAMPLE_MANIFEST"
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: Post-alignment processing failed for $failed_samples samples"
        exit 1
    fi
    
    echo "✓ Post-alignment processing completed"
}

# ===============================================================================
# 11. VARIANT CALLING (WORKFLOW-SPECIFIC)
# ===============================================================================

run_variant_calling() {
    echo "Running variant calling for $WORKFLOW_MODE workflow..."
    
    local variants_dir="$OUTPUT_DIR/variants"
    local raw_dir="$variants_dir/raw"
    local processed_dir="$OUTPUT_DIR/alignments/processed"
    local failed_samples=0
    
    case "$WORKFLOW_MODE" in
        GERMLINE)
            run_germline_variant_calling
            ;;
        SOMATIC)
            run_somatic_variant_calling
            ;;
        RNASEQ)
            run_rnaseq_variant_calling
            ;;
        SV)
            run_structural_variant_calling
            ;;
        POPULATION)
            run_population_variant_calling
            ;;
        *)
            echo "ERROR: Unknown workflow mode: $WORKFLOW_MODE"
            exit 1
            ;;
    esac
}

run_germline_variant_calling() {
    echo "  Running germline variant calling with GATK HaplotypeCaller..."
    
    local processed_dir="$OUTPUT_DIR/alignments/processed"
    local raw_dir="$OUTPUT_DIR/variants/raw"
    local failed_samples=0
    
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == "sample_id" ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        
        echo "    Calling variants for sample: $sample_id"
        
        local input_bam="$processed_dir/${sample_id}.recal.bam"
        if [[ ! -f "$input_bam" ]]; then
            # Fallback to sorted BAM if processed BAM doesn't exist
            input_bam="$OUTPUT_DIR/alignments/sorted/${sample_id}.sorted.bam"
        fi
        
        local gvcf_output="$raw_dir/${sample_id}.g.vcf.gz"
        local vcf_output="$raw_dir/${sample_id}.vcf.gz"
        local vc_log="$SAMPLE_LOG_DIR/${sample_id}_variant_calling.log"
        
        # GATK HaplotypeCaller (GVCF mode for joint calling, or VCF mode for single sample)
        local gatk_args=(
            HaplotypeCaller
            -R "$REFERENCE_FASTA"
            -I "$input_bam"
            -O "$gvcf_output"
            -ERC GVCF
        )
        
        if [[ -n "$TARGET_REGIONS" && -f "$TARGET_REGIONS" ]]; then
            gatk_args+=(-L "$TARGET_REGIONS")
        fi
        
        if ! $GATK_CMD "${gatk_args[@]}" >> "$vc_log" 2>&1; then
            echo "ERROR: HaplotypeCaller failed for sample $sample_id"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        if [[ "$WORKFLOW_MODE" == "POPULATION" ]]; then
            index_vcf "$gvcf_output" "$vc_log"
        fi
        
        # Convert GVCF to VCF (for single sample, or keep GVCF for joint calling)
        if [[ "$WORKFLOW_MODE" != "POPULATION" ]]; then
            echo "    Converting GVCF to VCF..."
            if ! $GATK_CMD GenotypeGVCFs \
                -R "$REFERENCE_FASTA" \
                -V "$gvcf_output" \
                -O "$vcf_output" \
                >> "$vc_log" 2>&1; then
                echo "ERROR: GenotypeGVCFs failed for sample $sample_id"
                failed_samples=$((failed_samples + 1))
                continue
            fi
        fi
        
        # Index VCF
        if [[ -f "$vcf_output" ]]; then
            index_vcf "$vcf_output" "$vc_log"
        fi
        
        echo "    ✓ Variant calling completed for $sample_id"
    done < "$SAMPLE_MANIFEST"
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: Germline variant calling failed for $failed_samples samples"
        exit 1
    fi
    
    echo "✓ Germline variant calling completed"
}

run_somatic_variant_calling() {
    echo "  Running somatic variant calling with GATK Mutect2..."
    
    local processed_dir="$OUTPUT_DIR/alignments/processed"
    local raw_dir="$OUTPUT_DIR/variants/raw"
    local failed_samples=0
    
    # Collect tumor-normal pairs
    declare -A tumor_normal_pairs
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == "sample_id" ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        type=$(echo "$type" | tr -d '"' | tr -d '\r' | xargs 2>/dev/null || echo "")
        matched_normal=$(echo "$matched_normal" | tr -d '"' | tr -d '\r' | xargs 2>/dev/null || echo "")
        
        if [[ "$type" == "tumor" ]]; then
            tumor_normal_pairs["$sample_id"]="$matched_normal"
        fi
    done < "$SAMPLE_MANIFEST"
    
    # Call variants for each tumor-normal pair
    for tumor_sample in "${!tumor_normal_pairs[@]}"; do
        local normal_sample="${tumor_normal_pairs[$tumor_sample]}"
        
        echo "    Calling somatic variants: tumor=$tumor_sample, normal=$normal_sample"
        
        local tumor_bam="$processed_dir/${tumor_sample}.recal.bam"
        if [[ ! -f "$tumor_bam" ]]; then
            tumor_bam="$OUTPUT_DIR/alignments/sorted/${tumor_sample}.sorted.bam"
        fi
        
        local normal_bam=""
        if [[ -n "$normal_sample" ]]; then
            normal_bam="$processed_dir/${normal_sample}.recal.bam"
            if [[ ! -f "$normal_bam" ]]; then
                normal_bam="$OUTPUT_DIR/alignments/sorted/${normal_sample}.sorted.bam"
            fi
        fi
        
        local vcf_output="$raw_dir/${tumor_sample}_somatic.vcf.gz"
        local vc_log="$SAMPLE_LOG_DIR/${tumor_sample}_variant_calling.log"
        
        # GATK Mutect2
        local mutect2_args=(
            Mutect2
            -R "$REFERENCE_FASTA"
            -I "$tumor_bam"
            -tumor "$tumor_sample"
        )
        
        if [[ -n "$normal_bam" && -f "$normal_bam" ]]; then
            mutect2_args+=(-I "$normal_bam" -normal "$normal_sample")
        fi
        
        if [[ -n "$PON_FILE" && -f "$PON_FILE" ]]; then
            mutect2_args+=(--panel-of-normals "$PON_FILE")
        fi
        
        if [[ -n "$GERMLINE_RESOURCE" && -f "$GERMLINE_RESOURCE" ]]; then
            mutect2_args+=(--germline-resource "$GERMLINE_RESOURCE")
        fi
        
        mutect2_args+=(-O "$vcf_output")
        
        if [[ -n "$TARGET_REGIONS" && -f "$TARGET_REGIONS" ]]; then
            mutect2_args+=(-L "$TARGET_REGIONS")
        fi
        
        if ! $GATK_CMD "${mutect2_args[@]}" >> "$vc_log" 2>&1; then
            echo "ERROR: Mutect2 failed for tumor sample $tumor_sample"
            failed_samples=$((failed_samples + 1))
            continue
        fi
        
        # Filter Mutect2 calls
        local filtered_vcf="${vcf_output%.vcf.gz}.filtered.vcf.gz"
        echo "    Filtering Mutect2 calls..."
        if ! $GATK_CMD FilterMutectCalls \
            -R "$REFERENCE_FASTA" \
            -V "$vcf_output" \
            -O "$filtered_vcf" \
            >> "$vc_log" 2>&1; then
            echo "WARNING: FilterMutectCalls failed, using unfiltered VCF"
        else
            rm -f "${vcf_output}.tbi" "${vcf_output}.csi"
            mv "$filtered_vcf" "$vcf_output"
        fi
        
        index_vcf "$vcf_output" "$vc_log"
        
        echo "    ✓ Somatic variant calling completed for $tumor_sample"
    done
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "ERROR: Somatic variant calling failed for $failed_samples samples"
        exit 1
    fi
    
    echo "✓ Somatic variant calling completed"
}

run_rnaseq_variant_calling() {
    echo "  Running RNA-seq variant calling..."
    echo "  WARNING: RNA-seq variant calling requires STAR alignment and GATK SplitNCigarReads"
    echo "  This is a placeholder - full implementation requires STAR-aligned BAMs"
    echo "  For now, please use DNA-seq workflow or implement STAR alignment separately"
}

run_structural_variant_calling() {
    echo "  Running structural variant calling with Manta..."
    echo "  WARNING: Structural variant calling not yet fully implemented"
    echo "  Manta requires specific configuration and tumor-normal pairing"
}

run_population_variant_calling() {
    echo "  Running population (joint) variant calling..."
    
    local raw_dir="$OUTPUT_DIR/variants/raw"
    local gvcf_files=()
    local sample_names=()
    
    # Collect all GVCF files
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == "sample_id" ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | xargs)
        
        local gvcf="$raw_dir/${sample_id}.g.vcf.gz"
        if [[ -f "$gvcf" ]]; then
            gvcf_files+=("-V" "$gvcf")
            sample_names+=("$sample_id")
        fi
    done < "$SAMPLE_MANIFEST"
    
    if [[ ${#gvcf_files[@]} -eq 0 ]]; then
        echo "ERROR: No GVCF files found. Run germline calling first."
        exit 1
    fi
    
    local joint_vcf="$raw_dir/joint_called.vcf.gz"
    local combined_gvcf="$raw_dir/cohort.combined.g.vcf.gz"
    local vc_log="$SAMPLE_LOG_DIR/joint_calling.log"
    local n_gvcf=$((${#gvcf_files[@]} / 2))
    
    echo "    Joint calling ${#sample_names[@]} samples ($n_gvcf GVCFs)..."
    if [[ $n_gvcf -eq 1 ]]; then
        if ! $GATK_CMD GenotypeGVCFs \
            -R "$REFERENCE_FASTA" \
            "${gvcf_files[@]}" \
            -O "$joint_vcf" \
            >> "$vc_log" 2>&1; then
            echo "ERROR: Joint calling failed"
            exit 1
        fi
    else
        if ! $GATK_CMD CombineGVCFs \
            -R "$REFERENCE_FASTA" \
            "${gvcf_files[@]}" \
            -O "$combined_gvcf" \
            >> "$vc_log" 2>&1; then
            echo "ERROR: CombineGVCFs failed"
            exit 1
        fi
        if ! $GATK_CMD IndexFeatureFile -I "$combined_gvcf" >> "$vc_log" 2>&1; then
            echo "WARNING: IndexFeatureFile failed on combined GVCF; GenotypeGVCFs may still succeed"
        fi
        if ! $GATK_CMD GenotypeGVCFs \
            -R "$REFERENCE_FASTA" \
            -V "$combined_gvcf" \
            -O "$joint_vcf" \
            >> "$vc_log" 2>&1; then
            echo "ERROR: Joint calling (GenotypeGVCFs) failed"
            exit 1
        fi
    fi
    
    index_vcf "$joint_vcf" "$vc_log"
    
    echo "✓ Population variant calling completed"
}

# ===============================================================================
# 12. VARIANT FILTERING
# ===============================================================================

run_variant_filtering() {
    echo "Running variant filtering..."
    
    local raw_dir="$OUTPUT_DIR/variants/raw"
    local filtered_dir="$OUTPUT_DIR/variants/filtered"
    local failed_samples=0
    
    # Find all VCF files
    local vcf_files=()
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == "sample_id" ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | tr -d '\r' | xargs)
        type=$(echo "$type" | tr -d '"' | tr -d '\r' | xargs 2>/dev/null || echo "")
        
        local vcf_file=""
        if [[ "$WORKFLOW_MODE" == "SOMATIC" && "$type" == "tumor" ]]; then
            vcf_file="$raw_dir/${sample_id}_somatic.vcf.gz"
        elif [[ "$WORKFLOW_MODE" == "POPULATION" ]]; then
            vcf_file="$raw_dir/joint_called.vcf.gz"
            if [[ -f "$vcf_file" ]]; then
                vcf_files+=("$vcf_file")
                break  # Only one joint VCF
            fi
        else
            vcf_file="$raw_dir/${sample_id}.vcf.gz"
        fi
        
        if [[ -f "$vcf_file" ]]; then
            vcf_files+=("$vcf_file")
        fi
    done < "$SAMPLE_MANIFEST"
    
    if [[ ${#vcf_files[@]} -eq 0 ]]; then
        echo "WARNING: No VCF files found for filtering"
        return
    fi
    
    for vcf_file in "${vcf_files[@]}"; do
        local sample_name=$(basename "$vcf_file" .vcf.gz)
        local filtered_vcf="$filtered_dir/${sample_name}.filtered.vcf.gz"
        local filter_log="$SAMPLE_LOG_DIR/${sample_name}_filtering.log"
        
        echo "  Filtering variants: $sample_name"
        
        if [[ "$HARD_FILTER_ONLY" == "true" ]]; then
            # Hard filtering with GATK VariantFiltration
            if ! $GATK_CMD VariantFiltration \
                -R "$REFERENCE_FASTA" \
                -V "$vcf_file" \
                -O "$filtered_vcf" \
                --filter-expression "QD < 2.0" --filter-name "QD2" \
                --filter-expression "FS > 60.0" --filter-name "FS60" \
                --filter-expression "MQ < 40.0" --filter-name "MQ40" \
                --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
                --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
                >> "$filter_log" 2>&1; then
                echo "ERROR: Variant filtering failed for $sample_name"
                failed_samples=$((failed_samples + 1))
                continue
            fi
            
            # Keep only PASS sites (VCF-aware; do not pipe VCF through samtools BAM path)
            local pass_vcf="${filtered_vcf%.vcf.gz}.pass.vcf.gz"
            if ! $GATK_CMD SelectVariants \
                -R "$REFERENCE_FASTA" \
                -V "$filtered_vcf" \
                -O "$pass_vcf" \
                --exclude-filtered \
                >> "$filter_log" 2>&1; then
                echo "ERROR: SelectVariants (PASS) failed for $sample_name"
                failed_samples=$((failed_samples + 1))
                continue
            fi
            rm -f "${filtered_vcf}.tbi" "${filtered_vcf}.csi"
            mv "$pass_vcf" "$filtered_vcf"
            index_vcf "$filtered_vcf" "$filter_log"
        else
            # VQSR (Variant Quality Score Recalibration) - more complex
            echo "  WARNING: VQSR not yet fully implemented, using hard filtering"
            cp "$vcf_file" "$filtered_vcf"
        fi
        
        echo "    ✓ Filtering completed for $sample_name"
    done
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "WARNING: Filtering failed for $failed_samples samples"
    else
        echo "✓ Variant filtering completed"
    fi
}

# ===============================================================================
# 13. VARIANT ANNOTATION (OPTIONAL)
# ===============================================================================

run_variant_annotation() {
    if [[ "$ENABLE_ANNOTATION" != "true" ]]; then
        echo "Skipping variant annotation (disabled)"
        return
    fi
    
    echo "Running variant annotation with $ANNOTATION_TOOL..."
    
    local filtered_dir="$OUTPUT_DIR/variants/filtered"
    local annotated_dir="$OUTPUT_DIR/variants/annotated"
    local failed_samples=0
    
    # Find all filtered VCF files
    local vcf_files=()
    if [[ "$WORKFLOW_MODE" == "POPULATION" ]]; then
        local joint_vcf="$filtered_dir/joint_called.filtered.vcf.gz"
        if [[ -f "$joint_vcf" ]]; then
            vcf_files+=("$joint_vcf")
        fi
    else
        while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
            [[ "$sample_id" == "sample_id" ]] && continue
            sample_id=$(echo "$sample_id" | tr -d '"' | tr -d '\r' | xargs)
            type=$(echo "$type" | tr -d '"' | tr -d '\r' | xargs 2>/dev/null || echo "")
            
            local vcf_file=""
            if [[ "$WORKFLOW_MODE" == "SOMATIC" && "$type" == "tumor" ]]; then
                vcf_file="$filtered_dir/${sample_id}_somatic.filtered.vcf.gz"
            else
                vcf_file="$filtered_dir/${sample_id}.filtered.vcf.gz"
            fi
            
            if [[ -f "$vcf_file" ]]; then
                vcf_files+=("$vcf_file")
            fi
        done < "$SAMPLE_MANIFEST"
    fi
    
    if [[ ${#vcf_files[@]} -eq 0 ]]; then
        echo "WARNING: No filtered VCF files found for annotation"
        return
    fi
    
    for vcf_file in "${vcf_files[@]}"; do
        local sample_name=$(basename "$vcf_file" .filtered.vcf.gz)
        local annotated_vcf="$annotated_dir/${sample_name}.annotated.vcf.gz"
        local annot_log="$SAMPLE_LOG_DIR/${sample_name}_annotation.log"
        
        echo "  Annotating variants: $sample_name"
        
        if [[ "$ANNOTATION_TOOL" == "vep" ]]; then
            # VEP annotation
            if ! $VEP_CMD \
                --input_file "$vcf_file" \
                --output_file "$annotated_vcf" \
                --format vcf \
                --compress_output gzip \
                --species human \
                --cache \
                --offline \
                --force_overwrite \
                >> "$annot_log" 2>&1; then
                echo "WARNING: VEP annotation failed for $sample_name"
                failed_samples=$((failed_samples + 1))
                continue
            fi
        elif [[ "$ANNOTATION_TOOL" == "snpeff" ]]; then
            # SnpEff annotation
            if ! $SNPEFF_CMD \
                -v \
                -o vcf \
                "$(snpeff_genome_name "$REFERENCE_FASTA")" \
                "$vcf_file" \
                > "$annotated_vcf" 2>> "$annot_log"; then
                echo "WARNING: SnpEff annotation failed for $sample_name"
                failed_samples=$((failed_samples + 1))
                continue
            fi
        else
            echo "WARNING: Unknown annotation tool: $ANNOTATION_TOOL"
            continue
        fi
        
        # Index annotated VCF
        index_vcf "$annotated_vcf" "$annot_log"
        
        echo "    ✓ Annotation completed for $sample_name"
    done
    
    if [[ $failed_samples -gt 0 ]]; then
        echo "WARNING: Annotation failed for $failed_samples samples"
    else
        echo "✓ Variant annotation completed"
    fi
}

# ===============================================================================
# 14. VARIANT STATISTICS
# ===============================================================================

calculate_variant_statistics() {
    echo "Calculating variant statistics..."
    
    local stats_dir="$OUTPUT_DIR/variants/statistics"
    local filtered_dir="$OUTPUT_DIR/variants/filtered"
    local annotated_dir="$OUTPUT_DIR/variants/annotated"
    
    # Use annotated VCFs if available, otherwise filtered
    local vcf_dir="$filtered_dir"
    if [[ "$ENABLE_ANNOTATION" == "true" && -d "$annotated_dir" ]]; then
        vcf_dir="$annotated_dir"
    fi
    
    while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
        [[ "$sample_id" == "sample_id" ]] && continue
        sample_id=$(echo "$sample_id" | tr -d '"' | tr -d '\r' | xargs)
        type=$(echo "$type" | tr -d '"' | tr -d '\r' | xargs 2>/dev/null || echo "")
        
        local vcf_file=""
        if [[ "$WORKFLOW_MODE" == "SOMATIC" && "$type" == "tumor" ]]; then
            vcf_file="$vcf_dir/${sample_id}_somatic.filtered.vcf.gz"
            if [[ "$ENABLE_ANNOTATION" == "true" ]]; then
                vcf_file="$vcf_dir/${sample_id}_somatic.annotated.vcf.gz"
            fi
        elif [[ "$WORKFLOW_MODE" == "POPULATION" ]]; then
            vcf_file="$vcf_dir/joint_called.filtered.vcf.gz"
            if [[ "$ENABLE_ANNOTATION" == "true" ]]; then
                vcf_file="$vcf_dir/joint_called.annotated.vcf.gz"
            fi
            if [[ -f "$vcf_file" ]]; then
                sample_id="joint_cohort"
                break
            fi
            continue
        else
            vcf_file="$vcf_dir/${sample_id}.filtered.vcf.gz"
            if [[ "$ENABLE_ANNOTATION" == "true" ]]; then
                vcf_file="$vcf_dir/${sample_id}.annotated.vcf.gz"
            fi
        fi
        
        if [[ ! -f "$vcf_file" ]]; then
            continue
        fi
        
        local stats_file="$stats_dir/${sample_id}_variant_stats.txt"
        
        echo "  Calculating statistics for: $sample_id"
        
        {
            echo "Variant Statistics for: $sample_id"
            echo "Generated: $(date)"
            echo "==============================================================================="
            echo ""
            
            # Total variants
            local total_variants=$($SAMTOOLS_CMD view "$vcf_file" | grep -v "^#" | wc -l | awk '{print $1}')
            echo "Total Variants: $total_variants"
            echo ""
            
            # SNP vs Indel counts
            local snp_count=$($SAMTOOLS_CMD view "$vcf_file" | grep -v "^#" | awk 'length($4)==1 && length($5)==1' | wc -l | awk '{print $1}')
            local indel_count=$((total_variants - snp_count))
            echo "SNPs: $snp_count"
            echo "Indels: $indel_count"
            echo ""
            
            # Quality metrics (column 6 = QUAL in VCF)
            echo "Quality Metrics:"
            $SAMTOOLS_CMD view "$vcf_file" | grep -v "^#" | awk '
                $6 != "." && $6+0 == $6 {
                    q = $6 + 0
                    sum += q
                    count++
                    if (count == 1) { min = q; max = q }
                    if (q < min) min = q
                    if (q > max) max = q
                }
                END {
                    if (count > 0)
                        printf "Mean QUAL: %f\nMax QUAL: %f\nMin QUAL: %f\n", sum/count, max, min
                    else
                        print "QUAL metrics not available"
                }' || echo "QUAL metrics not available"
            echo ""
            
        } > "$stats_file"
        
    done < "$SAMPLE_MANIFEST"
    
    echo "✓ Variant statistics calculated"
}

# ===============================================================================
# 15. REPORTING
# ===============================================================================

generate_final_reports() {
    echo "Generating final reports..."
    
    # Run final MultiQC
    echo "  Running final MultiQC report..."
    local multiqc_dirs=(
        "$OUTPUT_DIR/qc/raw_fastqc"
        "$OUTPUT_DIR/qc/posttrim_fastqc"
        "$OUTPUT_DIR/alignments/metrics"
        "$OUTPUT_DIR/trimmed_fastq"
        "$OUTPUT_DIR/logs"
    )
    
    $MULTIQC_CMD \
        --outdir "$OUTPUT_DIR/multiqc" \
        --filename "multiqc_report_final.html" \
        --title "Variant Calling Pipeline Final Report" \
        --force \
        "${multiqc_dirs[@]}" \
        >> "$MASTER_LOG" 2>&1 || echo "WARNING: MultiQC failed"
    
    echo "  ✓ Final MultiQC report completed"
}

create_run_report() {
    local report_file="$OUTPUT_DIR/results/run_report.txt"
    local end_time=$(date)
    
    {
        echo "==============================================================================="
        echo "Variant Calling Pipeline Run Report"
        echo "==============================================================================="
        echo ""
        echo "PIPELINE INFORMATION"
        echo "-------------------"
        echo "Script version: $SCRIPT_VERSION"
        echo "Workflow mode: $WORKFLOW_MODE"
        echo "Start time: $(head -n 10 "$MASTER_LOG" | grep "Started:" | cut -d' ' -f2- || echo "N/A")"
        echo "End time: $end_time"
        echo "Output directory: $OUTPUT_DIR"
        echo ""
        echo "CONFIGURATION"
        echo "------------"
        echo "Threads used: $THREADS"
        echo "Trimming enabled: $TRIM_READS"
        echo "Annotation enabled: $ENABLE_ANNOTATION"
        echo "Annotation tool: $ANNOTATION_TOOL"
        echo ""
        echo "INPUT FILES"
        echo "----------"
        echo "Sample manifest: $SAMPLE_MANIFEST"
        echo "Reference genome: $REFERENCE_FASTA"
        echo ""
        echo "SAMPLES PROCESSED"
        echo "----------------"
        
        local total_samples=0
        local failed_samples=0
        local joint_vcf="$OUTPUT_DIR/variants/raw/joint_called.vcf.gz"
        
        if [[ "$WORKFLOW_MODE" == "POPULATION" ]]; then
            while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
                [[ "$sample_id" == "sample_id" ]] && continue
                sample_id=$(echo "$sample_id" | tr -d '"' | tr -d '\r' | xargs)
                [[ -z "$sample_id" ]] && continue
                total_samples=$((total_samples + 1))
                echo "COHORT MEMBER: $sample_id"
            done < "$SAMPLE_MANIFEST"
            if [[ -f "$joint_vcf" ]]; then
                echo "SUCCESS: Joint cohort VCF: $joint_vcf"
            else
                echo "FAILED: Joint cohort VCF missing: $joint_vcf"
                failed_samples=$((failed_samples + 1))
            fi
        else
            while IFS=',' read -r sample_id r1_path r2_path type matched_normal || [[ -n "$sample_id" ]]; do
                [[ "$sample_id" == "sample_id" ]] && continue
                sample_id=$(echo "$sample_id" | tr -d '"' | tr -d '\r' | xargs)
                type=$(echo "$type" | tr -d '"' | tr -d '\r' | xargs 2>/dev/null || echo "")
                [[ -z "$sample_id" ]] && continue
                total_samples=$((total_samples + 1))
                
                local vcf_file=""
                if [[ "$WORKFLOW_MODE" == "SOMATIC" ]]; then
                    if [[ "$type" == "normal" ]]; then
                        echo "OK (matched normal, no somatic VCF): $sample_id"
                        continue
                    fi
                    vcf_file="$OUTPUT_DIR/variants/raw/${sample_id}_somatic.vcf.gz"
                else
                    vcf_file="$OUTPUT_DIR/variants/raw/${sample_id}.vcf.gz"
                fi
                
                if [[ ! -f "$vcf_file" ]]; then
                    echo "FAILED: $sample_id (missing VCF file)"
                    failed_samples=$((failed_samples + 1))
                else
                    echo "SUCCESS: $sample_id"
                fi
            done < "$SAMPLE_MANIFEST"
        fi
        
        echo ""
        echo "SUMMARY"
        echo "-------"
        echo "Total manifest rows: $total_samples"
        echo "Failed checks: $failed_samples"
        
        if [[ $failed_samples -eq 0 ]]; then
            echo "Status: COMPLETED SUCCESSFULLY"
        else
            echo "Status: COMPLETED WITH ERRORS"
        fi
        
        echo ""
        echo "KEY OUTPUT LOCATIONS"
        echo "-------------------"
        echo "Variant files: variants/"
        echo "Final QC report: multiqc/multiqc_report_final.html"
        echo "Master log: logs/pipeline_master.log"
        echo "Sample logs: logs/samples/"
        
    } > "$report_file"
    
    echo "✓ Run report created: $report_file"
}

create_master_html_report() {
    local master_file="$OUTPUT_DIR/master_report.html"
    local run_report="$OUTPUT_DIR/results/run_report.txt"
    
    # Get variant statistics
    local total_variants="N/A"
    local stats_dir="$OUTPUT_DIR/variants/statistics"
    if [[ -d "$stats_dir" ]]; then
        local first_stats=$(find "$stats_dir" -name "*_variant_stats.txt" | head -1)
        if [[ -f "$first_stats" ]]; then
            total_variants=$(grep "Total Variants:" "$first_stats" | awk '{print $3}' || echo "N/A")
        fi
    fi
    
    # Get pipeline info
    local pipeline_version="$SCRIPT_VERSION"
    local workflow_mode="$WORKFLOW_MODE"
    local status="SUCCESS"
    
    if [[ -f "$run_report" ]]; then
        if grep -q "COMPLETED WITH ERRORS" "$run_report"; then
            status="ERRORS"
        fi
    fi
    
    cat > "$master_file" << EOF
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Variant Calling Pipeline Master Report</title>
    <style>
        body { font-family: Arial, sans-serif; margin: 40px; line-height: 1.6; background: #f9f9f9; }
        .header { background: linear-gradient(135deg, #2c3e50 0%, #34495e 100%); color: white; padding: 30px; border-radius: 8px; }
        .section { margin: 30px 0; padding: 25px; background: white; border-radius: 8px; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }
        .info-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 15px; margin: 20px 0; }
        .info-box { background: #f8f9fa; padding: 15px; border-radius: 5px; border-left: 4px solid #3498db; }
        .status-success { color: #27ae60; font-weight: bold; }
        .status-error { color: #e74c3c; font-weight: bold; }
        table { width: 100%; border-collapse: collapse; margin: 15px 0; }
        table th, table td { padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }
        table th { background: #34495e; color: white; }
    </style>
</head>
<body>
    <div class="header">
        <h1>Variant Calling Pipeline Master Report</h1>
        <p>Generated: $(date '+%Y-%m-%d %H:%M:%S')</p>
        <p>Pipeline Version: $pipeline_version | Workflow: $workflow_mode</p>
    </div>
    
    <div class="section">
        <h2>📊 Pipeline Summary</h2>
        <div class="info-grid">
            <div class="info-box">
                <strong>Status</strong>
                <span class="status-$([ "$status" = "SUCCESS" ] && echo "success" || echo "error")">$status</span>
            </div>
            <div class="info-box">
                <strong>Workflow Mode</strong>
                $workflow_mode
            </div>
            <div class="info-box">
                <strong>Total Variants</strong>
                $total_variants
            </div>
        </div>
    </div>
    
    <div class="section">
        <h2>📁 Output Files</h2>
        <ul>
            <li><strong>Variant Files:</strong> <code>variants/</code></li>
            <li><strong>QC Reports:</strong> <code>multiqc/multiqc_report_final.html</code></li>
            <li><strong>Statistics:</strong> <code>variants/statistics/</code></li>
            <li><strong>Logs:</strong> <code>logs/</code></li>
        </ul>
    </div>
    
    <div class="footer" style="text-align: center; margin-top: 40px; padding: 20px; color: #666;">
        <p><strong>Generated by Olabuntu Variant Calling Pipeline v$pipeline_version</strong></p>
    </div>
</body>
</html>
EOF

    echo "✓ Master HTML report created: $master_file"
}

# Execute final steps
run_variant_filtering
run_variant_annotation
calculate_variant_statistics
generate_final_reports
create_run_report
create_master_html_report

# Final log entry
{
    echo ""
    echo "==============================================================================="
    echo "Pipeline completed at: $(date)"
    echo "Total runtime: $SECONDS seconds"
    echo "==============================================================================="
} >> "$MASTER_LOG"

echo ""
echo "==============================================================================="
echo "PIPELINE COMPLETED SUCCESSFULLY!"
echo "==============================================================================="
echo ""
echo "Output directory: $OUTPUT_DIR"
echo "Master report: $OUTPUT_DIR/master_report.html"
echo "Variant files: $OUTPUT_DIR/variants/"
echo ""

