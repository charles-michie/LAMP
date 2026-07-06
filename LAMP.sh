#!/bin/bash
# =============================================================================
# LAMP: LTR Annotation and Mining Pipeline
# =============================================================================
#
# DESCRIPTION:
#   This pipeline identifies and characterises LTR retrotransposons from a
#   RepeatMasker 2 (RM2) Stockholm file. It:
#     1. Parses the Stockholm file to find co-located LTR and INT repeat pairs
#     2. Classifies families by whether they have matched LTR-INT pairs
#     3. BLASTs LTR sequences against the genome to find all insertion sites
#     4. Searches flanking regions for conserved protein domains (CDD/rpstblastn)
#     5. Produces summary tables per category (ltr_matched, ltr_unmatched, int_unmatched)
#
# USAGE:
#   bash LAMP.sh <genome_fasta> <query_ltrs_fasta> <stockholm_file> <cdd_db_dir> [flank_size]
#
# ARGUMENTS:
#   genome_fasta      Path to the genome assembly FASTA file (any number/naming of chromosomes)
#   query_ltrs_fasta  Path to the RM2 LTR family consensus FASTA (query sequences for BLAST)
#   stockholm_file    Path to the RM2 Stockholm (.stk) alignment file
#   cdd_db_dir        Path to the CDD (Conserved Domain Database) for rpstblastn
#   flank_size        (Optional) Flanking region size in bp around BLAST hits (default: 14000)
#
# OUTPUTS (all in a timestamped directory):
#   matched_families_yes.txt                   - Families with at least one LTR-INT pair
#   matched_families_no.txt                    - Families with no LTR-INT pairs
#   family_match_results.csv                   - Family-level LTR-INT match pairs
#   ltr_int_matches.csv                        - Sequence-level LTR-INT pair details
#   all_ltr_int_sequences.csv                  - All sequences with match status (Yes/No)
#   all_hits_combined.txt                      - Combined BLAST results (all categories)
#   all_hits_deduplicated.txt                  - BLAST results after removing overlapping hits
#   *_blast_extended_sorted.txt                - Per-category BLAST hits with strand + sorted
#   *_flank_${flank_size}bp.bed/.fa            - Flanking regions in BED and FASTA format
#   *_cdd_summary.csv                          - CDD domain hits per flanking region (with human-readable names)
#   *_cdd_hits.gff                             - CDD hits in GFF format (for Geneious etc.)
#                                                GFF attributes include both the original BLAST hit
#                                                coordinates AND the human-readable domain name
#                                                (e.g. RVT_1, Gag_p24, rve, ENV) so ERV structure
#                                                (LTR-GAG-POL-ENV-LTR) is immediately visible
#   *_family_sequences_summary.tsv             - Per-category summary table
#
# DEPENDENCIES:
#   python3, blastn, makeblastdb, rpstblastn, samtools, bedtools, gawk, xargs
#
# NOTES:
#   - Sequences on scaffolds prefixed with "NW_" are excluded from BLAST analysis.
#     Edit the filter in run_blast_pipeline() if your genome uses different scaffold naming.
#   - The deduplication step keeps the highest percent-identity hit within a 10 bp window.
#   - CDD search uses e-value 1e-3 and up to 40 target sequences per query.
#   - rpstblastn calls are parallelised across families using xargs -P (default 8 parallel jobs).
#     Adjust CDD_PARALLEL_JOBS below if your system has more/fewer cores available.
#   - CDD accessions are mapped to human-readable names via a single batched query
#     to the NCBI CDD esummary API, made after all rpstblastn jobs have finished.
#     All unique accessions across all three categories are collected first, then
#     resolved in one go (batches of 200, rate-limited to 3 req/s). Domain names
#     like RVT_1 (reverse transcriptase), Gag_p24, rve (integrase), and ENV are
#     written into GFF attributes and the summary CSV so you can identify full-length
#     ERVs (LTR-GAG-POL-ENV-LTR) without cross-referencing an external database.
#     The resolved lookup is saved as cdd_name_lookup.tsv in the output directory.
#     If NCBI is unreachable, raw accessions (gnl|CDD|XXXXX) are used as fallback.
#   - Requires internet access for the CDD name lookup step only.
#   - The GFF OriginalBlastHit attribute records the exact BLAST hit coordinates (before flanking
#     extension) so you can see what the LTR query matched before the 14 kb window was added.
#
# =============================================================================

echo "Started at $(date +%T)"

# Exit immediately on error, treat unset variables as errors, propagate pipe failures
set -euo pipefail

# =============================================================================
# ARGUMENT PARSING AND VALIDATION
# =============================================================================

# Check minimum required arguments (4 required, 1 optional)
if [ "$#" -lt 4 ]; then
    echo "Usage: $0 <genome_fasta> <query_ltrs_fasta> <stockholm_file> <cdd_db_dir> [flank_size]"
    echo ""
    echo "  genome_fasta      - Genome assembly FASTA"
    echo "  query_ltrs_fasta  - RM2 LTR family consensus FASTA"
    echo "  stockholm_file    - RM2 Stockholm alignment file (.stk)"
    echo "  cdd_db_dir        - CDD database path for rpstblastn"
    echo "  flank_size        - (Optional) Flanking bp around BLAST hits (default: 14000)"
    exit 1
fi

genome_fasta="$1"       # Genome FASTA - used for BLAST database and coordinate extraction
query_fasta="$2"        # RM2 consensus FASTA - used as BLAST query sequences
stockholm_file="$3"     # RM2 Stockholm file - parsed by the Python section below
cdd_db="$4"             # CDD database directory - used by rpstblastn for domain annotation
flank_size="${5:-14000}"   # Flanking region size in bp (default 14000); passed as 5th arg if desired

# Number of parallel rpstblastn jobs to run simultaneously during CDD annotation.
# Each job already uses -num_threads 14 internally, so this multiplies CPU usage.
# Set lower (e.g. 2) if memory is constrained; higher (e.g. 16) on large servers.
CDD_PARALLEL_JOBS=8

# Validate that input files/directories actually exist before starting
if [ ! -f "$genome_fasta" ]; then
    echo "ERROR: Genome FASTA not found: $genome_fasta"
    exit 1
fi
if [ ! -f "$query_fasta" ]; then
    echo "ERROR: Query FASTA not found: $query_fasta"
    exit 1
fi
if [ ! -f "$stockholm_file" ]; then
    echo "ERROR: Stockholm file not found: $stockholm_file"
    exit 1
fi
if [ ! -d "$cdd_db" ] && [ ! -f "$cdd_db" ]; then
    echo "ERROR: CDD database not found: $cdd_db"
    exit 1
fi

echo "Flank size set to: ${flank_size} bp"

# =============================================================================
# OUTPUT DIRECTORY SETUP
# =============================================================================

# Create a unique timestamped output directory so repeated runs don't overwrite each other.
# The name is derived from the genome filename (stripping .fasta/.fa/.fna extensions).
timestamp=$(date +"%Y%m%d_%H%M%S")
base_name=$(basename "$genome_fasta")
base_name="${base_name%.fasta}"
base_name="${base_name%.fa}"
base_name="${base_name%.fna}"

out_dir="${base_name}_RM_pipe_${timestamp}"
mkdir -p "$out_dir"
echo "Output directory: $out_dir"

# Subdirectory to hold per-family FASTA files (split by family/subtype for CDD analysis)
fasta_out_dir="$out_dir/${base_name}_region_fastas"
mkdir -p "$fasta_out_dir"

# =============================================================================
# PYTHON SECTION: Stockholm parsing, LTR-INT matching, family classification
# =============================================================================
# This replaces type_matches.py and match_process.py as an embedded Python script.
# It is written to a temp file and executed, then removed.
#
# Overview of what the Python does:
#   1. parse_stockholm_blocks()   - Read the Stockholm file, extract type/ID/coordinates per block
#   2. match_ltr_int_families()   - Find LTR-INT pairs on the same chromosome within max_distance
#   3. write_matches_to_csv()     - Write matched pairs to ltr_int_matches.csv
#   4. write_all_sequences_csv()  - Write all sequences with Matched=Yes/No to all_ltr_int_sequences.csv
#   5. Family classification      - Group families into yes/no matched lists
#   6. family_match_results.csv   - Family-level source->target match map

python3 - "$stockholm_file" << 'PYTHON_SCRIPT'
import sys
import re
import csv
from collections import defaultdict

# ---- Read the Stockholm file path from command line argument ----
stockholm_file = sys.argv[1]

# -------------------------------------------------------------------------
# FUNCTION: parse_stockholm_blocks
# Reads the Stockholm file and splits it into per-family alignment blocks.
# Each block corresponds to one repeat family annotated by RM2.
#
# For each block, extracts:
#   - block_id:   the family identifier (from #=GF ID line)
#   - block_type: the repeat type (e.g. LTR or INT, from Type= in the header)
#   - sequences:  list of (chrom, start, end, strand, full_seq_id) tuples
#                 parsed from alignment rows matching chr:start-end_strand format
# -------------------------------------------------------------------------
def parse_stockholm_blocks(file_path):
    with open(file_path) as f:
        content = f.read()

    # Split file into individual blocks at each #=GF ID marker
    raw_blocks = re.split(r'(?=#=GF ID)', content)
    blocks = []

    for block in raw_blocks:
        if not block.strip():
            continue

        # Extract the repeat type (e.g. LTR, INT) from the block header annotation
        type_match = re.search(r'Type=([^,\]]+)', block)
        # Extract the family ID from the #=GF ID line
        id_match = re.search(r'#=GF ID\s+(\S+)', block)

        if not type_match or not id_match:
            continue  # Skip blocks that lack type or ID information

        block_type = type_match.group(1).strip()
        block_id = id_match.group(1).strip()

        # Parse all sequence alignment rows in this block.
        # Rows have format: chr:start-end_strand  <aligned_sequence>
        sequences = []
        for line in block.splitlines():
            if re.match(r'^[A-Za-z0-9_.]+:\d+-\d+[_+-]', line):
                seq_id = line.split()[0]
                chrom_match = re.match(r'^([^:]+):(\d+)-(\d+)_([+-])', seq_id)
                if chrom_match:
                    chrom, start, end, strand = chrom_match.groups()
                    sequences.append((chrom, int(start), int(end), strand, seq_id))

        blocks.append({
            "id": block_id,
            "type": block_type,
            "sequences": sequences
        })

    return blocks


# -------------------------------------------------------------------------
# FUNCTION: match_ltr_int_families
# Finds pairs of LTR and INT sequences that are co-located on the same chromosome.
#
# Two sequences are considered a pair if they are on the same chromosome AND:
#   - Their start coordinates are within max_distance bp of each other, OR
#   - Their end coordinates are within max_distance bp of each other, OR
#   - They physically overlap (i.e. one spans into the other)
#
# Returns:
#   matches          - list of dicts, one per LTR-INT pair with full coordinates
#   matched_ltr_ids  - set of LTR sequence IDs that were matched
#   matched_int_ids  - set of INT sequence IDs that were matched
# -------------------------------------------------------------------------
def match_ltr_int_families(blocks, max_distance=1000):
    ltr_blocks = [b for b in blocks if b["type"] == "LTR"]
    int_blocks = [b for b in blocks if b["type"] == "INT"]

    matches = []
    matched_ltr_ids = set()
    matched_int_ids = set()

    for ltr in ltr_blocks:
        for ltr_seq in ltr["sequences"]:
            ltr_chrom, ltr_start, ltr_end, ltr_strand, ltr_full_id = ltr_seq
            for int_blk in int_blocks:
                for int_seq in int_blk["sequences"]:
                    int_chrom, int_start, int_end, int_strand, int_full_id = int_seq

                    # Must be on the same chromosome
                    if ltr_chrom != int_chrom:
                        continue

                    # Check proximity: within max_distance at either end, or overlapping
                    if (
                        abs(ltr_start - int_start) < max_distance or
                        abs(ltr_end - int_end) < max_distance or
                        (ltr_start <= int_end and int_start <= ltr_end)
                    ):
                        matches.append({
                            "LTR_family": ltr["id"],
                            "LTR_seq_id": ltr_full_id,
                            "LTR_chrom": ltr_chrom,
                            "LTR_start": ltr_start,
                            "LTR_end": ltr_end,
                            "LTR_strand": ltr_strand,
                            "INT_family": int_blk["id"],
                            "INT_seq_id": int_full_id,
                            "INT_chrom": int_chrom,
                            "INT_start": int_start,
                            "INT_end": int_end,
                            "INT_strand": int_strand
                        })
                        matched_ltr_ids.add(ltr_full_id)
                        matched_int_ids.add(int_full_id)

    return matches, matched_ltr_ids, matched_int_ids


# -------------------------------------------------------------------------
# FUNCTION: write_matches_to_csv
# Writes all LTR-INT matched pairs to ltr_int_matches.csv.
# Each row contains the full coordinates for both members of the pair.
# -------------------------------------------------------------------------
def write_matches_to_csv(matches, output_file):
    fieldnames = [
        "LTR_family", "LTR_seq_id", "LTR_chrom", "LTR_start", "LTR_end", "LTR_strand",
        "INT_family", "INT_seq_id", "INT_chrom", "INT_start", "INT_end", "INT_strand"
    ]
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for match in matches:
            writer.writerow(match)


# -------------------------------------------------------------------------
# FUNCTION: write_all_sequences_csv
# Writes every sequence from every block to all_ltr_int_sequences.csv.
# Includes a Matched column (Yes/No) based on whether the sequence was part
# of a matched LTR-INT pair.
# -------------------------------------------------------------------------
def write_all_sequences_csv(blocks, matched_ltr_ids, matched_int_ids, output_file):
    fieldnames = ["Family", "Type", "Seq_ID", "Chrom", "Start", "End", "Strand", "Matched"]
    with open(output_file, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for block in blocks:
            for seq in block["sequences"]:
                chrom, start, end, strand, full_id = seq
                matched = (
                    "Yes" if (block["type"] == "LTR" and full_id in matched_ltr_ids)
                    or (block["type"] == "INT" and full_id in matched_int_ids)
                    else "No"
                )
                writer.writerow({
                    "Family": block["id"],
                    "Type": block["type"],
                    "Seq_ID": full_id,
                    "Chrom": chrom,
                    "Start": start,
                    "End": end,
                    "Strand": strand,
                    "Matched": matched
                })


# ---- Main execution ----

matched_output  = "ltr_int_matches.csv"
all_seqs_output = "all_ltr_int_sequences.csv"

# Step 1: Parse Stockholm file into structured blocks
blocks = parse_stockholm_blocks(stockholm_file)

# Step 2: Identify LTR-INT pairs within 1000 bp proximity
matches, matched_ltr_ids, matched_int_ids = match_ltr_int_families(blocks, max_distance=1000)

# Step 3: Write sequence-level outputs
write_matches_to_csv(matches, matched_output)
write_all_sequences_csv(blocks, matched_ltr_ids, matched_int_ids, all_seqs_output)

print(f"Matches written to: {matched_output}")
print(f"All sequences written to: {all_seqs_output}")

# ---- Family-level classification (from match_process.py) ----

# Read the all_sequences CSV and group by family name (first column).
# A family goes into yes_list if ANY of its sequences were matched,
# otherwise into no_list if ALL sequences are unmatched.
import csv as csv2
from collections import defaultdict

family_matched = defaultdict(list)  # family_name -> list of Matched values

with open(all_seqs_output, newline='') as f:
    reader = csv2.DictReader(f)
    for row in reader:
        family = row["Family"]
        matched_val = row["Matched"].strip().lower()
        family_matched[family].append(matched_val)

yes_list = []
no_list  = []

for family, matched_vals in family_matched.items():
    if any(v == "yes" for v in matched_vals):
        yes_list.append(family)
    else:
        no_list.append(family)

with open('matched_families_yes.txt', 'w') as f:
    for name in yes_list:
        f.write(f"{name}\n")

with open('matched_families_no.txt', 'w') as f:
    for name in no_list:
        f.write(f"{name}\n")

print(f"Matched families (yes): {len(yes_list)}")
print(f"Unmatched families (no): {len(no_list)}")

# ---- Build family-level match map from ltr_int_matches.csv ----
# Columns: LTR_family (col 0) and INT_family (col 6)
# Builds source_family -> set of matched INT families, deduplicating via sets.

matches_dict = defaultdict(set)

with open(matched_output, newline='') as f:
    reader = csv2.reader(f)
    next(reader)  # skip header
    for row in reader:
        if len(row) < 7:
            continue
        source = row[0].strip()   # LTR_family
        target = row[6].strip()   # INT_family
        matches_dict[source].add(target)

# Write family_match_results.csv: one row per unique source-target pair
match_output_rows = []
for source_family, matched_families in matches_dict.items():
    for target_family in matched_families:
        match_output_rows.append((source_family, target_family))

with open('family_match_results.csv', 'w', newline='') as f:
    writer = csv2.writer(f)
    writer.writerow(['Source_Family', 'Matched_Family'])
    writer.writerows(match_output_rows)

print("family_match_results.csv written.")
PYTHON_SCRIPT

# Move all Python outputs into the main output directory
mv matched_families_yes.txt "$out_dir/"
mv matched_families_no.txt  "$out_dir/"
mv family_match_results.csv "$out_dir/"
mv ltr_int_matches.csv      "$out_dir/"
mv all_ltr_int_sequences.csv "$out_dir/"

echo "✅ Python parsing and classification complete"

# =============================================================================
# VARIABLE SETUP FOR DOWNSTREAM STEPS
# =============================================================================

match_csv="$out_dir/family_match_results.csv"       # Family-level LTR-INT match pairs
match_no_csv="$out_dir/matched_families_no.txt"     # Families with no LTR-INT pairs
filtered_fasta="$out_dir/matched_source_families.fa" # FASTA of matched LTR families
family_blast="$out_dir/${base_name}_source_family_blast.txt" # BLAST output for matched LTRs

# Output files for unmatched families (split by type)
ltr_families_fa="$out_dir/ltr_matched_families_no.fa"   # FASTA: LTR-type unmatched families
int_families_fa="$out_dir/int_matched_families_no.fa"   # FASTA: INT-type unmatched families
blast_ltr_add="$out_dir/blast_ltr_from_no.txt"          # BLAST results: unmatched LTR families
blast_int="$out_dir/blast_int_from_no.txt"              # BLAST results: unmatched INT families

# =============================================================================
# STEP 0a: EXTRACT MATCHED SOURCE FAMILIES AND RUN BLAST
# =============================================================================
# Extract the unique source (LTR) family names from family_match_results.csv,
# then subset the query FASTA to only those families.
# Further filter to only #LTR/-annotated entries (exclude INT-type entries).
# Finally run BLASTn against the genome to find all insertion sites.

echo "Running match_process steps..."

# Pull unique Source_Family values (skip CSV header, sort for uniqueness)
cut -d',' -f1 "$match_csv" | tail -n +2 | sort -u > "$out_dir/source_families.txt"

# Filter query FASTA to only the matched source families.
# The FASTA header format is: >FamilyName#Type/Subtype_...
# We strip the #... suffix to get the bare family name for matching.
awk '
    BEGIN {
        while ((getline < "'"$out_dir/source_families.txt"'") > 0)
            keep[$1] = 1
    }
    /^>/ {
        id = substr($0, 2);
        fam = id;
        sub(/#.*/, "", fam);           # strip everything after # to get family name
        printing = (fam in keep);
    }
    printing
' "$query_fasta" > "$filtered_fasta"

# Further filter the matched families FASTA to only #LTR/-type sequences.
# This excludes any INT-domain sequences in the matched set.
ltr_only_fasta="$out_dir/matched_source_families_LTR_only.fa"
awk '/^>/ { keep = ($0 ~ /#LTR\//) } keep' "$filtered_fasta" > "$ltr_only_fasta"

# =============================================================================
# STEP 0b: BUILD BLAST DATABASE (reuse if already present)
# =============================================================================
# makeblastdb creates three index files (.nhr, .nin, .nsq).
# We check for those to avoid rebuilding unnecessarily on repeated runs.

blast_db="$out_dir/${base_name}_blastdb"

if [[ -f "${blast_db}.nhr" && -f "${blast_db}.nin" && -f "${blast_db}.nsq" ]]; then
    echo "✅ Using existing BLAST database: $blast_db"
else
    echo "Building BLAST database from $genome_fasta..."
    makeblastdb -in "$genome_fasta" -dbtype nucl -out "$blast_db"
fi

# =============================================================================
# STEP 0c: BLAST MATCHED LTR FAMILIES AGAINST GENOME
# =============================================================================
# BLASTn parameters:
#   -evalue 1e-5         : only report hits with E-value <= 1e-5
#   -qcov_hsp_perc 60    : query must cover at least 60% of the subject HSP
#   -task blastn         : standard nucleotide BLAST (not megablast)
#   -outfmt 6            : tabular format with standard columns

blastn -query "$ltr_only_fasta" -db "$blast_db" \
    -evalue 1e-5 -qcov_hsp_perc 60 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "$family_blast" -task "blastn"

echo "BLAST of matched LTR families done"

# =============================================================================
# STEP 0d: HANDLE UNMATCHED FAMILIES (from matched_families_no.txt)
# =============================================================================
# Families with no LTR-INT pair need to be split by type (LTR vs INT),
# each BLASTed separately, and tracked under their own category labels.

# Extract unique family IDs from the no-match list
cut -d',' -f1 "$match_no_csv" | tail -n +2 | sort -u > "$out_dir/families_no_ids.txt"

# Classify each unmatched family as LTR-type or INT-type using FASTA header annotations.
# Writes ltr_ids.txt and int_ids.txt into the output directory.
awk '
    BEGIN {
        while ((getline < "'"$out_dir/families_no_ids.txt"'") > 0)
            ids[$1] = 1
    }
    /^>/ {
        id = substr($0, 2);
        fam = id; sub(/#.*/, "", fam);
        printing = 0;
        if (fam in ids) {
            if ($0 ~ /Type=LTR/) ltr_ids[fam] = 1;
            else                 int_ids[fam] = 1;
        }
    }
    END {
        for (f in ltr_ids) print f > "'"$out_dir/ltr_ids.txt"'";
        for (f in int_ids) print f > "'"$out_dir/int_ids.txt"'";
    }
' "$query_fasta"

# Build the LTR-type unmatched FASTA
awk '
    BEGIN {
        while ((getline < "'"$out_dir/ltr_ids.txt"'") > 0)
            keep[$1] = 1
    }
    /^>/ {
        id = substr($0, 2);
        fam = id; sub(/#.*/, "", fam);
        printing = (fam in keep);
    }
    printing
' "$query_fasta" > "$ltr_families_fa"

# Build the INT-type unmatched FASTA
awk '
    BEGIN {
        while ((getline < "'"$out_dir/int_ids.txt"'") > 0)
            keep[$1] = 1
    }
    /^>/ {
        id = substr($0, 2);
        fam = id; sub(/#.*/, "", fam);
        printing = (fam in keep);
    }
    printing
' "$query_fasta" > "$int_families_fa"

# From both unmatched FASTAs, extract only the #LTR/-type entries for BLAST
ltr_only_fasta_ltr="$out_dir/matched_source_families_LTR_only_ltr_unmatched.fa"
ltr_only_fasta_int="$out_dir/matched_source_families_LTR_only_int_unmatched.fa"
awk '/^>/ { keep = ($0 ~ /#LTR\//) } keep' "$ltr_families_fa" > "$ltr_only_fasta_ltr"
awk '/^>/ { keep = ($0 ~ /#LTR\//) } keep' "$int_families_fa" > "$ltr_only_fasta_int"

# BLAST unmatched LTR families
blastn -query "$ltr_only_fasta_ltr" -db "$blast_db" \
    -evalue 1e-5 -qcov_hsp_perc 60 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "$blast_ltr_add" -task "blastn"
echo "BLAST of unmatched LTR families done"

# BLAST unmatched INT families
blastn -query "$ltr_only_fasta_int" -db "$blast_db" \
    -evalue 1e-5 -qcov_hsp_perc 60 \
    -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" \
    -out "$blast_int" -task "blastn"
echo "BLAST of unmatched INT families done"

# =============================================================================
# STEP 1: TAG AND COMBINE ALL BLAST RESULTS
# =============================================================================
# Append a category label (last column) to each BLAST result file so we can
# track which analysis category each hit came from after merging.

awk '{print $0 "\tltr_matched"}'   "$family_blast"  > "$out_dir/tmp_ltr_matched.txt"
awk '{print $0 "\tltr_unmatched"}' "$blast_ltr_add" > "$out_dir/tmp_ltr_unmatched.txt"
awk '{print $0 "\tint_unmatched"}' "$blast_int"     > "$out_dir/tmp_int_unmatched.txt"

combined_blast="$out_dir/all_hits_combined.txt"
cat "$out_dir"/tmp_*.txt > "$combined_blast"

# Remove temporary tagged files now that they are merged
rm "$out_dir"/tmp_*.txt

# =============================================================================
# STEP 2: DEDUPLICATE OVERLAPPING BLAST HITS
# =============================================================================
# When multiple BLAST hits overlap on the same chromosome (within a 10 bp
# tolerance), only the hit with the highest percent identity is kept.
# This prevents double-counting insertions detected by multiple query sequences.
#
# The gawk script:
#   - For each hit, compares it against all stored hits on the same chromosome
#   - If an overlap is found, keeps whichever has higher pident
#   - Replaced hits are written to a separate removed_duplicates file for audit
#   - Non-overlapping hits are stored and printed at END

dedup_blast="$out_dir/all_hits_deduplicated.txt"
removed_duplicates="$out_dir/all_hits_removed_duplicates.txt"

gawk -v tolerance=10 -v removed_out="$removed_duplicates" '
function overlap(a_start, a_end, b_start, b_end, tol) {
    return (a_end + tol >= b_start) && (b_end + tol >= a_start)
}
BEGIN { OFS = "\t" }
{
    qid    = $1
    chr    = $2
    pident = $3 + 0
    sstart = $9
    send   = $10
    # Normalise so start < end regardless of strand
    start = (sstart < send) ? sstart : send
    end   = (sstart > send) ? sstart : send

    found = 0
    for (i = 1; i <= n_hits[chr]; i++) {
        if (overlap(start, end, hit_start[chr,i], hit_end[chr,i], tolerance)) {
            found = 1
            if (pident > hit_pident[chr,i]) {
                # Current hit is better; replace the stored one
                print hit_line[chr,i], "REPLACED_BY=" $0 >> removed_out
                hit_line[chr,i]   = $0
                hit_pident[chr,i] = pident
                hit_start[chr,i]  = start
                hit_end[chr,i]    = end
            } else {
                # Stored hit is better; discard current
                print $0, "REPLACED_BY=" hit_line[chr,i] >> removed_out
            }
            break
        }
    }
    if (!found) {
        n = ++n_hits[chr]
        hit_line[chr,n]   = $0
        hit_pident[chr,n] = pident
        hit_start[chr,n]  = start
        hit_end[chr,n]    = end
    }
}
END {
    for (chr in n_hits) {
        for (i = 1; i <= n_hits[chr]; i++) {
            print hit_line[chr,i]
        }
    }
}
' "$combined_blast" > "$dedup_blast"

echo "✅ Deduplicated hits saved to $dedup_blast"
echo "✅ Step 2 complete: BLAST deduplication done - $(date +%T)"

# =============================================================================
# GENOME INDEX: required for coordinate clamping in flanking region step
# =============================================================================
# samtools faidx produces a .fai index with chromosome names and lengths.
# We use this to ensure extended coordinates don't exceed chromosome boundaries.

samtools faidx "$genome_fasta"
cut -f1,2 "${genome_fasta}.fai" > "$out_dir/genome.chrom.sizes"

# =============================================================================
# QUERY SEQUENCE LENGTHS: needed to infer missing query coverage
# =============================================================================
# For each query sequence in the LTR FASTA, compute total length.
# Multi-line FASTA sequences are summed per ID.
# This is used in run_blast_pipeline to determine left/right coverage gaps.

query_lengths="$out_dir/query_lengths.tsv"
awk '/^>/ { gsub(/^>/,""); id=$1; next } { print id, length($0) }' "$query_fasta" | \
awk '{ a[$1]+=$2 } END { for (k in a) print k, a[k] }' > "$query_lengths"

echo "✅ Step 3 complete - $(date +%T)"

# =============================================================================
# FUNCTION: run_blast_pipeline
# =============================================================================
# Runs the downstream analysis for one BLAST category (ltr_matched,
# ltr_unmatched, or int_unmatched).
#
# Steps:
#   1. Subset deduplicated BLAST hits for this category
#   2. Add strand column and filter out NW_ scaffolds; sort the result
#   3. Create flanking BED regions (size controlled by $flank_size)
#   4. Extract flanking FASTA sequences with bedtools
#   5. Run CDD domain annotation via process_fasta_group()
#
# Arguments:
#   $1  label - one of: ltr_matched, ltr_unmatched, int_unmatched
# =============================================================================
run_blast_pipeline () {
    local label="$1"
    local blast_input="$dedup_blast"

    echo "Running pipeline for category: $label"

    # Subset the combined deduplicated BLAST file to only hits from this category
    # (the category label is in the last column, appended earlier)
    local category_blast="$out_dir/${base_name}_${label}_blast.txt"
    awk -v cat="$label" '$NF==cat {print $0}' "$blast_input" > "$category_blast"
    local filtered_blast="$category_blast"

    # ---- Step 2: Add strand, filter scaffolds, extend coordinates, sort ----
    # - Strand is inferred from BLAST hit orientation (sstart < send => +)
    # - Sequences on NW_-prefixed scaffolds are skipped.
    #   NOTE: If your genome uses different scaffold naming conventions
    #   (e.g. "scaffold_", "SUPER_", "Un_"), edit the condition below.
    # - left_missing and right_missing are computed but currently not used
    #   for coordinate adjustment; they are available for future extension.
    # - Result is sorted by: query ID, chromosome, strand, start coordinate

    local blast_sorted="$out_dir/${base_name}_${label}_blast_extended_sorted.txt"

    gawk -v len_file="$query_lengths" '
    BEGIN {
        FS = OFS = "\t"
        # Load query sequence lengths into lookup table
        while ((getline < len_file) > 0) {
            qlen[$1] = $2
        }
    }
    {
        # Skip hits on NW_-prefixed scaffolds (NCBI unplaced contigs).
        # Edit this line if your genome uses a different scaffold prefix.
        if ($2 ~ /^NW/) next

        qid      = $1
        chr      = $2
        sstart   = $9
        send     = $10
        strand   = ($9 < $10) ? "+" : "-"
        qstart   = $7
        qend     = $8
        full_len = qlen[qid]

        # How much of the query is missing from each end of the BLAST hit
        left_missing  = (qstart > 1)         ? qstart - 1        : 0
        right_missing = (qend < full_len)     ? full_len - qend   : 0

        # Ensure sstart is non-negative
        if (sstart < 0) sstart = 0

        $9  = sstart
        $10 = send
        $13 = strand   # append strand as a new column

        print
    }' "$filtered_blast" | sort -k1,1 -k2,2 -k13,13 -k9,9n > "$blast_sorted"

    # ---- Step 3: Create flanking BED regions ----
    # For each BLAST hit, extend $flank_size bp in each direction.
    # Coordinates are clamped to [0, chromosome_length] using the genome .fai index.
    #
    # The BED name field encodes:
    #   queryID::chr:origStart-origEnd(strand)_extStart_extEnd
    #
    # The origStart/origEnd/strand are the RAW BLAST hit coordinates BEFORE flanking
    # was added. They are embedded in the name so that downstream steps (GFF writing)
    # can report exactly what the query sequence matched, independent of the window size.

    local flank_bed="$out_dir/${base_name}_${label}_flank_${flank_size}bp.bed"
    local flank_fasta="$out_dir/${base_name}_${label}_flank_${flank_size}bp.fa"

    gawk -v OFS="\t" -v flank="$flank_size" -v genome_index="${genome_fasta}.fai" '
    BEGIN {
        # Load chromosome sizes from the genome .fai index
        while ((getline < genome_index) > 0) {
            chrom_sizes[$1] = $2;
        }
    }
    {
        chr    = $2
        sstart = $9
        send   = $10
        start  = (sstart < send) ? sstart : send
        end    = (sstart > send) ? sstart : send
        strand = (sstart < send) ? "+" : "-"

        ext_start = start - flank
        ext_end   = end   + flank

        # Clamp to valid chromosome coordinates
        if (ext_start < 0)              ext_start = 0
        if (ext_end > chrom_sizes[chr]) ext_end   = chrom_sizes[chr]

        # Encode original BLAST hit coords into the name so the GFF can report them.
        # Format: queryID::chr:origStart-origEnd(strand)_extStart_extEnd
        name = $1 "::" chr ":" start "-" end "(" strand ")_" ext_start "_" ext_end

        print chr, ext_start, ext_end, name, ".", strand
    }
    ' "$filtered_blast" > "$flank_bed"

    # ---- Step 4: Extract FASTA for each flanking region ----
    bedtools getfasta -fi "$genome_fasta" -bed "$flank_bed" -name -fo "$flank_fasta"

    # ---- Step 5: CDD domain annotation ----
    # Pass the blast_sorted file so the GFF writer can annotate original hit coordinates
    process_fasta_group "$flank_fasta" "$label" "$blast_sorted"

    echo "✅ Pipeline for $label complete"
}

# =============================================================================
# FUNCTION: process_fasta_group
# =============================================================================
# Splits a multi-FASTA file of flanking regions into per-family FASTA files,
# then runs rpstblastn (CDD) on each to identify conserved protein domains
# (e.g. reverse transcriptase, integrase, gag, env).
#
#   1. PARALLEL EXECUTION: rpstblastn calls are dispatched simultaneously using
#      xargs -P $CDD_PARALLEL_JOBS rather than sequentially. On a genome with
#      many LTR families this can cut CDD annotation time by 4-8x.
#
#   2. HUMAN-READABLE DOMAIN NAMES: After all rpstblastn jobs finish, every
#      unique CDD accession across all results is collected and looked up in a
#      single batch query to the NCBI CDD esummary API. This returns the short
#      name (e.g. RVT_1, Gag_p24, rve, ENV) for each accession and writes them
#      to cdd_name_lookup.tsv. GFF files are then annotated in a second pass.
#      No .smp files or cddid.tbl are required.
#
#   3. ORIGINAL BLAST HIT IN GFF: The FASTA header encodes the raw BLAST hit
#      coordinates before flanking extension, written as OriginalBlastHit in GFF
#      attributes so you can see exactly what the LTR query matched.
#
# NOTE: This function only runs the rpstblastn step and splits the FASTA.
#       GFF/CSV writing happens in a separate pass after build_cdd_name_lookup()
#       has been called once all three categories have run their rpstblastn jobs.
#
# Arguments:
#   $1  input_fasta  - FASTA file of flanking regions
#   $2  label        - category label (ltr_matched / ltr_unmatched / int_unmatched)
#   $3  blast_sorted - the per-category blast_extended_sorted.txt (coords in headers)
# =============================================================================
process_fasta_group () {
    local input_fasta="$1"
    local label="$2"
    local blast_sorted="$3"
    local out_dir_base="$fasta_out_dir/$label"
    mkdir -p "$out_dir_base"

    # -------------------------------------------------------------------------
    # SPLIT FLANKING FASTA INTO PER-FAMILY FILES
    # -------------------------------------------------------------------------
    # FASTA headers have the format set in run_blast_pipeline:
    #   >queryID::chr:origStart-origEnd(strand)_extStart_extEnd
    # (bedtools getfasta appends ::extStart-extEnd after -name)
    # We parse out family name and type/subtype for the filename.

    gawk -v out_dir="$out_dir_base" '
    BEGIN { filename = "" }
    /^>/ {
        raw_header   = substr($0, 2);
        clean_header = raw_header;
        gsub(/[^A-Za-z0-9_.:()\-]/, "_", clean_header);

        # Family name is the part before the first "::"
        split(raw_header, nameparts, "::");
        query_part = nameparts[1];

        split(query_part, hashparts, "#");
        fam        = hashparts[1];
        after_hash = hashparts[2];
        split(after_hash, typeparts, "_");
        split(typeparts[1], ts, "/");
        type    = ts[1];
        subtype = ts[2];

        filetag = (type == "" || subtype == "") ? fam "_Unknown" : fam "_" type "_" subtype;
        gsub(/[^A-Za-z0-9_.-]/, "_", filetag);

        filename = out_dir "/" filetag ".fa";
        print ">" clean_header > filename;
        next;
    }
    { print > filename; }
    ' "$input_fasta"

    # -------------------------------------------------------------------------
    # PARALLEL rpstblastn
    # -------------------------------------------------------------------------
    # Run CDD search on each per-family FASTA concurrently.
    # Results are written to <family>.fa_cdd_results.txt beside the FASTA.
    # GFF and CSV writing happens later once the name lookup is available.
    # -num_threads 4 per job; total threads = CDD_PARALLEL_JOBS x 4.

    export cdd_db CDD_PARALLEL_JOBS

    fa_count=$(find "$out_dir_base" -maxdepth 1 -name "*.fa" | wc -l)
    if [ "$fa_count" -eq 0 ]; then
        echo "  WARNING: No per-family FASTA files found in $out_dir_base — skipping rpstblastn for $label"
        return 0
    fi

    echo "  Running rpstblastn on $fa_count family FASTAs (parallel jobs: $CDD_PARALLEL_JOBS)..."

    # Derive the RPS database stem for rpstblastn.
    #
    # Single-volume databases have one Cdd.rps file -> stem is Cdd/Cdd
    # Multi-volume databases (split into Cdd.00.rps, Cdd.01.rps etc.) have a
    # Cdd.pal alias file that lists all volumes. rpstblastn takes the stem of
    # the .pal file and searches all volumes automatically via the alias.
    # We prefer the .pal file if present; fall back to a single .rps stem.

    pal_file=$(find "$cdd_db" -maxdepth 1 -name "*.pal" | head -1)
    if [ -n "$pal_file" ]; then
        export cdd_db_stem="${pal_file%.pal}"
        echo "  Multi-volume RPS database detected — using alias: $cdd_db_stem"
    else
        rps_file=$(find "$cdd_db" -maxdepth 1 -name "*.rps" | head -1)
        if [ -z "$rps_file" ]; then
            echo "  ERROR: No .pal or .rps file found in $cdd_db"
            echo "  Check that \$cdd_db points to the directory containing the CDD files"
            return 1
        fi
        export cdd_db_stem="${rps_file%.rps}"
        echo "  Single-volume RPS database detected — using: $cdd_db_stem"
    fi

    # Disable set -e around the xargs block so a single failing rpstblastn job
    # does not silently kill the whole pipeline. Each job logs its own success
    # or failure explicitly. set -e is restored immediately after.
    set +e
    find "$out_dir_base" -maxdepth 1 -name "*.fa" | \
    xargs -P "$CDD_PARALLEL_JOBS" -I{} bash -c '
        fasta="$1"
        fname=$(basename "$fasta" .fa)
        cdd_output="${fasta%.fa}_cdd_results.txt"
        stderr_log="${fasta%.fa}_cdd_stderr.txt"

        # Skip empty FASTA files rather than letting rpstblastn error on them
        if [ ! -s "$fasta" ]; then
            echo "  SKIP (empty FASTA): $fname"
            exit 0
        fi

        echo "  Starting rpstblastn: $fname"

        rpstblastn -query "$fasta" -db "$cdd_db_stem" \
            -evalue 1e-3 \
            -outfmt "6 qseqid sseqid evalue qstart qend" \
            -num_threads 4 \
            -max_target_seqs 40 \
            -out "$cdd_output" 2>"$stderr_log"
        exit_code=$?

        if [ $exit_code -ne 0 ]; then
            echo "  WARNING: rpstblastn failed (exit $exit_code) for: $fname"
            echo "    stderr: $(cat $stderr_log)"
        else
            hit_count=$(wc -l < "$cdd_output" 2>/dev/null || echo 0)
            echo "  Done: $fname ($hit_count hits)"
            # Remove empty stderr log to keep directory tidy
            [ ! -s "$stderr_log" ] && rm -f "$stderr_log"
        fi
    ' -- {}
    xargs_exit=$?
    set -e   # restore strict error handling

    if [ $xargs_exit -ne 0 ]; then
        echo "  WARNING: one or more rpstblastn jobs returned non-zero for $label (exit $xargs_exit)"
        echo "  Check *_cdd_stderr.txt files in $out_dir_base for details"
        echo "  Pipeline continuing — partial results may exist"
    fi

    echo "  ✅ rpstblastn complete for: $label"
}

# =============================================================================
# FUNCTION: build_cdd_name_lookup
# =============================================================================
# Called ONCE after all three categories have finished their rpstblastn runs.
#
# Collects every unique CDD accession from all *_cdd_results.txt files across
# all categories, then queries the NCBI CDD esummary API in batches of 200
# accessions to retrieve human-readable short names (e.g. RVT_1, Gag_p24,
# rve, ENV_gp41).
#
# The accessions returned by rpstblastn have the form: gnl|CDD|XXXXX
# where XXXXX is the numeric PSSM ID. The esummary API accepts these IDs
# directly via the 'cdd' database.
#
# Output: $out_dir/cdd_name_lookup.tsv  (two columns: accession <TAB> short_name)
#
# If the network is unavailable or NCBI returns an error, the lookup file will
# be empty/partial and raw accessions will be used as a fallback in write_gff_pass.
#
# The lookup is cached in the output directory — if you re-run the GFF pass
# without re-running rpstblastn you can skip this step and reuse the file.
# =============================================================================
build_cdd_name_lookup () {
    local lookup_out="$out_dir/cdd_name_lookup.tsv"
    cdd_name_lookup="$lookup_out"
    export cdd_name_lookup

    echo "Collecting CDD accessions from rpstblastn results..."

    # Gather every unique accession from column 2 of all cdd_results.txt files.
    # Accession format from rpstblastn: gnl|CDD|XXXXX  where XXXXX is numeric PSSM ID.
    find "$fasta_out_dir" -name "*_cdd_results.txt" -exec cut -f2 {} \; \
        | sort -u > "$out_dir/all_cdd_accessions.txt"

    local n_acc
    n_acc=$(wc -l < "$out_dir/all_cdd_accessions.txt")
    echo "  Found $n_acc unique CDD accessions — querying NCBI esummary API..."

    # Query NCBI in batches and write the TSV lookup
    python3 - "$out_dir/all_cdd_accessions.txt" "$lookup_out" << 'NCBI_LOOKUP'
import sys
import re
import time
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET

accessions_file = sys.argv[1]
output_tsv      = sys.argv[2]

# Read all accessions; extract numeric PSSM ID from "gnl|CDD|XXXXX" format.
# If the accession is already numeric or in another format, use it directly.
acc_list = []
raw_to_id = {}   # maps original accession string -> numeric ID for the API
with open(accessions_file) as f:
    for line in f:
        acc = line.strip()
        if not acc:
            continue
        m = re.search(r'gnl\|CDD\|(\d+)', acc)
        if m:
            num_id = m.group(1)
        elif re.match(r'^\d+$', acc):
            num_id = acc
        else:
            # Accession string we can't resolve to a numeric ID; skip API,
            # store as-is so the fallback uses the raw string as name.
            raw_to_id[acc] = None
            acc_list.append(acc)
            continue
        raw_to_id[acc] = num_id
        acc_list.append(acc)

# Batch the numeric IDs for esummary (max 200 per request, with polite delay)
BATCH_SIZE = 200
ESUMMARY   = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
name_map   = {}   # numeric_id -> short_name

numeric_ids = [raw_to_id[a] for a in acc_list if raw_to_id.get(a) is not None]
numeric_ids = list(dict.fromkeys(numeric_ids))  # deduplicate, preserve order

for i in range(0, len(numeric_ids), BATCH_SIZE):
    batch = numeric_ids[i : i + BATCH_SIZE]
    params = urllib.parse.urlencode({
        'db':      'cdd',
        'id':      ','.join(batch),
        'retmode': 'xml',
    })
    url = f"{ESUMMARY}?{params}"
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            xml_data = resp.read()
        root = ET.fromstring(xml_data)
        # Each <DocSum> contains:
        #   <Id>          - numeric PSSM ID
        #   <Item Name="Accession">  - CDD accession e.g. pfam00429
        #   <Item Name="Title">      - short name e.g. TLV_coat, RVT_1, Gag_p24
        #   <Item Name="Subtitle">   - longer description e.g. "ENV polyprotein (coat polyprotein)"
        #
        # We use Title as the primary readable name (it contains the short domain
        # name directly, e.g. TLV_coat). Subtitle is stored as a bonus comment
        # in the lookup TSV so you can see the description if needed.
        # ShortName does not appear in this API version so we do not look for it.
        for docsum in root.findall('DocSum'):
            id_el = docsum.find('Id')
            if id_el is None:
                continue
            pssm_id = id_el.text.strip()
            title    = None
            subtitle = None
            for item in docsum.findall('Item'):
                iname = item.get('Name', '')
                if iname == 'Title' and item.text:
                    title = item.text.strip()
                elif iname == 'Subtitle' and item.text:
                    subtitle = item.text.strip()
            if title:
                # Store "ShortName\tSubtitle" so the lookup TSV has both;
                # the gawk reader only uses column 2 (title) but the full
                # subtitle is there if you want to inspect the TSV manually.
                name_map[pssm_id] = (title, subtitle or "")
    except Exception as e:
        print(f"  WARNING: NCBI API batch {i//BATCH_SIZE + 1} failed: {e}")
        print("  Remaining accessions in this batch will use raw IDs as fallback.")
    time.sleep(0.34)   # NCBI rate limit: max 3 requests/second without API key

# Write lookup TSV: original_accession <TAB> title <TAB> subtitle
# The gawk reader uses column 2 (title) as the display name.
# Column 3 (subtitle) is for human inspection of cdd_name_lookup.tsv only.
written = 0
with open(output_tsv, 'w') as out:
    for acc in acc_list:
        num_id = raw_to_id.get(acc)
        if num_id and num_id in name_map:
            title, subtitle = name_map[num_id]
        else:
            title    = acc   # fallback: raw accession as its own name
            subtitle = ""
        out.write(f"{acc}\t{title}\t{subtitle}\n")
        written += 1

resolved = sum(1 for a in acc_list if raw_to_id.get(a) and raw_to_id[a] in name_map)
print(f"  CDD name lookup written: {written} entries ({resolved} resolved via NCBI API)")
if len(name_map) == 0:
    print("  NOTE: No names resolved. Check internet connectivity.")
    print("  Raw accessions (gnl|CDD|XXXXX) will appear in GFF and CSV outputs.")
NCBI_LOOKUP

    echo "✅ CDD name lookup complete"
}

# =============================================================================
# FUNCTION: write_gff_pass
# =============================================================================
# Second pass over all *_cdd_results.txt files to write GFF and CSV outputs
# now that the name lookup is available.
#
# Called once per category AFTER build_cdd_name_lookup() has run.
# This separation means rpstblastn (slow) and name lookup (one API call) are
# fully decoupled — the lookup is done once across all categories together.
#
# Arguments:
#   $1  label - category label (ltr_matched / ltr_unmatched / int_unmatched)
# =============================================================================
write_gff_pass () {
    local label="$1"
    local out_dir_base="$fasta_out_dir/$label"
    local summary_csv_out="$out_dir/${base_name}_${label}_cdd_summary.csv"

    echo "Family,Region,DomainName,Accession,Evalue,QueryStart,QueryEnd" > "$summary_csv_out"

    local result_count
    result_count=$(find "$out_dir_base" -maxdepth 1 -name "*_cdd_results.txt" | wc -l)
    if [ "$result_count" -eq 0 ]; then
        echo "  WARNING: No CDD result files found for $label — skipping GFF/CSV writing"
        return 0
    fi

    for cdd_output in "$out_dir_base"/*_cdd_results.txt; do
        [ -s "$cdd_output" ] || continue   # skip empty result files

        local fasta="${cdd_output%_cdd_results.txt}.fa"
        local fname
        fname=$(basename "$fasta" .fa)
        local gff_output="${cdd_output%_cdd_results.txt}_cdd_hits.gff"

        # ---- Write GFF3 and append to summary CSV ----
        # Reads the name lookup once per file (small TSV, fast).
        # GFF attributes include:
        #   Name=             human-readable domain name (e.g. RVT_1, Gag_p24)
        #   DomainName=       same as Name, for explicit querying in Geneious
        #   Accession=        raw CDD accession (e.g. gnl|CDD|47234)
        #   Evalue=           rpstblastn e-value
        #   OriginalBlastHit= the LTR BLAST hit coords before flanking extension
        #                     format: chr:start-end(strand)
        #                     parsed from the encoded FASTA header
        gawk -v fam="$fname" \
             -v name_lookup="$cdd_name_lookup" \
             -v summary_out="$summary_csv_out" \
             -v OFS="\t" \
        '
        BEGIN {
            # Load lookup TSV: col1=accession, col2=Title (short name), col3=Subtitle (description)
            # e.g.: gnl|CDD|306850  TLV_coat  ENV polyprotein (coat polyprotein)
            while ((getline line < name_lookup) > 0) {
                n = split(line, f, "\t");
                if (n >= 2) {
                    name_map[f[1]]     = f[2];           # title -> display name
                    subtitle_map[f[1]] = (n >= 3) ? f[3] : "";  # subtitle -> description
                }
            }
        }
        {
            seq_id = $1;
            acc    = $2;
            evalue = $3;
            qstart = $4;
            qend   = $5;

            # Resolve human-readable name; fall back to raw accession if absent
            dom_name = (acc in name_map) ? name_map[acc] : acc;
            desc     = (acc in subtitle_map) ? subtitle_map[acc] : "";

            # ----------------------------------------------------------------
            # Parse original BLAST hit coords from the encoded sequence ID.
            # Header format (from run_blast_pipeline BED name field):
            #   queryFam::chr:origStart-origEnd(strand)_extStart_extEnd
            # bedtools getfasta appends ::extStart-extEnd after -name, giving:
            #   queryFam::chr:origStart-origEnd(strand)_extStart_extEnd::extStart-extEnd
            # ----------------------------------------------------------------
            orig_hit = "unknown";
            if (match(seq_id, /::([^:]+):([0-9]+)-([0-9]+)\(([+-])\)/, arr)) {
                orig_hit = arr[1] ":" arr[2] "-" arr[3] "(" arr[4] ")";
            }

            # Strand from CDD hit direction; normalise coords so start < end
            strand = (qstart <= qend) ? "+" : "-";
            s = (qstart <= qend) ? qstart : qend;
            e = (qstart <= qend) ? qend   : qstart;

            # Sanitise all fields for GFF attribute string
            safe_acc  = acc;      gsub(/[\|\/# \t;=,]/, "_", safe_acc);
            safe_name = dom_name; gsub(/[\|\/# \t;=,]/, "_", safe_name);
            safe_seq  = seq_id;   gsub(/[\t;=,]/,       "_", safe_seq);
            safe_desc = desc;     gsub(/[\t;=,]/,        " ", safe_desc);

            attributes = "ID="              safe_acc "_" NR  \
                       ";Name="             safe_name         \
                       ";DomainName="       safe_name         \
                       ";Description="      safe_desc         \
                       ";Accession="        safe_acc          \
                       ";Evalue="           evalue            \
                       ";OriginalBlastHit=" orig_hit          \
                       ";Region="           safe_seq;

            print safe_seq, "CDD", safe_name, s, e, evalue, strand, ".", attributes;

            # Append row to summary CSV
            print fam "," safe_seq "," safe_name "," safe_acc "," evalue "," qstart "," qend \
                >> summary_out;
        }
        ' "$cdd_output" > "$gff_output"
    done

    echo "  ✅ GFF/CSV annotation written for: $label"
}


# =============================================================================
# STEP 3: RUN THE FULL PIPELINE FOR EACH CATEGORY
# =============================================================================
# Two-pass approach:
#
#   Pass 1 (process_fasta_group): Run all rpstblastn jobs for all three
#     categories. This is the slow step. Running all categories before the
#     API call means we collect every accession in one go.
#
#   build_cdd_name_lookup: Collect all unique CDD accessions from every
#     *_cdd_results.txt file, then make a single batched query to the NCBI
#     CDD esummary API to resolve them to human-readable short names
#     (RVT_1, Gag_p24, rve, ENV etc.). Writes cdd_name_lookup.tsv once.
#
#   Pass 2 (write_gff_pass): Write GFF3 and summary CSV for each category
#     now that the name lookup is available.

echo "=== Pass 1: rpstblastn for all categories ==="
run_blast_pipeline "ltr_matched"
run_blast_pipeline "ltr_unmatched"
run_blast_pipeline "int_unmatched"

echo "=== Building CDD name lookup via NCBI API ==="
build_cdd_name_lookup

echo "=== Pass 2: writing GFF and CSV outputs ==="
write_gff_pass "ltr_matched"
write_gff_pass "ltr_unmatched"
write_gff_pass "int_unmatched"

# =============================================================================
# STEP 4: GENERATE PER-CATEGORY SUMMARY TABLES
# =============================================================================
# For each category, combine the BED coordinate information with BLAST alignment
# metrics to produce a human-readable TSV summary table.
#
# Columns: Index, Family, Query_ID, Chromosome, Coordinates, Orientation,
#          Percent_identity, Evalue, Bitscore
#
# The BED file provides: chromosome, extended start-end, and strand.
# The BLAST file provides: percent identity, e-value, and bitscore.
# They are joined on the query ID (column 1 of BLAST / parsed from column 4 of BED).

for category in ltr_matched ltr_unmatched int_unmatched; do
    bed_file="$out_dir/${base_name}_${category}_flank_${flank_size}bp.bed"
    blast_file="$out_dir/${base_name}_${category}_blast_extended_sorted.txt"
    summary_table="$out_dir/${base_name}_${category}_family_sequences_summary.tsv"

    echo -e "Index\tFamily\tQuery_ID\tChromosome\tCoordinates\tOrientation\tPercent_identity\tEvalue\tBitscore" > "$summary_table"

    awk -F"\t" '
    NR==FNR {
        # First file (BED): build lookup of query name -> chromosome, coord range, strand
        split($4, id_parts, "#")
        fam    = id_parts[1]
        coords = $2 "-" $3
        bed_coords[$4] = $1 "\t" coords "\t" $6
        next
    }
    {
        # Second file (BLAST extended): join with BED info
        qid = $1
        fam = qid
        sub(/#.*/, "", fam)

        count[fam]++
        idx = count[fam]   # sequential index within this family

        chrom_coord_orient = bed_coords[qid]
        if (chrom_coord_orient == "") {
            chrom = "NA"; coords = "NA"; orient = "NA"
        } else {
            split(chrom_coord_orient, fields, "\t")
            chrom  = fields[1]
            coords = fields[2]
            orient = fields[3]
        }

        pident   = $3
        evalue   = $12
        bitscore = $13

        print idx, fam, qid, chrom, coords, orient, pident, evalue, bitscore
    }' OFS="\t" "$bed_file" "$blast_file" >> "$summary_table"
done

# =============================================================================
# DONE
# =============================================================================
echo ""
echo "✅ Pipeline complete at $(date +%T)"
echo "• All output saved in: $out_dir/"
echo ""
echo "Key output files:"
echo "  $out_dir/matched_families_yes.txt"
echo "  $out_dir/matched_families_no.txt"
echo "  $out_dir/family_match_results.csv"
echo "  $out_dir/all_hits_deduplicated.txt"
for category in ltr_matched ltr_unmatched int_unmatched; do
    echo "  $out_dir/${base_name}_${category}_family_sequences_summary.tsv"
    echo "  $out_dir/${base_name}_${category}_cdd_summary.csv"
    echo "    (columns: Family, Region, DomainName, Accession, Evalue, QueryStart, QueryEnd)"
done
echo ""
echo "GFF files (per family, in $fasta_out_dir/):"
echo "  Attributes include:"
echo "    Name=            human-readable domain name (e.g. RVT_1, Gag_p24, rve, ENV)"
echo "    OriginalBlastHit the LTR BLAST hit before flanking extension (chr:start-end(strand))"
echo "    Evalue           rpstblastn e-value for the domain hit"
