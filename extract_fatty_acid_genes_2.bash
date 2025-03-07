#!/bin/bash
# Define input and output files
GFF_FILE="data/intestine_prokka/intestine_annotation.gff"
FFN_FILE="data/intestine_prokka/intestine_annotation.ffn"
OUTPUT_GENES="intestine_carbon_abR_genes.ffn"
# List of genes
GENE_LIST=("xylA" "bglx" "deoB" "araD" "tetA" "lmrB" "poxtA" "efrB" "pmrA" "bcrA" "arlA" "vanRF" "mdtG")
# Extract matching gene IDs from the GFF file based on 'locus_tag=' field
echo "Extracting gene IDs..."
grep -E "$(IFS=\|; echo "${GENE_LIST[*]}")" "$GFF_FILE" | grep -o "locus_tag=[^;]*" | sed 's/locus_tag=//' > gene_ids.txt
# If no locus tags
if [ ! -s gene_ids.txt ]; then
    echo "No matching locus tags found. Check Prokka annotation for alternative gene names."
fi
# Extract corresponding protein sequences from the FFN file
if [ -s gene_ids.txt ]; then
    echo "Extracting protein sequences..."
    seqkit grep -f gene_ids.txt "$FFN_FILE" > "$OUTPUT_GENES"
    echo "Extraction complete. Genes saved to $OUTPUT_GENES"
else
echo "No matching gene IDs found. Check Prokka annotation for alternative gene names."
fi
# Define input and output files
GFF_FILE="data/burongisda_prokka/burongisda.gff"
FFN_FILE="data/burongisda_prokka/burongisda.ffn"
OUTPUT_GENES="burongisda_carbon_abR_genes.ffn"
# List of genes
GENE_LIST=("xylA" "bglx" "deoB" "araD" "tetA" "lmrB" "poxtA" "efrB" "pmrA" "bcrA" "arlA" "vanRF" "mdtG")
# Extract matching gene IDs from the GFF file based on 'locus_tag=' field
echo "Extracting gene IDs..."
grep -E "$(IFS=\|; echo "${GENE_LIST[*]}")" "$GFF_FILE" | grep -o "locus_tag=[^;]*" | sed 's/locus_tag=//' > gene_ids.txt

# If no locus tags
if [ ! -s gene_ids.txt ]; then
    echo "No matching locus tags found. Check Prokka annotation for alternative gene names."
fi
# Extract corresponding protein sequences from the FFN file
if [ -s gene_ids.txt ]; then
    echo "Extracting protein sequences..."
    seqkit grep -f gene_ids.txt "$FFN_FILE" > "$OUTPUT_GENES"
    echo "Extraction complete. genes saved to $OUTPUT_GENES"
else
echo "No matching gene IDs found. Check Prokka annotation for alternative gene names."
fi
# Define input and output files
GFF_FILE="data/strainA_prokka/strainA.gff"
FFN_FILE="data/strainA_prokka/strainA.ffn"
OUTPUT_GENES="strainA_carbon_abR_genes.ffn"
# List of genes
GENE_LIST=("xylA" "bglx" "deoB" "araD" "tetA" "lmrB" "poxtA" "efrB" "pmrA" "bcrA" "arlA" "vanRF" "mdtG")
# Extract matching gene IDs from the GFF file based on 'locus_tag=' field
echo "Extracting gene IDs..."
grep -E "$(IFS=\|; echo "${GENE_LIST[*]}")" "$GFF_FILE" | grep -o "locus_tag=[^;]*" | sed 's/locus_tag=//' > gene_ids.txt
# If no locus tags
if [ ! -s gene_ids.txt ]; then
    echo "No matching locus tags found. Check Prokka annotation for alternative gene names."
fi
# Extract corresponding protein sequences from the FFN file
if [ -s gene_ids.txt ]; then
    echo "Extracting protein sequences..."
    seqkit grep -f gene_ids.txt "$FFN_FILE" > "$OUTPUT_GENES"
    echo "Extraction complete. genes saved to $OUTPUT_GENES"
else
echo "No matching gene IDs found. Check Prokka annotation for alternative gene names."
fi
# Define input and output files
GFF_FILE="data/strainB_prokka/strainB.gff"
FFN_FILE="data/strainB_prokka/strainB.ffn"
OUTPUT_GENES="strainB_carbon_abR_genes.ffn"
# List of genes
GENE_LIST=("xylA" "bglx" "deoB" "araD" "tetA" "lmrB" "poxtA" "efrB" "pmrA" "bcrA" "arlA" "vanRF" "mdtG")
# Extract matching gene IDs from the GFF file based on 'locus_tag=' field
echo "Extracting gene IDs..."
grep -E "$(IFS=\|; echo "${GENE_LIST[*]}")" "$GFF_FILE" | grep -o "locus_tag=[^;]*" | sed 's/locus_tag=//' > gene_ids.txt
# If no locus tags
if [ ! -s gene_ids.txt ]; then
    echo "No matching locus tags found. Check Prokka annotation for alternative gene names."
fi
# Extract corresponding protein sequences from the FFN file
if [ -s gene_ids.txt ]; then
    echo "Extracting protein sequences..."
    seqkit grep -f gene_ids.txt "$FFN_FILE" > "$OUTPUT_GENES"
    echo "Extraction complete. genes saved to $OUTPUT_GENES"
else
echo "No matching gene IDs found. Check Prokka annotation for alternative gene names."
fi
# Define input and output files
GFF_FILE="data/strainC_prokka/strainC.gff"
FFN_FILE="data/strainC_prokka/strainC.ffn"
OUTPUT_GENES="strainC_carbon_abR_genes.ffn"
# List of genes
GENE_LIST=("xylA" "bglx" "deoB" "araD" "tetA" "lmrB" "poxtA" "efrB" "pmrA" "bcrA" "arlA" "vanRF" "mdtG")
# Extract matching gene IDs from the GFF file based on 'locus_tag=' field
echo "Extracting gene IDs..."
grep -E "$(IFS=\|; echo "${GENE_LIST[*]}")" "$GFF_FILE" | grep -o "locus_tag=[^;]*" | sed 's/locus_tag=//' > gene_ids.txt
# If no locus tags
if [ ! -s gene_ids.txt ]; then
    echo "No matching locus tags found. Check Prokka annotation for alternative gene names."
fi
# Extract corresponding protein sequences from the FFN file
if [ -s gene_ids.txt ]; then
    echo "Extracting protein sequences..."
    seqkit grep -f gene_ids.txt "$FFN_FILE" > "$OUTPUT_GENES"
    echo "Extraction complete. genes saved to $OUTPUT_GENES"
else
echo "No matching gene IDs found. Check Prokka annotation for alternative gene names."
fi
