#!/bin/bash

# Script to index a reference genome

# Default values
GENOME=""
TOOL=""

# Usage function
usage() {
    echo "Usage: $0 -g <path_to_genome> -t <tool_name>"
    echo "  -g  Path to the reference genome file (required)"
    echo "  -t  Tool for indexing (e.g., bwa mem, bowtie2) (required)"
    exit 1
}

# Parse command-line arguments
while getopts ":g:t:" opt; do
    case "${opt}" in
        g)
            GENOME=${OPTARG}
            ;;
        t)
            TOOL=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

# Check if required arguments are provided
if [ -z "${GENOME}" ] || [ -z "${TOOL}" ]; then
    usage
fi

# Function to index with bwa
index_bwa() {
    echo "Indexing genome with bwa..."
    bwa mem index ${GENOME}
    echo "Indexing completed."
}

# Function to index with bowtie2
index_bowtie2() {
    echo "Indexing genome with bowtie2..."
    bowtie2-build ${GENOME} ${GENOME%.fasta}
    echo "Indexing completed."
}

# Perform indexing based on the tool specified
case "${TOOL}" in
    bwa)
        index_bwa
        ;;
    bowtie2)
        index_bowtie2
        ;;
    *)
        echo "Unsupported tool: ${TOOL}"
        exit 1
        ;;
esac
