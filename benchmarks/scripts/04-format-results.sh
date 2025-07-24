#!/bin/bash

set -e  # Exit on error

# Create results directory if it doesn't exist
mkdir -p ../results

# Function to convert CSV to markdown table
csv_to_markdown() {
    local csv_file=$1
    local output_file=$2
    
    # Get the base name without extension for the title
    local title=$(basename "$csv_file" .csv)
    
    # Create markdown table header
    echo "# $title" > "$output_file"
    echo "" >> "$output_file"
    
    # Convert CSV to markdown table
    awk -F',' '
    BEGIN {
        printf "| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |\n"
        printf "|:---|---:|---:|---:|---:|---:|\n"
    }
    NR == 2 {
        # First data row, initialize min_mean
        cmd = $1
        if (cmd ~ /merkurio/) cmd = "merkurio"
        else if (cmd ~ /st-/) cmd = "seqtool"
        else if (cmd ~ /fetch/) cmd = "fetch_reads"
        else if (cmd ~ /seqkit/) cmd = "seqkit"
        else if (cmd ~ /grep/) cmd = "grep"
        else if (cmd ~ /cookiecutter/) cmd = "cookiecutter"
        else if (cmd ~ /back_to_sequences/) cmd = "back_to_sequences"
        tool[0] = cmd
        mean[0] = $2
        stddev[0] = $3
        min[0] = $7
        max[0] = $8
        min_mean = $2
        n = 1
    }
    NR > 2 {
        cmd = $1
        if (cmd ~ /merkurio/) cmd = "merkurio"
        else if (cmd ~ /st-/) cmd = "seqtool"
        else if (cmd ~ /fetch/) cmd = "fetch_reads"
        else if (cmd ~ /seqkit/) cmd = "seqkit"
        else if (cmd ~ /grep/) cmd = "grep"
        else if (cmd ~ /cookiecutter/) cmd = "cookiecutter"
        else if (cmd ~ /back_to_sequences/) cmd = "back_to_sequences"
        tool[n] = cmd
        mean[n] = $2
        stddev[n] = $3
        min[n] = $7
        max[n] = $8
        if ($2 < min_mean) min_mean = $2
        n++
    }
    END {
        for (i = 0; i < n; i++) {
            printf "| %s | %.3f | %.3f | %.3f | %.3f | %.2fx |\n", tool[i], mean[i], stddev[i], min[i], max[i], mean[i]/min_mean
        }
    }
    ' "$csv_file" >> "$output_file"
    
    echo "" >> "$output_file"
}

# Create summary markdown file
echo "# Benchmark Results Summary" > ../results/summary.md
echo "" >> ../results/summary.md
echo "## Contents" >> ../results/summary.md
echo "" >> ../results/summary.md

# Process all CSV files in the results directory
for csv_file in ../results/*/*.csv; do
    if [ -f "$csv_file" ]; then
        # Get the benchmark name and category
        benchmark_name=$(basename "$csv_file" .csv)
        category=$(basename $(dirname "$csv_file"))
        
        # Add to contents
        echo "- [$benchmark_name](#$category-$benchmark_name)" >> ../results/summary.md
    fi
done

echo "" >> ../results/summary.md

# Add all benchmark results to summary
for csv_file in ../results/*/*.csv; do
    if [ -f "$csv_file" ]; then
        benchmark_name=$(basename "$csv_file" .csv)
        category=$(basename $(dirname "$csv_file"))
        
        # Add section header
        echo "## $category: $benchmark_name" >> ../results/summary.md
        echo "" >> ../results/summary.md
        
        # Add the table
        awk -F',' '
        BEGIN {
            printf "| Tool | Mean [s] | Stddev [s] | Min [s] | Max [s] | Relative |\n"
            printf "|:---|---:|---:|---:|---:|---:|\n"
        }
        NR == 2 {
            cmd = $1
            if (cmd ~ /merkurio/) cmd = "merkurio"
            else if (cmd ~ /st-/) cmd = "seqtool"
            else if (cmd ~ /fetch/) cmd = "fetch_reads"
            else if (cmd ~ /seqkit/) cmd = "seqkit"
            else if (cmd ~ /grep/) cmd = "grep"
            else if (cmd ~ /cookiecutter/) cmd = "cookiecutter"
            else if (cmd ~ /back_to_sequences/) cmd = "back_to_sequences"
            tool[0] = cmd
            mean[0] = $2
            stddev[0] = $3
            min[0] = $7
            max[0] = $8
            min_mean = $2
            n = 1
        }
        NR > 2 {
            cmd = $1
            if (cmd ~ /merkurio/) cmd = "merkurio"
            else if (cmd ~ /st-/) cmd = "seqtool"
            else if (cmd ~ /fetch/) cmd = "fetch_reads"
            else if (cmd ~ /seqkit/) cmd = "seqkit"
            else if (cmd ~ /grep/) cmd = "grep"
            else if (cmd ~ /cookiecutter/) cmd = "cookiecutter"
            else if (cmd ~ /back_to_sequences/) cmd = "back_to_sequences"
            tool[n] = cmd
            mean[n] = $2
            stddev[n] = $3
            min[n] = $7
            max[n] = $8
            if ($2 < min_mean) min_mean = $2
            n++
        }
        END {
            for (i = 0; i < n; i++) {
                printf "| %s | %.3f | %.3f | %.3f | %.3f | %.2fx |\n", tool[i], mean[i], stddev[i], min[i], max[i], mean[i]/min_mean
            }
        }
        ' "$csv_file" >> ../results/summary.md
        
        echo "" >> ../results/summary.md
    fi
done

echo "Results processing complete! Check ../results/summary.md for the full report." 