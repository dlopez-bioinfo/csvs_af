#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process PROCESS_INTERVALS {
    container "${params.container__pandas}"
    label 'high_mem'

    input:        
        path intervals

        
    output:
        path("intervals.tsv")

    script:        
        """
python -c '
import pandas as pd

print(f"Loading input file: ${intervals}")
intervals = pd.read_csv("${intervals}", sep="\t", dtype=str)

required_columns = {"chrom", "start", "end", "list"}
if not required_columns.issubset(intervals.columns):
    missing = required_columns - set(intervals.columns)
    raise ValueError(f"Missing required columns {missing}")

print("Building interval IDs")
intervals["id"] = intervals["chrom"] + "_" + intervals["start"] + "_" + intervals["end"]

print("Exploding bed_list and grouping IDs")
intervals["list"] = intervals["list"].str.split(",")
intervals = intervals.explode("list")

# Agrupar sin separar por cromosoma
grouped = intervals.groupby("list", as_index=False).agg({"id": lambda x: ",".join(x)})

output_file = "intervals.tsv"
grouped.to_csv(output_file, sep="\t", index=False)
print(f"Output written to: {output_file}")
'
        """
}
