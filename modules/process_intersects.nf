#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process PROCESS_INTERSECTS {
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

intervals = pd.read_csv("${intervals}", sep="\t", dtype=str, usecols=range(5))
intervals["regions"] = intervals["chrom"] + "_" + intervals["start"] + "_" + intervals["end"]
intervals = intervals.groupby("list").agg({"regions": list}).reset_index(drop=False)
intervals.to_csv("intervals.tsv", sep="\t", index=True, header=False)
'
        """
}