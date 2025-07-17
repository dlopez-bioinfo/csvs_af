#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl = 2

process PREPARE_INTERSECTIONS {
    container "${params.container__pandas}"
    label 'high_mem'

    input:        
        tuple val(intersect_id), val(bed_list), val(region_list), path(sample_bed_file)
        path (vcf_list)
        path (vcf_index_list)

    output:
        tuple path("regions.bed"), path("vcfs")

    script:        
        """
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd
import sys
import os
import glob

id = "${intersect_id}"
bed_list = "${bed_list}".split(",")
region_list = ${region_list}

sample_bed_df = pd.read_csv(
    "${sample_bed_file}", 
    dtype=str, 
    header=None, 
    names=["sample_id", "vcf_path", "bed_hash", "bed_path", "gender"])

# Prepare regions file
with open("regions.bed", "w") as f:
    for region in region_list:
        print(region)
        line = "\\t".join(region.split("_"))
        f.write(f"{line}\\n")

# Create output directory
outdir = "vcfs"
os.makedirs(outdir, exist_ok=True)

# select VCF to be copied
filtered_ids = sample_bed_df.loc[sample_bed_df['bed_hash'].isin(bed_list), 'sample_id']

# Copy VCFs and indices to output directory
for sample_id in filtered_ids:
    pattern = f"{sample_id}*.vcf.gz*"
    for filepath in glob.glob(pattern):
        print(f"Copiando {filepath} → {outdir}")
        os.symlink(os.path.abspath(filepath), outdir + "/" + os.path.basename(filepath))
        """
}