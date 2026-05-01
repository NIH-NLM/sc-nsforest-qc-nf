"""
Generate master_s3_manifest.csv listing all staged output files and their S3 paths.
"""
import csv
import os
from .common_utils import log_section, logger

def run_generate_s3_manifest(s3_base: str):
    log_section("Generate S3 Manifest")
    files = sorted(
        f for f in os.listdir('.')
        if not f.startswith('.') and f != 'master_s3_manifest.csv'
    )
    output = 'master_s3_manifest.csv'
    with open(output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['filename', 's3_path'])
        writer.writeheader()
        for fname in files:
            writer.writerow({'filename': fname, 's3_path': f"{s3_base}/{fname}"})
    logger.info(f"Saved: {output} ({len(files)} entries)")
    
