"""
Generate master_s3_manifest.csv listing all staged output files and their S3 paths.
"""
import csv
def run_generate_s3_manifest(s3_base: str):
    log_section("Generate S3 Manifest")
    with open('file_names.txt') as f:
        files = sorted(set(line.strip() for line in f if line.strip()))
    output = 'master_s3_manifest.csv'
    with open(output, 'w', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=['filename', 's3_path'])
        writer.writeheader()
        for fname in files:
            writer.writerow({'filename': fname, 's3_path': f"{s3_base}/{fname}"})
    logger.info(f"Saved: {output} ({len(files)} entries)")

