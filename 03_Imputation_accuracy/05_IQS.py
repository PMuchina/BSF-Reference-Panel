import cyvcf2
import pandas as pd
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import os
import json

# 1. Function to load sample names from file
def load_sample_names(sample_file):
    with open(sample_file, 'r') as f:
        samples = [line.strip() for line in f.readlines()]
    return samples

# 2. Function to extract true genotypes (GT) and imputed genotype probabilities (GP) from VCF chunk pair
def process_vcf_chunk_pair(true_vcf_file, imputed_vcf_file, samples, start, end, output_dir, chunk_id):
    true_vcf = cyvcf2.VCF(true_vcf_file)
    imputed_vcf = cyvcf2.VCF(imputed_vcf_file)
    chunk_results = []

    for i, (true_record, imputed_record) in enumerate(zip(true_vcf, imputed_vcf)):
        if i < start:
            continue
        if i >= end:
            break

        if true_record.CHROM != imputed_record.CHROM or true_record.POS != imputed_record.POS:
            raise ValueError(f"Records do not match: {true_record.CHROM}:{true_record.POS} != {imputed_record.CHROM}:{imputed_record.POS}")

        chrom = true_record.CHROM
        pos = true_record.POS

        # Extract True Genotypes (GT)
        true_genotypes = {}
        for sample, gt in zip(true_vcf.samples, true_record.genotypes):
            if sample in samples:
                gt_str = '/'.join(map(str, gt[:2]))
                dosage = (
                    2 if gt_str == "0/0" else
                    1 if gt_str in ["0/1", "1/0"] else
                    0 if gt_str == "1/1" else -1
                )
                true_genotypes[sample] = {'GT': gt_str, 'Dosage': dosage}

        # Extract Imputed Genotype Probabilities (GP)
        if "GP" in imputed_record.FORMAT:
            imputed_genotypes = {}
            for sample, gp in zip(imputed_vcf.samples, imputed_record.format("GP")):
                if sample in samples:
                    imputed_genotypes[sample] = gp.tolist()  # Ensure it's JSON serializable
        else:
            continue

        # Extract MAF
        maf = true_record.INFO.get('MAF', None)

        # Combine data for this SNP
        if true_genotypes and imputed_genotypes:
            chunk_results.append({
                "SNP": (chrom, pos),
                "True": json.dumps(true_genotypes),  # Save as JSON string
                "Imputed": json.dumps(imputed_genotypes),  # Save as JSON string
                "MAF": maf
            })

    # Save chunk results to a file
    output_file = os.path.join(output_dir, f"chunk_{chunk_id}.csv")
    pd.DataFrame(chunk_results).to_csv(output_file, index=False)

# 3. Function to calculate IQS for a single SNP
def calculate_iqs_for_snp(snp_data, samples):
    true_info = snp_data["True"]
    imputed_info = snp_data["Imputed"]
    maf = snp_data["MAF"]

    P0 = Pc = 0  # Observed and expected agreement
    p11 = p22 = p33 = 0
    N1_i = N2_i = N3_i = 0
    N1_j = N2_j = N3_j = 0
    N = len(samples)

    for sample in samples:
        if sample not in true_info or sample not in imputed_info:
            continue

        true_gt = true_info[sample]['Dosage']
        imputed_probs = imputed_info[sample]

        true_vec = [0, 0, 0]  # [AA, AB, BB]
        if true_gt == 2:
            true_vec[0] = 1
        elif true_gt == 1:
            true_vec[1] = 1
        elif true_gt == 0:
            true_vec[2] = 1

        p11 += true_vec[0] * imputed_probs[0]
        p22 += true_vec[1] * imputed_probs[1]
        p33 += true_vec[2] * imputed_probs[2]

        N1_i += true_vec[0]
        N2_i += true_vec[1]
        N3_i += true_vec[2]

        N1_j += imputed_probs[0]
        N2_j += imputed_probs[1]
        N3_j += imputed_probs[2]

    P0 = (p11 + p22 + p33) / N
    Pc = ((N1_i * N1_j) + (N2_i * N2_j) + (N3_i * N3_j)) / (N ** 2)

    IQS = (P0 - Pc) / (1 - Pc) if Pc != 1 else 0

    return {"SNP": snp_data["SNP"], "MAF": maf, "IQS": max(0, min(IQS, 1))}

# 4. Main function to process the VCF in chunks and calculate IQS
def process_vcf_in_chunks(true_vcf_file, imputed_vcf_file, samples, output_dir, chunk_size=1000):
    true_vcf = cyvcf2.VCF(true_vcf_file)
    total_records = sum(1 for _ in true_vcf)
    true_vcf.close()  # Close and reopen to reset iterator
    true_vcf = cyvcf2.VCF(true_vcf_file)

    chunks = [(i, min(i + chunk_size, total_records)) for i in range(0, total_records, chunk_size)]

    with ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(process_vcf_chunk_pair, true_vcf_file, imputed_vcf_file, samples, start, end, output_dir, chunk_id)
            for chunk_id, (start, end) in enumerate(chunks)
        ]

        for future in futures:
            future.result()

# 5. Function to combine chunk results and calculate IQS
def combine_and_calculate_iqs(output_dir, samples):
    all_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.startswith("chunk_")]
    results = []

    for file in all_files:
        chunk_data = pd.read_csv(file).to_dict('records')
        for snp_data in chunk_data:
            # Parse JSON strings back into dictionaries
            snp_data["True"] = json.loads(snp_data["True"])
            snp_data["Imputed"] = json.loads(snp_data["Imputed"])
            results.append(calculate_iqs_for_snp(snp_data, samples))

    return results

# 6. Function to bin MAF values and calculate average IQS
def bin_maf_and_calculate_iqs(results):
    maf_bins = [(0, 0.05), (0.05, 0.1), (0.1, 0.15), (0.15, 0.2), (0.2, 0.25),
                (0.25, 0.3), (0.3, 0.35), (0.35, 0.4), (0.4, 0.45), (0.45, 0.5)]

    bin_results = {f"{bin[0]} - {bin[1]}": [] for bin in maf_bins}

    for result in results:
        maf = result['MAF']
        iq = result['IQS']

        if maf is not None:
            for bin_range in maf_bins:
                if bin_range[0] <= maf < bin_range[1]:
                    bin_key = f"{bin_range[0]} - {bin_range[1]}"
                    bin_results[bin_key].append(iq)
                    break

    bin_avg_iqs = {}
    for bin_key, values in bin_results.items():
        if values:
            bin_avg_iqs[bin_key] = np.mean(values)

    return bin_avg_iqs

# Main function
# Adjust/change files
def main():
    sample_file = "sample_names.txt"
    true_vcf_file = "True_data.bcf"
    imputed_vcf_file = "0.5x_imputed.Glimpse2.bcf"
    output_dir = "output_chunks"
    chunk_size = 1000000

    # Create output directory if not exists
    os.makedirs(output_dir, exist_ok=True)

    # Load sample names
    samples = load_sample_names(sample_file)

    # Step 1: Process VCF in chunks and save intermediate results
    print("Processing VCFs in chunks...")
    process_vcf_in_chunks(true_vcf_file, imputed_vcf_file, samples, output_dir, chunk_size)

    # Step 2: Combine chunk results and calculate IQS
    print("Combining results and calculating IQS...")
    iqs_results = combine_and_calculate_iqs(output_dir, samples)

    # Step 3: Bin MAF values and calculate average IQS
    print("Binning MAF values and calculating average IQS...")
    maf_bins_avg_iqs = bin_maf_and_calculate_iqs(iqs_results)

    # Print the results
    print("MAF bins and average IQS values:")
    for bin_range, avg_iqs in maf_bins_avg_iqs.items():
        print(f"MAF range {bin_range}: Average IQS = {avg_iqs}")

    # Optionally save the IQS results and MAF bin averages
    results_file = os.path.join(output_dir, "iqs_results.csv")
    pd.DataFrame(iqs_results).to_csv(results_file, index=False)
    print(f"Saved IQS results to {results_file}")

    maf_bins_file = os.path.join(output_dir, "maf_bins_avg_iqs.csv")
    pd.DataFrame(list(maf_bins_avg_iqs.items()), columns=["MAF Range", "Average IQS"]).to_csv(maf_bins_file, index=False)
    print(f"Saved MAF bins average IQS to {maf_bins_file}")

if __name__ == "__main__":
    main()
