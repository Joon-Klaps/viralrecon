#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import glob
import json
import errno
import logging
import argparse
import pandas as pd

logger = logging.getLogger()
pd.set_option("display.max_columns", None)
pd.set_option("display.max_rows", None)


def parser_args(args=None):
    Description = "Parse Sierra-local JSON reports and generate a long table with resistance information."
    Epilog = """Example usage:
    python resistance_report.py --sierralocal_dir . --output_file sierralocal_resistance_table.csv
    """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument(
        "-sd",
        "--sierralocal_dir",
        type=str,
        default=".",
        help="Directory containing Sierra-local JSON files (default: '.').",
    )
    parser.add_argument(
        "-cd",
        "--codfreq_dir",
        type=str,
        default=".",
        help="Directory containing codfreq files (default: '.').",
    )
    parser.add_argument(
        "-fs",
        "--file_suffix",
        type=str,
        default="_resistance.json",
        help="Suffix to trim off JSON file name to obtain sample name (default: '_resistance.json').",
    )
    parser.add_argument(
        "-of",
        "--output_file",
        type=str,
        default="sierralocal_variants_table.csv",
        help="Full path to output CSV file (default: 'sierralocal_variants_table.csv').",
    )

    return parser.parse_args(args)


def make_dir(path):
    if not len(path) == 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


def get_file_dict(file_dir, file_suffix):
    """Return dictionary {sample_name: filepath} for all JSON reports in directory."""
    files = glob.glob(os.path.join(file_dir, f"*{file_suffix}"))
    samples = [os.path.basename(x).removesuffix(file_suffix) for x in files]
    return dict(zip(samples, files))

def parse_codfreq(codfreq_path):
    """Parse a .codfreq file into a pandas DataFrame."""
    try:
        df = pd.read_csv(codfreq_path)
        return df
    except Exception as e:
        logger.warning(f"Could not parse codfreq file {codfreq_path}: {e}")
        return pd.DataFrame()


def parse_sierra_json(sample_name, json_path):
    """Parse one Sierra-local JSON and return a pandas DataFrame with all mutations."""
    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)[0]

    rows = []

    if "alignedGeneSequences" not in data:
        logger.warning(f"No 'alignedGeneSequences' field found in {json_path}")
        return pd.DataFrame()

    for gene_entry in data["alignedGeneSequences"]:
        gene_name = gene_entry.get("gene", {}).get("name", "NA")
        aas_pattern = gene_entry.get("AAs", "NA")

        for mut in gene_entry.get("mutations", []):
            consensus = mut.get("consensus", "")
            aas = mut.get("AAs", "")
            pos = mut.get("position", "")
            mut_type = mut.get("primaryType", "NA")

            # Comentarios concatenados
            comments_list = mut.get("comments", [])
            comments_text = "; ".join([c.get("text", "") for c in comments_list]) if comments_list else ""

            # Si hay varios posibles aminoácidos, crear una fila por cada uno
            if len(aas) > 1:
                for aa in aas:
                    mut_text = f"{consensus}{pos}{aa}"
                    row = {
                        "Sample_name": sample_name,
                        "Gene_name": gene_name,
                        "Mutations": mut_text,
                        "Mutations_type": mut_type,
                        "Mutations_comments": comments_text,
                        "isInsertion": mut.get("isInsertion", False),
                        "isDeletion": mut.get("isDeletion", False),
                        "isApobecMutation": mut.get("isApobecMutation", False),
                        "isApobecDRM": mut.get("isApobecDRM", False),
                        "isUnusual": mut.get("isUnusual", False),
                        "isSDRM": mut.get("isSDRM", False),
                        "hasStop": mut.get("hasStop", False),
                        "Mutation_AF": "NA",
                        "Mutation_coverage": "NA",
                        "Sequenced": aas_pattern,
                    }
                    rows.append(row)
            else:
                # Si solo hay un aminoácido, dejar una sola fila
                mut_text = mut.get("text", "NA")
                row = {
                    "Sample_name": sample_name,
                    "Gene_name": gene_name,
                    "Mutations": mut_text,
                    "Mutations_type": mut_type,
                    "Mutations_comments": comments_text,
                    "isInsertion": mut.get("isInsertion", False),
                    "isDeletion": mut.get("isDeletion", False),
                    "isApobecMutation": mut.get("isApobecMutation", False),
                    "isApobecDRM": mut.get("isApobecDRM", False),
                    "isUnusual": mut.get("isUnusual", False),
                    "isSDRM": mut.get("isSDRM", False),
                    "hasStop": mut.get("hasStop", False),
                    "Mutation_AF": "NA",
                    "Mutation_coverage": "NA",
                    "Sequenced": aas_pattern,
                }
                rows.append(row)

    df = pd.DataFrame(rows)
    return df

def integrate_codfreq_info(df_json, codfreq_df):
    """Update df_json with Mutation_AF and Mutation_coverage from codfreq_df."""
    if codfreq_df.empty:
        return df_json

    updated_rows = []

    for _, row in df_json.iterrows():
        gene = row["Gene_name"]
        pos = row["Mutations"][1:-1]  # Extract position from format X123Y
        aa = row["Mutations"][-1]  # Extract mutated amino acid

        # Search in codfreq_df
        subset = codfreq_df[
            (codfreq_df["gene"] == gene) &
            (codfreq_df["position"] == int(pos)) &
            (codfreq_df["aa_codon"] == aa)
        ]

        if not subset.empty:
            total = subset["total"].iloc[0]
            coverage = subset["count"].sum()
            af = coverage / total if total > 0 else 0
            row["Mutation_AF"] = round(af, 5)
            row["Mutation_coverage"] = int(coverage)

        updated_rows.append(row)

    return pd.DataFrame(updated_rows)

def main(args=None):
    args = parser_args(args)

    # Create output directory if it doesn't exist
    out_dir = os.path.dirname(args.output_file)
    make_dir(out_dir)

    # Search JSON files
    json_files = get_file_dict(args.sierralocal_dir, args.file_suffix)

    if not json_files:
        logger.error(f"No JSON files found in directory: {args.sierralocal_dir}")
        sys.exit(1)

    # Process each sample
    all_samples = []
    for sample, json_path in sorted(json_files.items()):
        df_json = parse_sierra_json(sample, json_path)

        # Search for corresponding codfreq
        codfreq_path = os.path.join(args.codfreq_dir, f"{sample}.codfreq")
        if os.path.exists(codfreq_path):
            codfreq_df = parse_codfreq(codfreq_path)
            resistance_df = integrate_codfreq_info(df_json, codfreq_df)
        else:
            logger.warning(f"No codfreq file found for sample {sample}")

        if not resistance_df.empty:
            all_samples.append(resistance_df)
        else:
            logger.warning(f"No mutations found for sample {sample}")

    # Combinar todas las tablas
    if all_samples:
        merged_df = pd.concat(all_samples, ignore_index=True)
        merged_df.to_csv(args.output_file, index=False, encoding="utf-8-sig")
        print(f"✅ Sierra-local variants table with codfreq data saved to: {args.output_file}")
    else:
        logger.error("No valid data found in any JSON file.")
        sys.exit(1)


if __name__ == "__main__":
    sys.exit(main())
