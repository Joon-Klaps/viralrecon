#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import json
import argparse
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape


# ---------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------

def parser_args(args=None):
    Description = "Parse Sierra-local JSON reports and corresponding resistance and mutation tables to generate an HTML report."
    Epilog = """Example usage:
    python resistance_report.py --sierralocal_folder resistance_jsons --mutation_folder mutation_tables --resistance_folder resistance_tables --nextclade_folder nextclade_folder  --consensus_folder consensus --ivar_consensus_params "-t 0.8 -q 30 -m 50 -n N" --output_html resistance_report.html --template hiv_template_report.html --css hiv_template_report.css
    """
    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)

    parser.add_argument(
        "-s",
        "--sierralocal_folder",
        type=str,
        default="./sierralocal_json",
        help="Folder containing sierra-local JSON reports (default: ./sierralocal_json)",
    )
    parser.add_argument(
        "-m",
        "--mutation_folder",
        type=str,
        default="./mutation_tables",
        help="Folder containing mutation CSV files (default: ./mutation_tables)",
    )
    parser.add_argument(
        "-r",
        "--resistance_folder",
        type=str,
        default="./resistance_tables",
        help="Folder containing resistance CSV files (default: ./resistance_tables)",
    )
    parser.add_argument(
        "-n",
        "--nextclade_folder",
        type=str,
        default="./nextclade_reports",
        help="Folder containing nextclade CSV files (default: ./nextclade_reports)",
    )
    parser.add_argument(
        "-cn",
        "--consensus_folder",
        type=str,
        default="./consensus",
        help="Folder containing consensus files (default: ./consensus)",
    )
    parser.add_argument(
        "-i",
        "--ivar_consensus_params",
        type=str,
        help="Parameters used for ivar consensus calling",
    )
    parser.add_argument(
        "-o",
        "--output_html",
        type=str,
        help="Full path to output HTML report file.",
    )
    parser.add_argument(
        "-t",
        "--template",
        type=str,
        default="hiv_template_report.html",
        help="Full path to template HTML report file.",
    )
    parser.add_argument(
        "-cs",
        "--css",
        type=str,
        default="hiv_template_report.css",
        help="Full path to CSS file for the HTML report.",
    )
    return parser.parse_args(args)

def get_sample_number(sample_name, all_samples):
    """Return the sample index based on alphabetical order."""
    sorted_samples = sorted(all_samples)
    return sorted_samples.index(sample_name) + 1


def parse_sequence_summary(json_path, subtype_info=None, ivar_params=None):
    """
    Extract sequence summary information from Sierra-local JSON.
    - Lists each gene present (PR, RT, IN)
    - Detects missing nucleotide ranges from warnings
    - Detects missing codons at the start if firstAA != 1
    - Integrates subtype from Nextclade if available
    """

    with open(json_path, "r", encoding="utf-8") as f:
        data = json.load(f)[0]

    summary_lines = []
    warnings = {}

    # --- Parse validation warnings to find missing nucleotide ranges
    for warning in data.get("validationResults", []):
        msg = warning.get("message", "")
        match = re.search(r"\('([A-Z]{2})'.*?,\s*(\d+),\s*(\d+)\)", msg)
        if match:
            gene = match.group(1)
            start_nt = match.group(2)
            end_nt = match.group(3)
            warnings[gene] = (start_nt, end_nt)

    # --- Parse alignedGeneSequences
    for gene_entry in data.get("alignedGeneSequences", []):
        gene = gene_entry.get("gene", {}).get("name", "NA")
        first_aa = gene_entry.get("firstAA", "NA")
        last_aa = gene_entry.get("lastAA", "NA")

        # Build description for each gene sequence
        line = f"Sequence includes {gene}: codons {first_aa} - {last_aa}"

        missing_parts = []

        # If sequence does not start at codon 1, mention missing codons
        if isinstance(first_aa, int) or (isinstance(first_aa, str) and first_aa.isdigit()):
            first_aa = int(first_aa)
            if first_aa > 1:
                missing_parts.append(f"missing codons: 1-{first_aa - 1}")

        # If we have nucleotide warning info, include it too
        if gene in warnings:
            x, y = warnings[gene]
            missing_parts.append(f"missing nucleotides: {x} - {y}")

        # Append parentheses only if we have missing parts
        if missing_parts:
            line += " (" + ", ".join(missing_parts) + ")"

        summary_lines.append(line)

    # --- Add subtype info from Nextclade if provided
    summary_lines.append(f"Subtype: {subtype_info}")

    # --- Parse ivar consensus parameters if provided
    # Extract numeric values with regex
    match_t = re.search(r"-t\s*([\d.]+)", ivar_params)
    match_q = re.search(r"-q\s*(\d+)", ivar_params)
    match_m = re.search(r"-m\s*(\d+)", ivar_params)

    # Default values in case something is missing
    t_val = float(match_t.group(1))
    q_val = int(match_q.group(1))
    m_val = int(match_m.group(1))

    summary_lines.append(f"Minimum read depth: ≥{m_val}")
    summary_lines.append(f"Minimum allele frequency represented: ≥{t_val * 100:.0f}%")
    summary_lines.append(f"Minimum quality threshold: {q_val}")

    return summary_lines

def get_nextclade_subtype(nextclade_file, sample_name):
    """
    Extract the subtype (clade) assigned by Nextclade for the given sample.
    """
    if not nextclade_file or not os.path.exists(nextclade_file):
        return None
    try:
        df_next = pd.read_csv(nextclade_file, sep=';')
        # Attempt to find matching sample
        row = df_next[df_next["seqName"].str.contains(sample_name, case=False, na=False)]
        if not row.empty:
            subtype = row.iloc[0]["clade"]
            return subtype
    except Exception as e:
        print(f"⚠️ Could not parse Nextclade file for {sample_name}: {e}")
    return None


# ---------------------------------------------------------------------
# Batch processing
# ---------------------------------------------------------------------

def main():
    args = parser_args()

    # Detect all files
    mutation_files = sorted(glob.glob(os.path.join(args.mutation_folder, "*_mutation_table.csv")))
    resistance_files = sorted(glob.glob(os.path.join(args.resistance_folder, "*_resistance_table.csv")))
    json_files = sorted(glob.glob(os.path.join(args.sierralocal_folder, "*_resistance.json")))
    nextclade_files = sorted(glob.glob(os.path.join(args.nextclade_folder, "*.csv")))
    consensus_files = sorted(glob.glob(os.path.join(args.consensus_folder, "*.fa")))

    # --- Load CSS content
    with open(args.css, "r", encoding="utf-8") as css_file:
        css_content = css_file.read()

    all_samples_data = []

    for mut_file in mutation_files:
        df_mut = pd.read_csv(mut_file)
        if df_mut.empty:
            continue
        sample_name = df_mut["Sample_name"].iloc[0]

        # search corresponding resistance and JSON files
        res_file = next((r for r in resistance_files if sample_name in r), None)
        json_file = next((j for j in json_files if sample_name in j), None)
        nextclade_file = next((n for n in nextclade_files if sample_name in n), None)
        consensus_file = next((c for c in consensus_files if sample_name in c), None)
        consensus_seq = None
        if consensus_file and os.path.exists(consensus_file):
            with open(consensus_file, "r", encoding="utf-8") as f:
                consensus_seq = f.read().strip()

        if not res_file or not json_file:
            print(f"⚠️ Skipping {sample_name}: missing resistance or JSON file")
            continue

        # --- Extract subtype from Nextclade
        subtype = get_nextclade_subtype(nextclade_file, sample_name)

        # --- Parse sequence summary
        seq_summary = parse_sequence_summary(json_file, subtype_info=subtype, ivar_params=args.ivar_consensus_params)

        # Guardar toda la info en un dict
        all_samples_data.append({
            "sample_name": sample_name,
            "sequence_summary": seq_summary,
            "mutation_data": df_mut.to_dict(orient="records"),
            "resistance_data": pd.read_csv(res_file).to_dict(orient="records"),
            "consensus_genome": consensus_seq
        })

    # Ordenar alfabéticamente por nombre de muestra
    all_samples_data.sort(key=lambda x: x["sample_name"])

    # --- Load Jinja2 environment from the template's folder
    template_dir = os.path.dirname(os.path.abspath(args.template))
    template_name = os.path.basename(args.template)

    env = Environment(
        loader=FileSystemLoader(template_dir),
        autoescape=select_autoescape(['html', 'xml'])
    )
    template = env.get_template(template_name)

    # --- Render full HTML report
    html_content = template.render(
        all_samples=all_samples_data,
        css_content=css_content
    )

    output_html = args.output_html or "all_samples_report.html"
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_content)

    print(f"✅ Report generated: {output_html}")

if __name__ == "__main__":
    main()
