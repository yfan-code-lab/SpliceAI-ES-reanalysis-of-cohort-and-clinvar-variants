# SpliceAI re-analysis of exome sequencing (ES) data

Companion code for a "splicing-first" re-analysis framework: SpliceAI-based
re-prioritization of variants in exome sequencing (ES) data to (1) rescue
previously unsolved cases and (2) correct misinterpreted variant mechanisms
in already-diagnosed cases.


## Contents

| File | Description |
|---|---|
| `spliceAI_reA_ES.sh` | Cohort-level SpliceAI annotation and filtering pipeline (bash) |
| `plot_revel_ecdf.py` | REVEL score distribution comparison, Pathogenic vs VUS (Figure 4B) |

## Requirements

### `spliceAI_reA_ES.sh`
- Bash, `bcftools`, `bgzip`, `tabix`, `vcfanno`
- A conda environment named `vcfanno` containing the above tools
- `vcfanno_spliceai_conf.toml` (SpliceAI annotation configuration) and a
  `header` file (column header template) in the working directory

### `plot_revel_ecdf.py`
- Python 3 with `pandas`, `numpy`, `matplotlib`, `scipy`

## Usage

### 1. Cohort SpliceAI annotation — `spliceAI_reA_ES.sh`

The script expects a sample list file (passed as the first argument) with one
line per sample:

```
<SAMPLE_ID> <path/to/<SAMPLE_ID>.gst.vcf>
```

For each sample it:

1. Extracts OMIM-annotated variants from the gene-specific-target (gst) VCF;
2. Filters to `PASS`, `MSigDb != "."`, and rare alleles (`G.AF < 0.01` or missing);
3. Annotates with SpliceAI via `vcfanno`;
4. Extracts SpliceAI scores to `<SAMPLE_ID>_spliceai.tsv`;
5. Keeps variants with any SpliceAI delta score > 0.1 (`<SAMPLE_ID>_spliceai_sig.tsv`);
6. Merges all per-sample files into `combined.tsv`, then filters to any delta
   score > 0.5 (`combined_spliceai05.tsv`).

Before running, set the configuration variables at the top of the script:

| Variable | Meaning |
|---|---|
| `WORK_DIR` | Directory containing the input VCFs, vcfanno config, and header file |
| `OUTPUT_DIR` | Where per-sample result files are collected |
| `REF_SAMPLE_ID` | Any one of your sample IDs, used to build the header of `combined.tsv` |

```bash
bash spliceAI_reA_ES.sh sample_list.txt
```

### 2. REVEL distribution — `plot_revel_ecdf.py`

Expected repository layout:

```
.
├── input/
│   ├── clinvar_pathogenic_missense_spliceai_revel.tsv   # Pathogenic missense, SpliceAI + REVEL
│   └── clinvar_spiceai_filtered05_revel.tsv             # VUS, SpliceAI delta >= 0.5 + REVEL
└── output/
    ├── revel_density_pathogenic_vs_vus.png
    └── revel_stats_pathogenic_vs_vus.tsv
```

Run from the repository root:

```bash
python plot_revel_ecdf.py
```

Outputs a KDE density plot of REVEL scores with the BP4 recommended cutoff
(`REVEL = 0.644`) overlaid, plus summary statistics and a Mann-Whitney U test.


## Citation

If you use this code, please cite: *[citation — to be added]*

## License

*[License — to be added]* (e.g., MIT for the code; see the LICENSE file once added).
