# RISST

Rapid In Silico SeroTyping predicts capsular serotypes from assembled genomes
of *Streptococcus suis*, *Glaesserella parasuis*, and *Actinobacillus
pleuropneumoniae*. Choose one GenBank reference database for each run.

Serotyping is useful in diagnostics, epidemiological surveillance, and
microbiology research. For these three swine pathogens, the loci used for
typing encode capsular polysaccharide (CPS). RISST uses curated CPS references
to screen genome assemblies in batches and reports the evidence behind each
candidate call, including coverage, identity, and the genes found in the locus.

## Why I wrote RISST

I started by putting together three reference databases for use with
[Kaptive](https://github.com/katholt/Kaptive). As I learned more Python, I wanted
to see whether I could turn those databases into a standalone tool. I also
liked the name RISST, and it was satisfying to see the idea become a program
I could use.

## Dependencies

- Python 3.10 or later
- minimap2 available on `PATH`
- Biopython 1.85 or later and NumPy
- Matplotlib 3.6 or later for `-f` gene maps

```sh
python -m pip install "biopython>=1.85" numpy
# For gene maps
python -m pip install "matplotlib>=3.6"
```

RISST 3 replaces BLAST+ and Prodigal with minimap2. The `risst.tr` training
file is obsolete.

## Usage

```text
python RISST.py -i INPUT [INPUT ...] -r REFERENCE [-o OUTPUT] [-t THREADS]
                [--min_locus_cov PCT] [--min_locus_id PCT] [--min-cov PCT]
                [--n-best N] [--percent-expected PCT] [--no_cps_sequence]
                [-f] [-v] [-h]
```

| Option | Description |
| --- | --- |
| `-i`, `--input` | One or more FASTA or FASTA.gz assemblies. Quoted glob patterns also work. |
| `-r`, `--reference` | GenBank file with CDS features and a unique `K locus:`, `O locus:`, or `Serovar:` source note per record. |
| `-o`, `--output` | Output TSV. Defaults to `RISST_output.txt`. |
| `-t`, `--threads` | Threads for minimap2. Defaults to the available CPU count, capped at 4. |
| `--min_locus_cov`, `--min_locus_id` | Minimum expected-gene nucleotide coverage and mean protein identity, respectively. Both default to 95%. |
| `--min-cov` | Minimum CDS alignment coverage for scoring candidate loci. Defaults to 50%. |
| `--n-best` | Maximum number of candidates to align as full loci. Defaults to 2. |
| `--percent-expected` | A typeable call requires more than this percentage of distinct expected CDS inside the locus. Defaults to 50%. |
| `--no_cps_sequence` | Skip CPS FASTA export. |
| `-f`, `--figure` | Export a 300 dpi PNG and an editable SVG for each sample with reconstructed fragments. |
| `-v`, `--version` | Print the version and exit. |
| `-h`, `--help` | Show command-line help. |

The old `--minimum_piece`, `--min_gene_cov`, and `--min_gene_id` flags have
been retired. Candidate screening uses `--min-cov`. The final call uses the
locus and expected-gene thresholds above.

### Reference databases

| Species | Reference |
| --- | --- |
| *Streptococcus suis* | [CPS reference](Streptococcus_suis_cps_locus_reference.gbk) |
| *Glaesserella parasuis* | [Primary CPS reference](Glaesserella_parasuis/Glaesserella_parasuis_cps_locus_primary_reference.gbk) |
| *Actinobacillus pleuropneumoniae* | [CPS reference](A_pleuropneumoniae_cps_locus_reference.gbk) |

### Quick start

From the RISST directory, replace `path/to/assemblies/` with the directory
containing your *G. parasuis* FASTA files. For a different species, select its
reference from the table above.

```sh
python RISST.py -i "path/to/assemblies/*.fasta" -r Glaesserella_parasuis/Glaesserella_parasuis_cps_locus_primary_reference.gbk -o results/serotypes.tsv -f
```

## Serotype prediction

Minimap2 first aligns the reference CDS sequences to each assembly. RISST ranks
candidate loci by the sum of their CDS alignment scores weighted by the fraction
of expected genes found. It then aligns the top candidates as full loci,
reconstructs fragments from the best match, and evaluates CDS hits by protein
alignment.

## Output

RISST prints one call per sample to the terminal. The examples below come
from a run with the *G. parasuis* primary reference and default typing
thresholds. The repository includes the example report excerpts and gene map.

```text
GCA_000021885.1 : serotype 5
GCA_009796465.1 : serotype 5
GCA_000478405.1 : serotype 4?
```

### TSV report

RISST rewrites the TSV on each run and writes one row per input assembly. The
first seven columns of the example output are shown below.

| Sample | Best serotype | Confidence | Expected gene nucleotide coverage (%) | Expected gene mean protein identity (%) | Locus pieces | Problems |
| --- | --- | --- | ---: | ---: | ---: | --- |
| GCA_000021885.1 | 5 | Typeable | 99.79 | 98.96 | 1 | |
| GCA_009796465.1 | 5 | Typeable | 99.80 | 98.72 | 2 | `?2!` |
| GCA_000478405.1 | 4? | Untypeable | 99.90 | 93.11 | 1 | |

A question mark on the serotype marks a candidate that fell below a typing
threshold. In the last row, 93.11% mean protein identity is below the default
95% threshold. The second row shows that a locus can be typeable even when
it spans two fragments. RISST reports `NA` and `Untypeable` when it finds no
usable match, then continues with the next assembly.

Expected-gene nucleotide coverage divides the summed lengths of expected CDS
hits inside the reconstructed fragments by the total expected reference CDS
length. Multiple hits can raise this value above 100%. Mean protein identity
averages the protein alignment identities of those hits. RISST 2 reported
whole-locus BLAST percentages.

The remaining columns list expected and other-locus gene hits inside and
outside the reconstructed fragments, plus missing expected genes. Hit entries
include protein identity, protein coverage, assembly coordinates, and strand.

| Problem | Meaning |
| --- | --- |
| `?N` | The locus has `N` reconstructed fragments, where `N` differs from 1. |
| `-` | At least one expected gene is missing. |
| `+` | At least one hit to another reference locus lies inside the reconstructed fragments. |
| `*` | At least one expected hit inside the locus has protein identity below 80% and protein coverage of at least 90%. |
| `!` | At least one gene hit is flagged as truncated. |

### CPS sequences and gene maps

RISST writes a FASTA per sample with reconstructed fragments in
`result_cps_dict/` alongside the TSV. Each fragment has its own FASTA record
in inferred locus orientation. FASTA headers and TSV gene hits use 1-based
inclusive assembly coordinates. For `GCA_009796465.1`, the two FASTA headers are:

```text
>WCJJ01000005.1:440-1158(+)
>WCJJ01000006.1:1-14755(+)
```

With `-f`, RISST also writes PNG and SVG maps in `gene_structure_map/`
alongside the TSV. The map arranges fragments on one track in reference-locus
order. A `//` marks a break of unknown spacing. Blue arrows show expected
genes and orange arrows show hits to genes from other reference loci within
the plotted fragments.

## Further

You can build a GenBank reference for another species or a different gene
cluster with known types. Use one record per reference locus, following the
format of the supplied databases. Each record needs the full locus nucleotide
sequence, annotated CDS features, and a unique type label in the `source`
feature, such as `/note="K locus: 1"`. RISST derives the protein sequences
from the annotated CDS. If you download annotated records from NCBI, choose
the full GenBank format and ensure each record covers the intended locus.
Pass the resulting file to RISST with `-r`.
