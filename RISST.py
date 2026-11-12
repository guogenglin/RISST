"""Rapid In Silico SeroTyping."""

import argparse
import csv
import glob
import gzip
import os
from pathlib import Path
import shutil
import sys
from tempfile import TemporaryDirectory

from Bio import SeqIO

from risst_alignment import Assembly, RISSTError
from risst_typing import load_database, type_assembly

VERSION = '3.0.0'
HEADER = (
    'Sample', 'Best serotype', 'Confidence', 'Expected gene nucleotide coverage (%)',
    'Expected gene mean protein identity (%)', 'Locus pieces', 'Problems',
    'Expected genes inside locus', 'Expected genes outside locus', 'Missing genes',
    'Other genes inside locus', 'Other genes outside locus',
)


def percentage(text):
    value = float(text)
    if not 0 <= value <= 100:
        raise argparse.ArgumentTypeError('expected a percentage from 0 to 100')
    return value


def positive_int(text):
    value = int(text)
    if value < 1:
        raise argparse.ArgumentTypeError('expected an integer greater than zero')
    return value


def parse_args(argv):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', '--input', nargs='+', required=True, help='FASTA assemblies (plain or gzip)')
    parser.add_argument('-r', '--reference', required=True, type=Path, help='GenBank locus reference')
    parser.add_argument('-o', '--output', type=Path, default=Path('RISST_output.txt'), help='TSV output')
    parser.add_argument('-t', '--threads', type=positive_int, default=min(os.cpu_count() or 1, 4))
    parser.add_argument('--min_locus_cov', '--min-locus-cov', type=percentage, default=95.0)
    parser.add_argument('--min_locus_id', '--min-locus-id', type=percentage, default=95.0)
    parser.add_argument('--min-cov', type=percentage, default=50.0, help='CDS coverage for candidate scoring')
    parser.add_argument('--n-best', type=positive_int, default=2, help='number of loci to align')
    parser.add_argument('--percent-expected', type=percentage, default=50.0,
                        help='require more than this percentage of expected genes inside the locus')
    parser.add_argument('--no_cps_sequence', action='store_true', help='suppress CPS FASTA export')
    parser.add_argument('-f', '--figure', action='store_true', help='export a CPS gene map')
    parser.add_argument('-v', '--version', action='version', version=f'RISST v{VERSION}')
    return parser.parse_args(argv)


def input_samples(patterns):
    samples = {}
    for pattern in patterns:
        paths = [pattern] if Path(pattern).is_file() else sorted(glob.glob(pattern))
        if not paths:
            raise RISSTError(f'No FASTA input matched {pattern!r}')
        for entry in paths:
            path = Path(entry)
            if not path.is_file():
                raise RISSTError(f'Not a FASTA file: {path}')
            basename = path.stem if path.suffix.lower() == '.gz' else path.name
            name = Path(basename).stem
            key = name.casefold()
            if key in samples:
                raise RISSTError(f'Duplicate sample name {name!r}: {samples[key][1]} and {path}')
            samples[key] = (name, path)
    return samples.values()


def prepare_assembly(name, path, temporary, threads):
    opener = gzip.open if path.suffix.lower() == '.gz' else open
    contigs = {}
    with opener(path, 'rt', encoding='utf-8') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            if not record.id or not record.seq:
                raise RISSTError(f'{path}: missing contig ID or empty sequence')
            if record.id in contigs:
                raise RISSTError(f'{path}: duplicate contig ID {record.id!r}')
            contigs[record.id] = record.seq
    if not contigs:
        raise RISSTError(f'{path}: no FASTA sequences')
    # Normalize contig headers and avoid passing original paths to the aligner.
    fasta = temporary / 'assembly.fasta'
    with fasta.open('w', encoding='utf-8') as handle:
        for contig, sequence in contigs.items():
            handle.write(f'>{contig}\n{sequence}\n')
    assembly = Assembly(name, temporary / 'assembly.mmi', contigs, threads)
    assembly.build_index(fasta)
    return assembly


def gene_details(hits):
    details = []
    for hit in sorted(hits, key=lambda hit: (hit.gene.start, hit.alignment.contig, hit.start)):
        alignment = hit.alignment
        detail = (f'{hit.gene.label},{hit.identity:.2f}%,{hit.coverage:.2f}%,'
                  f'{alignment.contig}:{hit.start + 1}-{hit.end}({alignment.strand})')
        if alignment.partial:
            detail += ',partial'
        if hit.truncated:
            detail += ',truncated'
        if hit.identity < 80:
            detail += ',below_id_threshold'
        details.append(detail)
    return ';'.join(details)


def result_row(name, result, args):
    if result is None:
        return (name, 'NA', 'Untypeable', '0.00', '0.00', 0, '', '', '', '', '', '')
    typeable = result.typeable(args.min_locus_cov, args.min_locus_id, args.percent_expected)
    return (
        name, result.locus.label if typeable else result.locus.label + '?',
        'Typeable' if typeable else 'Untypeable', f'{result.coverage:.2f}', f'{result.identity:.2f}',
        len(result.pieces), result.problems, gene_details(result.expected_inside),
        gene_details(result.expected_outside), ';'.join(result.missing),
        gene_details(result.unexpected_inside), gene_details(result.unexpected_outside),
    )


def write_cps(path, assembly, pieces):
    path.parent.mkdir(exist_ok=True)
    with path.open('w', encoding='utf-8') as handle:
        for piece in pieces:
            sequence = assembly.sequence(piece.contig, piece.start, piece.end, piece.strand)
            handle.write(f'>{piece.contig}:{piece.start + 1}-{piece.end}({piece.strand})\n')
            for start in range(0, len(sequence), 80):
                handle.write(f'{sequence[start:start + 80]}\n')


def main(argv=None):
    args = parse_args(argv)
    try:
        samples = input_samples(args.input)
        inputs = {args.reference.resolve(), *(path.resolve() for _, path in samples)}
        if args.output.resolve() in inputs:
            raise RISSTError(f'Output path would overwrite an input: {args.output}')
        for name, _ in samples:
            exports = []
            if not args.no_cps_sequence:
                exports.append(args.output.parent / 'result_cps_dict' / f'{name}.fasta')
            if args.figure:
                figure_path = args.output.parent / 'gene_structure_map' / f'{name}.png'
                exports.extend((figure_path, figure_path.with_suffix('.svg')))
            for path in exports:
                if path.resolve() in inputs or path.resolve() == args.output.resolve():
                    raise RISSTError(f'Export path conflicts with an input or the result table: {path}')
        database = load_database(args.reference)
        if shutil.which('minimap2') is None:
            raise RISSTError('minimap2 is required on PATH')
        if args.figure:
            try:
                import matplotlib
                matplotlib.use('Agg')
                from risst_figure import write_figure
            except ImportError as error:
                raise RISSTError('Matplotlib is required for --figure') from error
        args.output.parent.mkdir(parents=True, exist_ok=True)
        with args.output.open('w', encoding='utf-8', newline='') as handle:
            writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
            writer.writerow(HEADER)
            for name, path in samples:
                with TemporaryDirectory(prefix='risst-') as temporary:
                    assembly = prepare_assembly(name, path, Path(temporary), args.threads)
                    result = type_assembly(assembly, database, args.min_cov, args.n_best)
                    row = result_row(name, result, args)
                    if result is not None and result.pieces:
                        if not args.no_cps_sequence:
                            write_cps(args.output.parent / 'result_cps_dict' / f'{name}.fasta',
                                      assembly, result.pieces)
                        if args.figure:
                            write_figure(args.output.parent / 'gene_structure_map' / f'{name}.png',
                                         result, name, row[1])
                    writer.writerow(row)
                    print(f'{name} : serotype {row[1]}')
    except (RISSTError, OSError, ValueError) as error:
        print(f'RISST: {error}', file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
