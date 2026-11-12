"""Two-stage locus typing adapted from leptyper by Genglin Guo (GPL-3.0).

The original workflow draws on Kaptive.
"""

from collections import defaultdict
from dataclasses import dataclass, field
from functools import cached_property
import re
import warnings

import numpy as np
from Bio import BiopythonWarning, SeqIO
from Bio.Align import PairwiseAligner
from Bio.Seq import Seq

from risst_alignment import Alignment, Assembly, RISSTError


@dataclass
class Gene:
    label: str
    start: int
    strand: str
    sequence: Seq
    symbol: str = ''

    @cached_property
    def protein(self):
        return translate(self.sequence)


@dataclass
class Locus:
    label: str
    sequence: Seq
    genes: dict[str, Gene]


@dataclass
class Database:
    loci: dict[str, Locus]
    gene_loci: dict[str, str]

    @cached_property
    def gene_fasta(self):
        return ''.join(f'>{key}\n{gene.sequence}\n'
                       for locus in self.loci.values() for key, gene in locus.genes.items())


def load_database(path):
    loci, gene_loci, labels = {}, {}, set()
    for index, record in enumerate(SeqIO.parse(path, 'genbank'), 1):
        location = f'{path}, record {record.id}'
        notes = (note for feature in record.features if feature.type == 'source'
                 for note in feature.qualifiers.get('note', []))
        names = {match[1].strip() for note in notes
                 if (match := re.fullmatch(r'(?:[KO] locus|Serovar):\s*(.+)', note.strip(), re.I))}
        if len(names) != 1 or not next(iter(names)):
            raise RISSTError(f'{location}: expected one K locus, O locus or Serovar label')
        label = names.pop()
        if label in labels:
            raise RISSTError(f'{location}: duplicate locus label {label!r}')
        labels.add(label)
        key = f'L{index}'
        genes = {}
        for number, feature in enumerate((f for f in record.features if f.type == 'CDS'), 1):
            if feature.location is None:
                raise RISSTError(f'{location}: CDS {number} has no valid location')
            sequence = feature.extract(record.seq)
            if not sequence:
                raise RISSTError(f'{location}: CDS {number} has no sequence')
            gene_name = feature.qualifiers.get('gene', [''])[0]
            display = f'{label}_{number:04d}' + (f'_{gene_name}' if gene_name else '')
            gene_key = f'{key}_G{number}'
            genes[gene_key] = Gene(display, int(feature.location.start),
                                   '+' if feature.location.strand == 1 else '-', sequence, gene_name)
            gene_loci[gene_key] = key
        if not genes:
            raise RISSTError(f'{location}: no CDS features')
        loci[key] = Locus(label, record.seq, genes)
    if not loci:
        raise RISSTError(f'{path}: no GenBank records')
    return Database(loci, gene_loci)


@dataclass
class GeneHit:
    key: str
    gene: Gene
    alignment: Alignment
    start: int
    end: int
    identity: float
    coverage: float
    truncated: bool


@dataclass
class Piece:
    contig: str
    start: int
    end: int
    strand: str = '+'
    expected: list[GeneHit] = field(default_factory=list)
    unexpected: list[GeneHit] = field(default_factory=list)


@dataclass
class TypingResult:
    locus: Locus
    pieces: list[Piece]
    expected_outside: list[GeneHit]
    unexpected_outside: list[GeneHit]
    missing: list[str]

    @cached_property
    def expected_inside(self):
        return [hit for piece in self.pieces for hit in piece.expected]

    @cached_property
    def unexpected_inside(self):
        return [hit for piece in self.pieces for hit in piece.unexpected]

    @property
    def identity(self):
        return float(np.mean([hit.identity for hit in self.expected_inside])) if self.expected_inside else 0.0

    @property
    def coverage(self):
        return 100 * sum(hit.end - hit.start for hit in self.expected_inside) / sum(
            len(gene.sequence) for gene in self.locus.genes.values())

    def typeable(self, min_coverage, min_identity, percent_expected):
        found = 100 * len({hit.key for hit in self.expected_inside}) / len(self.locus.genes)
        return self.coverage >= min_coverage and self.identity >= min_identity and found > percent_expected

    @property
    def problems(self):
        return ''.join((
            f'?{len(self.pieces)}' if len(self.pieces) != 1 else '',
            '-' if self.missing else '',
            '+' if self.unexpected_inside else '',
            '*' if any(hit.coverage >= 90 and hit.identity < 80 for hit in self.expected_inside) else '',
            '!' if any(hit.truncated for hits in (self.expected_inside, self.unexpected_inside,
                                                 self.expected_outside, self.unexpected_outside)
                       for hit in hits) else '',
        ))


def translate(sequence):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', BiopythonWarning)
        return sequence.translate(table=11, to_stop=True)


def compare_gene(assembly, alignment, gene, expected, inside, aligner):
    sequence = assembly.sequence(alignment.contig, alignment.start, alignment.end, alignment.strand)
    frame, protein = max(enumerate(translate(sequence[i:]) for i in range(3)),
                         key=lambda item: len(item[1]))
    start = alignment.start + frame if alignment.strand == '+' else alignment.start
    end = alignment.end - frame if alignment.strand == '-' else alignment.end
    identity = coverage = 0.0
    truncated = False
    if len(protein) > 1 and len(gene.protein) > 1:
        match = next(iter(aligner.align(gene.protein, protein)), None)
        if match is not None:
            identity = 100 * match.counts().identities / match.shape[1]
            coverage = 100 * len(protein) / len(gene.protein)
            truncated = coverage < 95 and not alignment.partial and (expected or inside)
    return GeneHit(alignment.query, gene, alignment, start, end, identity, coverage, truncated)


def overlap(first_start, first_end, second_start, second_end):
    return max(0, min(first_end, second_end) - max(first_start, second_start))


def conflicts(hit, kept):
    return (hit.contig == kept.contig
            and overlap(hit.start, hit.end, kept.start, kept.end) >= 0.1 * hit.block_length)


def selected_gene_hits(alignments, locus):
    expected = [hit for hit in alignments if hit.query in locus.genes]
    other = sorted((hit for hit in alignments if hit.query not in locus.genes),
                   key=lambda hit: hit.matches, reverse=True)
    retained = []
    for hit in other:
        if not any(conflicts(hit, previous) for previous in retained):
            retained.append(hit)
    yield from expected
    yield from (hit for hit in retained if not any(conflicts(hit, previous) for previous in expected))


def reconstruct_pieces(alignments, locus):
    by_contig = defaultdict(list)
    for hit in alignments:
        by_contig[hit.contig].append((hit.start, hit.end))
    pieces = defaultdict(list)
    for contig, ranges in by_contig.items():
        for start, end in sorted(ranges):
            if pieces[contig] and start <= pieces[contig][-1].end + len(locus.sequence):
                pieces[contig][-1].end = max(pieces[contig][-1].end, end)
            else:
                pieces[contig].append(Piece(contig, start, end))
    return pieces


def type_assembly(assembly: Assembly, database: Database, min_cov: float, n_best: int):
    by_gene = defaultdict(list)
    for hit in assembly.map(database.gene_fasta):
        by_gene[hit.query].append(hit)

    scores = defaultdict(int)
    found = defaultdict(int)
    for query, hits in by_gene.items():
        hit = max(hits, key=lambda item: item.matches)
        if 100 * hit.block_length / hit.query_length >= min_cov:
            key = database.gene_loci[query]
            scores[key] += hit.score
            found[key] += 1
    keys = list(database.loci)
    rankings = np.array([scores[key] * found[key] / len(database.loci[key].genes) for key in keys])
    candidates = [keys[index] for index in np.argsort(rankings)[::-1][:n_best] if rankings[index] > 0]
    if not candidates:
        return None

    full_hits = defaultdict(list)
    query = ''.join(f'>{key}\n{database.loci[key].sequence}\n' for key in candidates)
    for hit in assembly.map(query):
        full_hits[hit.query].append(hit)
    if not full_hits:
        return None
    best = max((key for key in candidates if full_hits[key]), key=lambda key: sum(
        hit.score * (hit.matches / hit.block_length) * ((hit.query_end - hit.query_start) / hit.query_length)
        for hit in full_hits[key]))
    locus = database.loci[best]
    pieces_by_contig = reconstruct_pieces(full_hits[best], locus)

    expected_outside, unexpected_outside = [], []
    aligner = PairwiseAligner(scoring='blastp', mode='local')
    for alignment in selected_gene_hits((hit for hits in by_gene.values() for hit in hits), locus):
        key = database.gene_loci[alignment.query]
        expected = key == best
        piece = next((piece for piece in pieces_by_contig.get(alignment.contig, [])
                      if overlap(piece.start, piece.end, alignment.start, alignment.end)), None)
        hit = compare_gene(assembly, alignment, database.loci[key].genes[alignment.query],
                           expected, piece is not None, aligner)
        if piece is not None:
            piece.start = min(piece.start, hit.start)
            piece.end = max(piece.end, hit.end)
            (piece.expected if expected else piece.unexpected).append(hit)
        elif hit.identity >= 80:
            (expected_outside if expected else unexpected_outside).append(hit)

    pieces = [piece for contig_pieces in pieces_by_contig.values()
              for piece in contig_pieces if piece.expected]
    for piece in pieces:
        piece.strand = '+' if any(hit.alignment.strand == hit.gene.strand
                                  for hit in piece.expected) else '-'
    pieces.sort(key=lambda piece: min(hit.gene.start for hit in piece.expected))
    detected = {hit.key for piece in pieces for hit in piece.expected} | {hit.key for hit in expected_outside}
    missing = [gene.label for key, gene in locus.genes.items() if key not in detected]
    return TypingResult(locus, pieces, expected_outside, unexpected_outside, missing)
