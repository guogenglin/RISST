"""Minimap2 adapter based on leptyper by Genglin Guo (GPL-3.0)."""

from dataclasses import dataclass
from pathlib import Path
import subprocess

from Bio.Seq import Seq


class RISSTError(RuntimeError):
    pass


@dataclass
class Alignment:
    query: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    contig: str
    contig_length: int
    start: int
    end: int
    matches: int
    block_length: int
    score: int

    @classmethod
    def from_paf(cls, line: str):
        fields = line.split('\t')
        score = next(tag[5:] for tag in fields[12:] if tag.startswith('AS:i:'))
        return cls(fields[0], *map(int, fields[1:4]), fields[4], fields[5],
                   *map(int, fields[6:11]), int(score))

    @property
    def partial(self):
        left = self.query_start
        right = self.query_length - self.query_end
        if self.strand == '-':
            left, right = right, left
        return ((left or right) and (self.start == 0 or self.end == self.contig_length)
                or left > self.start or right > self.contig_length - self.end)


@dataclass
class Assembly:
    name: str
    index: Path
    contigs: dict[str, Seq]
    threads: int

    def run(self, *args, query=None):
        try:
            return subprocess.run(
                ['minimap2', '-t', str(self.threads), *map(str, args)],
                input=query, text=True, encoding='utf-8', errors='replace',
                capture_output=True, check=True,
            ).stdout
        except subprocess.CalledProcessError as error:
            raise RISSTError(f'{self.name}: minimap2 failed: {error.stderr.strip()}') from error

    def build_index(self, fasta: Path):
        self.run('-d', self.index, fasta)

    def map(self, query: str):
        for line in self.run('-c', self.index, '-', query=query).splitlines():
            try:
                alignment = Alignment.from_paf(line)
            except (IndexError, ValueError, StopIteration) as error:
                raise RISSTError(f'{self.name}: invalid minimap2 PAF: {line}') from error
            if alignment.block_length > 0 and alignment.query_length > 0:
                yield alignment

    def sequence(self, contig, start, end, strand):
        sequence = self.contigs[contig][start:end]
        return sequence if strand == '+' else sequence.reverse_complement()
