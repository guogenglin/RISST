# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 17:04:05 2023

@author: 86182
"""

import argparse
import sys
import pathlib
import multiprocessing
import subprocess
import time
import shutil
import re
import random
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord

__version__ = '2.1'

def get_argument():
    # Parsers
    parser = argparse.ArgumentParser(description = 'RISST', formatter_class = argparse.ArgumentDefaultsHelpFormatter)
    parser_group_1 = parser.add_argument_group('Input and Output')
    parser_group_2 = parser.add_argument_group('Parameters')

    # Input and output
    parser_group_1.add_argument('-i', '--input', required = True, nargs = '+', type = str, 
                                help = 'Input FASTA file')
    parser_group_1.add_argument('-r', '--reference', required = True, type = str, 
                                help = 'Reference cps locus file')
    parser_group_1.add_argument('-o', '--output', required = False, type = str, default = 'RISST_output.txt',
                              help = 'Output file')
    # Parameters
    parser_group_2.add_argument('-t', '--threads', required = False, type = int, default = min(multiprocessing.cpu_count(), 4), 
                        help = 'Threads to use for BLAST searches')
    parser_group_2.add_argument('--minimum_piece', type = int, required = False, default = 200,
                        help='minimum cps match in input sequence, fragments lower than this threshold will be ignore')
    parser_group_2.add_argument('--min_gene_cov', required = False, type = float, default = 80.0, 
                               help = 'Minimum percentage coverage to consider a single gene complete. [default: 80.0%]')
    parser_group_2.add_argument('--min_gene_id', required = False, type = float, default = 70.0, 
                               help = 'Minimum percentage identity to consider a single gene complete. [default: 70.0%]')
    parser_group_2.add_argument('--no_cps_sequence', action = 'store_true', help = 'Suppress output of cps sequence file')
    parser_group_2.add_argument('-f', '--figure', action = 'store_true', help = 'Export the gene structure map of cps locus in inputfile')
    parser_group_2.add_argument('-v', '--version', action = 'version', version = 'BacSpecies v' + __version__, 
                        help = 'Show version number and exit')
    return parser

def check_dependencies():
    # Checks dependencies are available
    dependencies = ['makeblastdb', 'blastn', 'blastx', 'prodigal']
    for i in dependencies:
        try:
            subprocess.check_call(['which', i], stdout = subprocess.DEVNULL)
        except subprocess.CalledProcessError:
            print('Error: could not find {} tool'.format(i), file=sys.stderr)
            sys.exit(1)

def create_temp_dir(reference):
    tempdict = pathlib.Path(reference).parent / 'tempdict'
    result_cps_dict = pathlib.Path(reference).parent / 'result_cps_dict'
    if pathlib.Path(tempdict).is_dir():
        shutil.rmtree(tempdict)
    pathlib.Path.mkdir(tempdict)
    if not pathlib.Path(result_cps_dict).is_dir():
        pathlib.Path.mkdir(result_cps_dict)
    return pathlib.Path(tempdict).resolve(), pathlib.Path(result_cps_dict).resolve()
    
def parse_inputfile(inputfile):
    # Generate a dict, with contig_name as key and contig_sequence as value
    input_seq = {}
    for contig in SeqIO.parse(inputfile, 'fasta'):
        input_seq[contig.name] = contig.seq
    if not input_seq:
        print('invalid FASTA file: {}'.format(inputfile))
        sys.exit(1)
    return input_seq
            
def parse_genbank(reference, tempdict):
    # Parse genbank file to fasta, and generate a dict, with serotype as key and correspond cps sequence as value
    ref_seq_name = reference.split('.')[0] + '.fasta'
    ref_seq_file = pathlib.Path(tempdict) / ref_seq_name
    path = pathlib.Path(reference).resolve()
    ref_seqs = open(ref_seq_file, 'wt')
    for serotype in SeqIO.parse(path, 'genbank'):
        cps_num = ''
        for feature in serotype.features:
            if feature.type == 'source' and 'note' in feature.qualifiers:
                cps_num = feature.qualifiers['note'][0][9:].strip()
        sequence_for_write = ''
        sequence = str(serotype.seq)
        while len(sequence) > 60:
            sequence_for_write += sequence[:60] + '\n'
            sequence = sequence[60:]
        if sequence:
            sequence_for_write += sequence
            sequence_for_write += '\n'
        ref_seqs.write('>' + cps_num + '\n')
        ref_seqs.write(sequence_for_write)
    ref_seqs.close()
    ref_dict = {}
    for serotype in SeqIO.parse(ref_seq_file, 'fasta'):
        ref_dict[serotype.name] = serotype.seq
    return ref_seq_file, ref_dict

def get_serotype(inputfile, repa, ref_dict, threads):
    # Find the best serotype of inputfile by sequence alignment
    inpa = pathlib.Path(inputfile).resolve()
    blast_hits = run_blast(inpa, repa, threads, 'n')
    blast_result = {}
    for hit in blast_hits:
        if not hit.qseqid in blast_result:
            blast_result[hit.qseqid] = []
        if hit.qend > hit.qstart:
            blast_result[hit.qseqid].append((hit.qstart, hit.qend, hit.pident))
        else:
            blast_result[hit.qseqid].append((hit.qend, hit.qstart, hit.pident))
    simplified_result = simplify(blast_result)
    best_serotype = ''
    best_coverage = 0.0
    best_identity = 0.0
    for key, value in simplified_result.items():
        length = get_total_length(value)
        length_identity = get_total_length_identity(value)
        coverage = 100.0 * length / len(ref_dict[key])
        identity = 100.0 * length_identity / len(ref_dict[key])
        if coverage > best_coverage:
            best_serotype = key
            best_coverage = coverage
            best_identity = identity
        elif coverage == best_coverage and identity > best_identity:
            best_serotype = key
            best_identity = identity
    if not best_serotype:
        best_serotype = 'NA'
    if best_coverage <= 95.0:
        best_serotype += '?'
    if best_identity > 100.00:
        best_identity = 100.00
    return best_serotype, round(best_coverage, 2), round(best_identity, 2)
        
def get_total_length(value):
    # Sum the length of all non-redundant blast hit of one single cps locus sequence
    cover_length = 0
    for hit in value:
        length = hit[1] - hit[0]
        cover_length += length
    return cover_length

def get_total_length_identity(value):
    # Sum the identified length of all non-redundant blast hit of one single cps locus sequence
    length_identity = 0
    for hit in value:
        length_identity += hit[2]
    return length_identity
        
def simplify(blast_result):
    # Blast result may have multiple result and overlap, simplify the result, and convert the percent_identity to length_identity
    simplified_result = {}
    for key, value in blast_result.items():
        simplified_list = []
        ranks = sorted(value, key = lambda v: v[0])
        start = 0
        end = 0
        length_ident = 0
        for rank in ranks:
            if rank[0] <= end:
                if rank[1] > end:
                    end = rank[1]
                    length_ident += (rank[1] - rank[0]) * rank[2] / 100
                else:
                    continue
            else:
                if end == 0:
                    length_ident += (rank[1] - rank[0]) * rank[2] / 100
                    start, end = rank[0], rank[1]
                else:
                    length_ident += (rank[1] - rank[0]) * rank[2] / 100
                    simplified_list.append((start, end, length_ident))
                    start, end, length_ident = rank[0], rank[1], 0
        simplified_list.append((start, end, length_ident))
        simplified_result[key] = simplified_list
    return simplified_result     

def find_subject(inputfile, best_serotype, tempdict, ref_seqs_name, threads, result_cps_dict, input_seq, minimum_piece):
    # Get the cps sequence of inputfile
    inpa = pathlib.Path(inputfile).resolve()
    current_cps = pathlib.Path(tempdict) / 'curr_cps.fasta'
    if best_serotype[-1] =='?':
        best_serotype = best_serotype[:-1]
    with open(current_cps, 'wt') as file:
        for serotype in SeqIO.parse(ref_seqs_name, 'fasta'):
            if serotype.name == best_serotype:
                file.write('>')
                file.write(best_serotype)
                file.write('\n')
                file.write(str(serotype.seq))
    blast_hits = run_blast(inpa, current_cps, threads, 'n')
    blast_hits = sorted(blast_hits, key = lambda x: x.qstart)
    if not blast_hits:
        print('Error: No homologous cps found, please check the species of {}'.format(inputfile))
        sys.exit(1)
    subject_contig_range = find_subject_location(blast_hits, input_seq, minimum_piece)
    subject_sequence_path, sub_cps_length = get_subject_cps_sequence(inputfile, subject_contig_range, result_cps_dict, input_seq)
    return subject_sequence_path, sub_cps_length
    
def find_subject_location(blast_hits, input_seq, minimum_piece):
    # Find the cps location in inputfile, consice the blast results to get whole length of cps sequence
    subject_contig_range = []
    if blast_hits[0].sseqid == blast_hits[-1].sseqid:
        cps_range = [blast_hits[0].qstart, blast_hits[0].qend]
        final_sstart = 0
        final_send = 0
        for hit in blast_hits:
            if hit.length <= minimum_piece:
                continue
            elif hit.qstart > min(cps_range) and hit.qend < max(cps_range):
                continue
            else:
                cps_range.append(hit.qstart)
                cps_range.append(hit.qend)
                final_sstart = hit.sstart
                final_send = hit.send
        subject_contig_range = [[blast_hits[0].sseqid, min(blast_hits[0].sstart, final_sstart), max(blast_hits[0].send, final_send), blast_hits[0].strand]]
    else:
        cps_range = [blast_hits[0].qstart, blast_hits[0].qend]
        sseqid, sstart, send, strand = [], [], [], []
        for hit in blast_hits:
            if hit.length <= minimum_piece:
                continue
            elif hit.qstart > min(cps_range) and hit.qend < max(cps_range):
                continue
            else:
                cps_range.append(hit.qstart)
                cps_range.append(hit.qend)
                sseqid.append(hit.sseqid)
                sstart.append(hit.sstart)
                send.append(hit.send)
                strand.append(hit.strand)
        for i in range(len(sseqid)):
            subject_contig_range.append([sseqid[i], sstart[i], send[i], strand[i]])
    return subject_contig_range
        
def get_subject_cps_sequence(inputfile, subject_contig_range, result_cps_dict, input_seq):
    # Save the matched cps sequence to a file and saved in "result_cps_dict" folder
    cps_name = 'cps_locus_of_' + inputfile
    sub_cps_length = 0
    cps_path = pathlib.Path(result_cps_dict) / cps_name
    with open(cps_path, 'wt') as file:
        for i in subject_contig_range:
            file.write('>')
            file.write(str(i[0]) + '_' + str(i[1]) + '_' + str(i[2]) + '_' + str(i[3]))
            file.write('\n')
            cps_gene = str(input_seq[i[0]][i[1]:i[2]])
            sub_cps_length += len(cps_gene)
            gene_for_write = ''
            while len(cps_gene) > 60:
                gene_for_write += cps_gene[:60] + '\n'
                cps_gene = cps_gene[60:]
            if cps_gene:
                gene_for_write += cps_gene
                gene_for_write += '\n'
            file.write(gene_for_write)
            file.write('\n')
    return pathlib.Path(cps_path).resolve(), sub_cps_length
    
def gene_alignment(tempdict, reference, best_serotype, subject_sequence_path, min_gene_id, min_gene_cov, threads):
    # Screen the gene distribution of cps of inputfile
    gene_expected, other_genes = [], []
    orfs = prodigal(tempdict, subject_sequence_path)
    cur_cps_ref, other_cps_ref = parse_ref_gene(tempdict, reference, best_serotype)
    for orf in orfs:
        with open('temp_blast_inpa.fasta', 'wt') as file:
            file.write('>')
            file.write(str(orf.contig) + '_' + str(orf.start) + '_' + str(orf.end) + '_' + str(orf.strand))
            file.write('\n')
            file.write(str(orf.sequence))
        inpa = pathlib.Path('temp_blast_inpa.fasta').resolve()
        blast_hits = run_blast(cur_cps_ref, inpa, threads, 'x')
        best_match = ''
        best_cov = min_gene_cov
        best_pident = min_gene_id
        if not blast_hits:
            continue
        for hit in blast_hits:
            if hit.length < 100:
                continue
            if hit.pident >= best_pident and hit.x_query_cov >= best_cov:
                best_pident = hit.pident
                best_cov = hit.x_query_cov
                best_match = hit.sseqid
        if best_match and best_match not in gene_expected:
            gene_expected.append(best_match)
        elif not best_match:
            blast_hits = run_blast(other_cps_ref, inpa, threads, 'x')
            if not blast_hits:
                continue
            for hit in blast_hits:
                if hit.length < 100:
                    continue
                if hit.pident >= best_pident and hit.x_query_cov >= best_cov:
                    best_pident = hit.pident
                    best_cov = hit.x_query_cov
                    best_match = hit.sseqid
            if best_match and best_match not in other_genes:
                other_genes.append(best_match)
            else:
                continue
        orf.mapname = best_match
        orf.serotype = best_match.strip('cps').split('_')[0]
    pathlib.Path(inpa).unlink()
    other_genes_counts = len(other_genes)
    if best_serotype[-1] =='?':
        best_serotype = best_serotype[:-1]
    pathlib.Path(cur_cps_ref).unlink()
    pathlib.Path(other_cps_ref).unlink()
    return gene_expected, other_genes, other_genes_counts, orfs

def process_gene_result(gene_expected, other_genes):
    # Link the gene hits to one str for output in table
    expected, others = '', ''
    for i in gene_expected:
        expected += i
        expected += ';'
    for i in other_genes:
        others += i
        others +=';'
    return expected, others
    
def parse_ref_gene(tempdict, reference, best_serotype):
    # Separate the matched serotype reference cps and other reference sequences, prepared for gene screen
    if best_serotype[-1] =='?':
        best_serotype = best_serotype[:-1]
    cur_cps_ref = pathlib.Path(tempdict).resolve() / 'cur_cps_ref.fasta'
    other_cps_ref = pathlib.Path(tempdict).resolve() / 'other_cps_ref.fasta'
    path = pathlib.Path(reference).resolve()
    for serotype in SeqIO.parse(path, 'genbank'):
        cps_num = ''
        gene_count = 1
        gene_seqs = open(cur_cps_ref, 'at')
        for feature in serotype.features:
            if feature.type == 'source' and 'note' in feature.qualifiers:
                gene_seqs.close()
                for note in feature.qualifiers['note']:
                    cps_num = feature.qualifiers['note'][0][9:].strip()
                gene_count = 1
            if cps_num == best_serotype:
                gene_seqs = open(cur_cps_ref, 'at')
            else:
                gene_seqs = open(other_cps_ref, 'at')
            if feature.type == 'CDS' and 'translation' in feature.qualifiers:
                gene_name = 'cps' + cps_num + '_' + str(gene_count).zfill(3)
                gene_count += 1
                gene_for_write = ''
                gene_prot = str(feature.qualifiers['translation']).strip("[]'")
                while len(gene_prot) > 60:
                    gene_for_write += gene_prot[:60] + '\n'
                    gene_prot = gene_prot[60:]
                if gene_prot:
                    gene_for_write += gene_prot
                    gene_for_write += '\n'
                gene_seqs.write('>' + gene_name + '\n')
                gene_seqs.write(gene_for_write)
    gene_seqs.close()
    return cur_cps_ref, other_cps_ref
      
def prodigal(tempdict, subject_sequence_path):
    # Generate a list with the informations and sequence of each orf in the inputfile
    command = 'prodigal -f sco -i {} -t {} -m'.format(subject_sequence_path, 'risst.tr')
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    raw_orfs = result.stdout.decode().replace('\r', '')
    orfs = group_prodigal(raw_orfs)
    input_seq = parse_inputfile(subject_sequence_path)
    for orf in orfs:
        orf.sequence = input_seq[orf.contig][orf.start-1:orf.end]
    return orfs

def group_prodigal(raw_orfs):
    # group every orf to a list
    orfs = []
    contig = ''
    RESULT_RE = re.compile(r'^>[0-9]+_([0-9]+)_([0-9]+)_([-+])$')
    CONTIG_RE = re.compile(r'^# Sequence.+?seqhdr="(.+?)"(?:;|$)')
    for line in raw_orfs.rstrip().split('\n'):
        if line.startswith('# Sequence Data'):
            name = CONTIG_RE.match(line).group(1)
            contig = name.split(' ')[0]
        elif line.startswith('# Model Data'):
            continue
        else:
            result = RESULT_RE.match(line).groups()
            orfs.append(Orf(contig, *result))
    return orfs

class Orf(object):
    # Parse the prodigal results
    def __init__(self, contig, start, end, strand):
        self.contig = contig
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.sequence = str()
        self.mapname = ''
        self.serotype = ''
    
def run_blast(inpa, repa, threads, setting):
    # Do blast, iterator the result to a list
    blast_hits = []
    if setting == 'n':
        command = ['blastn', '-query', repa, '-subject', inpa, '-num_threads', str(threads), '-outfmt', 
                   '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    else:
        command = ['blastx', '-query', repa, '-subject', inpa, '-outfmt', 
                   '6 qseqid sseqid qstart qend sstart send evalue bitscore length pident qlen qseq']
    process = subprocess.run(command, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out = process.stdout.decode()
    for line in line_iterator(out):
        blast_hits.append(BlastResult(line))
    return blast_hits

def line_iterator(line_breaks):
    # Handle the BLAST output and remove the line breaks 
    line = -1
    while True:
        nextline = line_breaks.find('\n', line + 1)
        if nextline < 0:
            break
        yield line_breaks[line + 1:nextline]
        line = nextline

class BlastResult(object):
    # Handle the BLAST output
    def __init__(self, hit_string):
        parts = hit_string.split('\t')
        self.qseqid = parts[0].split('_')[-1]
        self.sqseqid = parts[0]
        self.sseqid = parts[1]
        self.qstart = int(parts[2]) - 1
        self.qend = int(parts[3])
        self.sstart = int(parts[4]) - 1
        self.send = int(parts[5])
        if self.sstart <= self.send:
            self.strand = '+'
        else:
            self.sstart, self.send = self.send, self.sstart
            self.strand = '-'
        self.length = int(parts[8])
        self.pident = float(parts[9])
        self.query_cov = 100.0 * len(parts[11]) / float(parts[10])
        self.x_query_cov = 300.0 * len(parts[11]) / float(parts[10])

def draw_gene_map(inputfile, best_serotype, orfs, sub_cps_length):
    # Export a gene structure map
    map_dict = pathlib.Path(inputfile).parent / 'gene_structure_map'
    if not pathlib.Path(map_dict).is_dir():
        pathlib.Path.mkdir(map_dict)
    inputfile = strip_suffix(inputfile)
    map_name = 'gene_structure_map_of_' + inputfile + '.png'
    map_path = pathlib.Path(map_dict) / map_name
    if best_serotype[-1] =='?':
        best_serotype = best_serotype[:-1]
    soft_color_database = ['#99CCCC', '#FFFFCC', '#FFCC99', '#CCFFFF', '#CCCCFF', '#FFCCCC', '#CCFFCC']
    light_color_database = ['#FF0033', '#CCFF00', '#0099CC', '#009966', '#CC3399']
    features = []
    orf_contig = ''
    orf_length_all = 0
    orf_length_for_one_contig = 0
    direction = orfs[0].strand
    for orf in orfs:
        orf_color = ''
        if best_serotype == orf.serotype:
            orf_color = random.choice(soft_color_database)
        else:
            orf_color = random.choice(light_color_database)
        if orf_contig and orf.contig != orf_contig:
            orf_length_all += orf_length_for_one_contig
        orf_length_for_one_contig = orf.end
        orf_contig = orf.contig
        orf.start += orf_length_all
        orf.end += orf_length_all
        if direction == '+':
            if orf.strand == '+':
                orf_strand = +1
            else:
                orf_strand = -1
        elif direction == '-':
            if orf.strand == '+':
                orf_strand = -1
            else:
                orf_strand = +1
        features.append(GraphicFeature(start = orf.start, end = orf.end, strand = orf_strand, color = orf_color, label = orf.mapname))
    record = GraphicRecord(sequence_length = sub_cps_length, features = features)
    ax, _ = record.plot(figure_width = 20)
    ax.figure.savefig(map_path, dpi = 600)

def strip_suffix(inputfile):
    # Strip the suffix of inputfile
    if inputfile.lower().endswith('.fa'):
        inputfile = inputfile[:-3]
    elif inputfile.lower().endswith('.fna'):
        inputfile = inputfile[:-4]
    elif inputfile.lower().endswith('.fas'):
        inputfile = inputfile[:-4]
    elif inputfile.lower().endswith('.fasta'):
        inputfile = inputfile[:-6]
    return inputfile

def generate_output(output):
    # Generate a blank output table file
    if pathlib.Path(output).is_file():
        return
    headers = ['Isolate', 'Serotype', 'Coverage', 'Identity', 'Identified_gene_details', 'Other_gene_counts', 'Other_gene_details']
    with open(output, 'wt') as file:
        file.write('\t'.join(headers))
        file.write('\n')

def output(output, inputfile, serotype, sero_coverage, sero_identity, expected_genes, other_genes_counts, genes_from_other_cps):
    # Generate output
    simple_output = inputfile + ' : '  + ' ' + 'serotype' + ' ' + serotype
    line = [inputfile, serotype, str(sero_coverage), str(sero_identity), expected_genes, str(other_genes_counts), genes_from_other_cps]
    print(simple_output)
    with open(output, 'at') as file:
        file.write('\t'.join(line))
        file.write('\n')

def main():
    starttime = time.perf_counter()
    # Initialize
    args = get_argument().parse_args()
    check_dependencies()
    tempdict, result_cps_dict = create_temp_dir(args.reference)
    ref_seqs_name, ref_dict = parse_genbank(args.reference, tempdict)
    # Run this pipeline for each single input genome
    for inputfile in args.input:
        input_seq = parse_inputfile(inputfile)
        serotype, sero_coverage, sero_identity = get_serotype(inputfile, ref_seqs_name, ref_dict, args.threads)
        subject_sequence_path, sub_cps_length = find_subject(inputfile, serotype, tempdict, ref_seqs_name, args.threads, result_cps_dict, input_seq, args.minimum_piece)
        gene_expected, other_genes, other_genes_counts, orfs = gene_alignment(tempdict, args.reference, serotype, subject_sequence_path, args.min_gene_id, args.min_gene_cov, args.threads)
        expected_genes, genes_from_other_cps = process_gene_result(gene_expected, other_genes)
    # Generate output
        if args.figure:
            draw_gene_map(inputfile, serotype, orfs, sub_cps_length)
        generate_output(args.output)
        output(args.output, inputfile, serotype, sero_coverage, sero_identity, expected_genes, other_genes_counts, genes_from_other_cps)
    shutil.rmtree(tempdict)
    if args.no_cps_sequence:
        shutil.rmtree(result_cps_dict)
    # Total time count
    endtime = time.perf_counter() - starttime
    per_genome_time = endtime / len(args.input)
    print('{:.1f}h{:.1f}m{:.1f}s for one genome'.format(per_genome_time // 3600, per_genome_time % 3600 // 60, per_genome_time % 60))
    print('Total time consumed : {:.1f}h{:.1f}m{:.1f}s'.format(endtime // 3600, endtime % 3600 // 60, endtime % 60))
   
main()
