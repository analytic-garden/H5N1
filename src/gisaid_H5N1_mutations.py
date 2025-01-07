#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gisaid_H%N1_mutatians.py - analyze H5N1 mutations
author: Bill Thompson
license: GPL 3
copyright: 2024-12-26
"""
import sys
import argparse
from Bio import AlignIO
from Bio import Align
from Bio import SeqRecord
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from collections import Counter
import pandas as pd
from gisaid_columns import VaryingCol, MutationCounts

def GetArgs() -> argparse.Namespace:
    """Get command line arguments

    Returns
    -------
    argparse.Namespace
        The command line arguments
    """
    def ParseArgs(parser):
        class Parser(argparse.ArgumentParser):
            def error(self, message):
                sys.stderr.write('error: %s\n' % message)
                self.print_help()
                sys.exit(2)

        parser = Parser(description='Analyze H5N1 mutationss.')
        parser.add_argument('-i', '--msa_file',
                            required = True,
                            help = 'MSA FASTA file created by mafft. Refereence sequence first in file. Heaaders must be ISOLATE_ID',
                            type = str)
        parser.add_argument('-d', '--H5N1_CSV',
                            required = True,
                            help = 'A CSV from a dataframe created by the R function fasta2dataframe.',
                            type = str)
        parser.add_argument('-g', '--GenBank_Ref',
                            required = True,
                            help = 'GenBank file of reference sequence',
                            type = str)
        parser.add_argument('-o', '--output',
                            required = True,
                            help = 'Output CSV file',
                            type = str)

        return parser.parse_args()

    parser = argparse.ArgumentParser()
    args = ParseArgs(parser)
    
    return args

def read_alignment_file(filename: str) -> tuple[dict[str, SeqRecord.SeqRecord], Align.MultipleSeqAlignment]:
    """ Read an MSA FASTA file

    Parameters
    ----------
    filename : str
        MSA file

    Returns
    -------
    tuple[dict[str, SeqRecord.SeqRecord], Align.MultipleSeqAlignment]
        A tuple:
            a dictionary with record ID as key and SeqRecords as value
            a BioAlign Multiple Sequence Alignment
    """
    alignment = AlignIO.read(open(filename), "fasta")

    align_dict = dict()
    for record in alignment:
        align_dict[record.id] = record
        
    return (align_dict, alignment)

def ref_pos_to_alignment(align: dict[str, SeqRecord.SeqRecord], ref_id: str) -> dict[int, int]:
    """Map aligned positions to reference sequence

    Parameters
    ----------
    align : Dictionary returned by read_alignment_file
        A BioPython MSA
    ref_id : str
        ID of reference seq

    Returns
    -------
    dict[int, int]
        A dictionary
            key: alignment position
            value: position in reference sequence, or -1 id a gap position
    """
    ref_seq = align[ref_id].seq

    pos_map = dict()
    ref_pos = 0    # 0 based
    for pos in range(len(ref_seq)):
        if ref_seq[pos] in ['A', 'C', 'G', 'T']:
            pos_map[pos] = ref_pos
            ref_pos += 1
        else:
            pos_map[pos] = -1

    return pos_map

def get_varying_columns(align: Align.MultipleSeqAlignment,
                        consensus_cutoff: float = 1.0,
                        start: int = 0, end: int = -1) -> dict[int, VaryingCol]:
    """Iterate across columns and find columns in alignment that vary more than consensus_cutoff

    Parameters
    ----------
    align : Align.MultipleSeqAlignment
        A BioPython MSA
    consensus_cutoff : float, optional
        Cutoff for consensus. Only columns with consensus freq less that this value are accepted, by default 1.0
    start : int, optional
        Start searching at this columns, by default 0
    end : int, optional
        last column to seach. If -1, end is last column. , by default -1

    Returns
    -------
    dict[int, VaryingCol]
        A dictionary
            key: column cposition
            value: a VaryingCol object
    """
    if end == -1:
        end = align.get_alignment_length()
        
    variant_cols = dict()
    for col in range(start, end):
        c = Counter(align[:, col])
        common = c.most_common(1)
        if common[0][0] in ['A', 'C', 'G', 'T']:
            denom = sum([c[k] for k in ['A', 'C', 'G', 'T', '-']])
            freq = common[0][1] / denom
            if freq < consensus_cutoff:
                # variant_cols[col] = (freq, tuple(list(align[:, col])))
                variant_cols[col] = VaryingCol(freq, tuple(list(align[:, col])))

    return variant_cols

def count_mutations(variant_cols: dict[int, tuple[float, list[float, list[str]]]],
                    pos_map: dict[int, int],
                    align_length: int) -> dict[int, MutationCounts]:
    """Count the common nucleotides and alternates

    Parameters
    ----------
    variant_cols : dict[int, tuple[float, list[float, list[str]]]]
        A dictionary returned by et_varying_column
    pos_map : dict[int, int]
        A dictionary returned by ref_pos_to_alignment
    align_length : int
        The number of rows in the alignment

    Returns
    -------
    dict[str, dict[str, int | float | str]]
        A dictionary
            key: reference position
            value: a MutationCounts object
    """
    mutations = dict()
    for pos, rec in variant_cols.items():
        ref_pos = pos_map[pos]
        nucleotides = rec.data
        c = Counter(nucleotides)
        ref_nuc = c.most_common(2)[0][0]
        ref_freq = c.most_common(2)[0][1] / align_length
        alt_nuc = c.most_common(2)[1][0]
        alt_freq = c.most_common(2)[1][1] / align_length
        mutations[ref_pos] = MutationCounts(align_pos = pos, ref_pos = ref_pos,
                                            ref_nuc = ref_nuc, ref_freq = ref_freq, 
                                            alt_nuc = alt_nuc, alt_freq = alt_freq)
        
    return mutations

def get_codons(gb_ref: SeqRecord, mutations:  dict[int, MutationCounts]) -> None:
    """get codons and alternate codons for mutated positions

    Parameters
    ----------
    gb_ref : SeqRecord
        Genbank record of reference sequnec
    mutations : dict[int, MutationCounts]
        a dictionary of MutationCounts for mutatated columns returned by count_mutations
        key: reference seq column
        value: a MutationCounts object
        
    Note: mutations is altered by this function
    """
    seq = list(gb_ref.seq)
    for pos in range(len(seq)):
        if pos in mutations.keys():
            mutations[pos].codon_pos = int(pos / 3) * 3
            
            mutations[pos].codon = "".join(seq[mutations[pos].codon_pos:mutations[pos].codon_pos+3])
            mutations[pos].aa = str(Seq(mutations[pos].codon).translate())
            mutations[pos].aa_name = seq3(mutations[pos].aa)
            
            pos_in_codon = pos % 3
            s = list(mutations[pos].codon)
            s[pos_in_codon] = mutations[pos].alt_nuc
            mutations[pos].alt_codon = "".join(s)
            mutations[pos].alt_aa = str(Seq(mutations[pos].alt_codon).translate())
            mutations[pos].alt_aa_name = seq3(mutations[pos].alt_aa)
     
def mutation_df(mutations:  dict[int, MutationCounts]) -> pd.DataFrame:
    """Make a DataFrame from mutations dictionary

    Parameters
    ----------
    mutations : dict[int, MutationCounts]
        A dictionary from get_codons

    Returns
    -------
    pd.DataFrame
        A dataframe of mutations and their positions
    """
    df = pd.DataFrame(columns = ['reference_position', 'alignment_pos',
                                 'ref_freq', 'alt_freq',
                                 'ref_nucleotide', 'codon_pos',
                                 'ref_aa', 'aa_name', 'ref_codon',
                                 'alt_nucleotide', 'alt_aa', 'alt_aa_name', 'alt_codon',
                                 'mutation', 'synonomous'])
    for pos, mut in mutations.items():
        df = pd.concat([df, pd.DataFrame([mut.to_list()], columns = df.columns)], ignore_index = True)
        
    return df
            
def main():
    args = GetArgs()
    
    # read the alinment and GenBank referenc file
    align_dict, alignment = read_alignment_file(args.msa_file)
    ref_id = alignment[0].id
    
    with open(args.GenBank_Ref, "r") as f:
        gb_ref = SeqIO.read(f,'genbank')
        
    # get columns with changes
    pos_map = ref_pos_to_alignment(align_dict, ref_id)
    variant_cols = get_varying_columns(alignment, consensus_cutoff=0.9)
    
    mutations = count_mutations(variant_cols, pos_map, len(alignment))
    get_codons(gb_ref, mutations)
    
    df = mutation_df(mutations)
    df.to_csv(args.output, index = False)
    
    print()

if __name__ == "__main__":
    main()
