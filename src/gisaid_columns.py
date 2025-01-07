#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
gisaid_columns.py - Data classes to hold column information
author: Bill Thompson
license: GPL 3
copyright: 2024-12-30
"""
from dataclasses import dataclass

@dataclass(frozen=True)
class VaryingCol:
    """
    A dataclass for holding info about an aligned column with variation
    """
    pct:    float
    data:   tuple   # tuple so that it is immutable

@dataclass
class MutationCounts:
    """
    A data class for holding mutation data
    """
    align_pos:      int     # positions are 0 indexed
    ref_pos:        int
    ref_nuc:        str
    ref_freq:       float
    alt_nuc:        str
    alt_freq:       float
    codon_pos:      int = -1 
    codon:          str = ""
    aa:             str = ""
    aa_name:        str = ""
    alt_codon:      str = ""
    alt_aa:         str = ""
    alt_aa_name:    str = ""
    
    def to_list(self) -> list:
        """convert self to a list

        Returns
        -------
        list
            A list of elements from self
        """
        return [self.ref_pos + 1, self.align_pos + 1,
                self.ref_freq, self.alt_freq,
                self.ref_nuc, self.codon_pos + 1,
                self.aa, self.aa_name, self.codon,
                self.alt_nuc, self.alt_aa, self.alt_aa_name, self.alt_codon,
                self.aa + str(self.codon_pos + 1) + self.alt_aa,
                'syn' if self.aa == self.alt_aa else 'non_syn']