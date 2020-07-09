# -*- coding: utf-8 -*-

__prog__="CROPS"
__description__="Cropping and Renumbering Operations for PDB structure and Sequence files"
__author__ = "J. Javier Burgos-Mármol"
__date__ = "Jul 2020"
__version__ = "0.3.1"

_reslist =	{
  "ALA": "A",
  "ARG": "R",
  "ASN": "N",
  "ASP": "D",
  "CYS": "C",
  "GLN": "Q",
  "GLU": "E",
  "GLY": "G",
  "HIS": "H",
  "ILE": "I",
  "LEU": "L",
  "LYS": "K",
  "MSE": "M",
  "MET": "M",
  "PHE": "F",
  "PRO": "P",
  "SER": "S",
  "THR": "T",
  "TRP": "W",
  "TYR": "Y",
  "VAL": "V",
  "SEC": "U",
  "PYL": "O",
  "XAA": "X",
  "ASX": "B",
  "GLX": "Z",
  "XLE": "J"
}

_nuclist =	{
  "A":"A",
  "T":"T",
  "C":"C",
  "G":"G",
  "U":"U",
  "R":"R",
  "Y":"Y",
  "K":"K",
  "M":"M",
  "S":"S",
  "W":"W",
  "B":"B",
  "D":"D",
  "H":"H",
  "V":"V",
  "N":"N",
  "I":"I",
  "DA": "A",
  "DT": "T",
  "DC": "C",
  "DG": "G",
  "DI": "I"
}

def ressymbol(name):
    """Conversion from residue 3-letter symbol to 1-letter symbol.
    
    :param name: Residue symbol (3-letter convention) or Nucleotide symbol (1,2-letter convention).
    :type name: str
    :return oneletter: Residue/Nucleotide symbol (1-letter convention).
    :rtype oneletter: str

    """
    wholelist={**_reslist,**_nuclist}
    try:
        oneletter=wholelist[name]
    except:
        oneletter=0

    return oneletter