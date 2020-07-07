# -*- coding: utf-8 -*-

__prog__="CROPS"
__description__="Cropping and Renumbering Operations for PDB structure and Sequence files"
__author__ = "J. Javier Burgos-Mármol"
__date__ = "May 2020"
__version__ = "0.3.0"

reslist =	{
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

nuclist =	{
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
    """
    Conversion from residue 3-letter symbol to 1-letter symbol and nucleotide recognition

    Parameters
    ----------
    name : str
        Residue symbol (3-letter convention) or Nucleotide symbol (1-letter convention)

    Returns
    -------
    oneletter : str
        Residue/Nucleotide symbol (1-letter convention)

    """
    wholelist={**reslist,**nuclist}
    try:
        oneletter=wholelist[name]
    except:
        oneletter=0

    return oneletter

def fastachains(fasta_seq):
    outdict={}
    for i in range(len(fasta_seq)):
        namechar=False
        nameseq=''
        for j in range(len(fasta_seq[i].id)):
            if fasta_seq[i].id[j]==":":
                namechar=True
            elif fasta_seq[i].id[j]=="|":
                if namechar:
                    break
            else:
                if namechar:
                    nameseq += fasta_seq[i].id[j]
        outdict[nameseq]=i
    return outdict
