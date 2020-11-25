# -*- coding: utf-8 -*-

from crops.about import __prog__, __description__, __author__, __date__, __version__

reslist =	{
  "ALA": "A",
  "ARG": "R",
  "ASN": "N",
  "ASP": "D",
  "PHD": "D",
  "CYS": "C",
  "CCC": "C",
  "CME": "C",
  "CAS": "C",
  "OMC": "C",
  "CSD": "C",
  "GLN": "Q",
  "GLU": "E",
  "GLY": "G",
  "TGP": "G",
  "GDP": "G",
  "OMG": "G",
  "HIS": "H",
  "ILE": "I",
  "DIL": "I",
  "LEU": "L",
  "LYS": "K",
  "MSE": "M",
  "MET": "M",
  "FME": "M",
  "MEN": "M",
  "PHE": "F",
  "PRO": "P",
  "SER": "S",
  "THR": "T",
  "TPO": "T",
  "TRP": "W",
  "TYR": "Y",
  "VAL": "V",
  "MVA": "V",
  "SEC": "U",
  "PSU": "U",
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
  "YG": "G",
  "DI": "I"
}

def ressymbol(name):
    """Conversion from residue 3-letter symbol to 1-letter symbol.

    :param name: Residue symbol (3-letter convention) or Nucleotide symbol (1,2-letter convention).
    :type name: str
    :return: Residue/Nucleotide symbol (1-letter convention).
    :rtype: str

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
