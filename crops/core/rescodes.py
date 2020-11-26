# -*- coding: utf-8 -*-

from crops.about import __prog__, __description__, __author__, __date__, __version__

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
  "MET": "M",
  "PHE": "F",
  "PRO": "P",
  "SER": "S",
  "THR": "T",
  "TRP": "W",
  "TYR": "Y",
  "VAL": "V",
  "SEC": "U",
  "PYL": "O"
}

mod_reslist =	{
  "PHD": "D",
  "CCC": "C",
  "CME": "C",
  "CAS": "C",
  "OMC": "C",
  "CSD": "C",
  "TGP": "G",
  "GDP": "G",
  "OMG": "G",
  "DIL": "I",
  "MSE": "M",
  "FME": "M",
  "MEN": "M",
  "TPO": "T",
  "MVA": "V",
  "PSU": "U"
}

multiple_reslist =	{
  "XAA": ["X"]+list(reslist.values()),
  "ASX": ["B","D","N"],
  "GLX": ["Z","E","Q"],
  "XLE": ["J","I","L"]
}

nuclist =	{
  "A": "A",
  "T": "T",
  "C": "C",
  "G": "G",
  "U": "U",
  "I": "I",
  "DA": "A",
  "DT": "T",
  "DC": "C",
  "DG": "G",
  "DI": "I"
}


mod_nuclist =	{
  "YG": "G"
}

multiple_nuclist =	{
  "R": ["R","A","G"],
  "Y": ["Y","C","T","U"],
  "K": ["K","G","T","U"],
  "M": ["M","A","C"],
  "S": ["S","C","G"],
  "W": ["W","A","T","U"],
  "B": ["B","C","G","T","U"],
  "D": ["D","A","G","T","U"],
  "H": ["H","A","C","T","U"],
  "V": ["V","A","C","G"],
  "N": ["N","A","C","G","T","U"]
}

def ressymbol(name,pick=None):
    """Conversion from residue 3-letter symbol to 1-letter symbol.

    :param name: Residue symbol (3-letter convention) or Nucleotide symbol (1,2-letter convention).
    :type name: str
    :param pick: If three-letter code yields multiple results, pick this one (if among results), defaults to None
    :type pick: str, optional
    :return: Residue/Nucleotide symbol (1-letter convention).
    :rtype: str

    """    
    wholelist={**reslist,**mod_reslist,**multiple_reslist,
               **nuclist,**mod_nuclist,**multiple_nuclist}
    try:
        oneletter=wholelist[name]
        if isinstance(oneletter,list):
            if pick is not None:
                oneletter=pick if pick in oneletter else oneletter[0]
    except:
        oneletter=0

    return oneletter