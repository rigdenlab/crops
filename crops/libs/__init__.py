"""This is CROPS: Cropping and Renumbering Operations for PDB structure and Sequence files"""

from crops.about import __prog__, __description__, __author__, __date__, __version__, __copyright__

from crops.libs import rescodes as rc

def ressymbol(name,pick=None):
    """Conversion from residue 3-letter symbol to 1-letter symbol.

    :param name: Residue symbol (3-letter convention) or Nucleotide symbol (1,2-letter convention).
    :type name: str
    :param pick: If three-letter code yields multiple results, pick this one (if among results), defaults to None
    :type pick: str, optional
    :return: Residue/Nucleotide symbol (1-letter convention).
    :rtype: str

    """
    wholelist={**rc.reslist,**rc.mod_reslist,**rc.multiple_reslist,
               **rc.nuclist,**rc.mod_nuclist,**rc.multiple_nuclist}
    try:
        oneletter=wholelist[name]
        if isinstance(oneletter,list):
            if pick is not None:
                oneletter=pick if pick in oneletter else oneletter[0]
            else:
                oneletter=oneletter[0]
    except:
        oneletter=0

    return oneletter