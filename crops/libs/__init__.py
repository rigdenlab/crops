"""Libraries and functions related to them."""

from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

from crops.libs import rescodes as rc
import logging


def ressymbol(name, pick=None):
    """Convert residue name from 3-letter formt to 1-letter format.

    :param name: Residue symbol (3-letter convention) or Nucleotide symbol (1,2-letter convention).
    :type name: str
    :param pick: If three-letter code yields multiple results, pick this one (if among results), otherwise the most standard is returned, defaults to None.
    :type pick: str, optional

    :return: Residue/Nucleotide symbol (1-letter convention).
    :rtype: str

    """
    if len(name) > 3:
        logging.warning('Chemical ' + name + ' not a residue. '
                        'Assigning 1 to it (ligand).')
        oneletter = 1
        return

    wholelist = {**rc.reslist, **rc.mod_reslist, **rc.multiple_reslist,
                 **rc.nuclist, **rc.mod_nuclist, **rc.multiple_nuclist}
    try:
        oneletter = wholelist[name]
        if isinstance(oneletter, list) is True:
            if pick is not None:
                oneletter = pick if pick in oneletter else oneletter[0]
            else:
                oneletter = oneletter[0]
    except Exception:
        logging.warning('Residue ' + name + ' not found in library. '
                        'Assigning 0 to it (unknown).')
        oneletter = 0

    return oneletter
