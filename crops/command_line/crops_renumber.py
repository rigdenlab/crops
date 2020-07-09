# -*- coding: utf-8 -*-

"""==========
This script will renumber a structure file in agreement with the
residue positions in the sequence file corresponding to that structure.

!!! IMPORTANT
=============
If the input sequence and the input structure files are not from the
same source (e.g. RCSB PDB) a source conflict might occur making the
renumbering operation unsuccessful even if the program does not crash.
"""

__prog__="CROPS"
__description__="Cropping and Renumbering Operations for PDB structure and Sequence files"
__author__ = "J. Javier Burgos-Mármol"
__date__ = "Jul 2020"
__version__ = "0.3.1"

import argparse
import os

from crops.core import cio
from crops.core import ops as cop


def main():
    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__+'\n'+__doc__)

    parser.add_argument("input_seqpath",nargs=1, metavar="Sequence filepath",
                        help="Input sequence filepath.")
    parser.add_argument("input_strpath",nargs=1, metavar="Structure filepath",
                        help="Input structure filepath or dir. If a directory is inserted, it will act on all structure files in such directory.")

    #parser.add_argument("-b","--bulk",action='store_true',
    #                    help="This option will interpret sequence file as containing sequences with more than one protein ID and, optionally, more than one pdb files in a pdb path.")
    parser.add_argument("-o","--outdir",nargs=1,metavar="Output Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence.")
    parser.add_argument('--version', action='version', version='%(prog)s '+ __version__)
    args = parser.parse_args()

    inseq=cio.check_path(args.input_seqpath,'file')
    instr=cio.check_path(args.input_strpath)

    if args.outdir is None:
        outdir=cio.check_path(os.path.dirname(inseq),'dir')
    else:
        outdir=cio.check_path(os.path.dirname(args.outdir[0]),'dir')

    infixlbl=".crops.seq"

    seqset=cio.parseseqfile(inseq)
    strset, fileset=cio.parsestrfile(instr)

    for pdbid, structure in strset.items():
        if pdbid in seqset:
            newstructure=cop.renumberpdb(seqset[pdbid],structure)
            outstr=cio.outpath(outdir,subdir=pdbid,filename=pdbid+infixlbl["renumber"]+os.path.splitext(instr)[1],mksubdir=True)
            newstructure.write_pdb(outstr)

    return

if __name__ == "__main__":
    import sys
    #import traceback

    try:
        main()
        sys.exit(0)
    except:
        sys.exit(1)
    #except Exception as e:
    #    if not isinstance(e, SystemExit):
    #        msg = "".join(traceback.format_exception(*sys.exc_info()))
    #        logger.critical(msg)
    #    sys.exit(1)