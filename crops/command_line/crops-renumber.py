"""==========
This script will renumber a structure file in agreement with the
residue positions in the sequence file corresponding to that structure.
Non-polymer elements are numbered starting right after the final (TER) residue.
IMPORTANT: If the input sequence and the input structure files are not from the
same source (e.g. RCSB PDB) a source conflict might occur making the
renumbering operation unsuccessful even if the program does not crash.
"""

from crops.about import __prog__, __description__, __author__, __date__, __version__

import argparse
import os

from crops.io import check_path
from crops.io import outpathgen
from crops.io import parsers as cin
from crops.core import ops as cop
from crops import command_line as ccl

logger=None

def create_argument_parser():
    """Create a parser for the command line arguments used in crops-renumber"""

    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__+'\n'+__doc__)

    parser.add_argument("input_seqpath",nargs=1, metavar="Sequence_filepath",
                        help="Input sequence filepath.")
    parser.add_argument("input_strpath",nargs=1, metavar="Structure_filepath",
                        help="Input structure filepath or dir. If a directory is inserted, it will act on all structure files in such directory.")

    parser.add_argument("-o","--outdir",nargs=1,metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence.")
    parser.add_argument('--version', action='version', version='%(prog)s '+ __version__)

    return parser

def main():

    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = ccl.crops_logger(level="info")
    logger.info(ccl.welcome())

    inseq=check_path(args.input_seqpath[0],'file')
    instr=check_path(args.input_strpath[0])

    if args.outdir is None:
        outdir=check_path(os.path.dirname(inseq),'dir')
    else:
        outdir=check_path(os.path.join(args.outdir[0],''),'dir')
    infixlbl=".crops.seq"

    logger.info('Parsing sequence file '+inseq)
    seqset=cin.parseseqfile(inseq)
    logger.info('Done')

    logger.info('Parsing structure file '+instr)
    strset, fileset=cin.parsestrfile(instr)
    logger.info('Done')

    logger.info('Renumbering structure(s)...')
    for pdbid, structure in strset.items():
        if pdbid in seqset:
            newstructure=cop.renumber_pdb(seqset[pdbid],structure)
            outstr=outpathgen(outdir,subdir=pdbid,filename=pdbid+infixlbl+os.path.splitext(instr)[1],mksubdir=True)
            #newstructure.write_pdb(outstr)
            newstructure.write_minimal_pdb(outstr)
    logger.info('Done\n')

    return

if __name__ == "__main__":
    import sys
    import traceback

    try:
        main()
        logger.info(ccl.ok())
        sys.exit(0)
    except Exception as e:
        if not isinstance(e, SystemExit):
            msg = "".join(traceback.format_exception(*sys.exc_info()))
            logger.critical(msg)
        sys.exit(1)
