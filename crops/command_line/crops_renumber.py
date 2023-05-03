"""This script will renumber a structure file in agreement with the residue positions in the sequence file corresponding to that structure.

Non-polymer elements are numbered starting right after the final (TER) residue.

WARNING: If the input sequence and the input structure files are not from the
same source (e.g. RCSB PDB) a source conflict might occur making the
renumbering operation unsuccessful even if the program does not crash.
"""

from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

import argparse
import os

from crops.iomod import check_path
from crops.iomod import outpathgen
from crops.iomod import parsers as cin
from crops.core import ops as cop
from crops import command_line as ccl

logger = None


def create_argument_parser():
    """Create a parser for the command line arguments used in crops-renumber."""
    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__+os.linesep+__doc__)

    parser.add_argument("input_seqpath", nargs=1, metavar="Sequence_filepath",
                        help="Input sequence filepath.")
    parser.add_argument("input_strpath", nargs=1, metavar="Structure_filepath",
                        help="Input structure filepath or dir. If a directory is inserted, it will act on all structure files in such directory.")

    parser.add_argument("-f", "--force_alignment", action='store_true', default=False,
                        help="Use Needleman-Wunsch algorithm to try to bypass small disagreements between fasta and pdb sequences.")

    parser.add_argument("-p", "--preselect", nargs='+', metavar="Oligoseq_ids", default=None,
                        help="From all the sequences in the input sequence file, just print out this preselected subset.")

    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence.")
    parser.add_argument('--version', action='version', version='%(prog)s '+ __version__)

    return parser


def main():
    """Renumber a structure file in agreement with the residue positions in the sequence file corresponding to that structure.

    Non-polymer elements are numbered starting right after the final (TER) residue.
    IMPORTANT: If the input sequence and the input structure files are not from the
    same source (e.g. RCSB PDB) a source conflict might occur making the
    renumbering operation unsuccessful even if the program does not crash.

    """
    # INITIALISE AND PARSE ARGUMENTS FROM COMMAND LINE
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = ccl.crops_logger(level="info")
    logger.info(ccl._welcome())

    inseq = check_path(args.input_seqpath[0], 'file')
    instr = check_path(args.input_strpath[0])

    if args.outdir is None:
        outdir = check_path(os.path.dirname(inseq), 'dir')
    else:
        outdir = check_path(os.path.join(args.outdir[0], ''), 'dir')
    infixlbl = ".crops.seq"

    # PARSE INPUT FILES
    logger.info('Parsing sequence file '+inseq)
    if args.preselect is not None:
        subset = set(args.preselect)
    else:
        subset = None
    seqset = cin.parseseqfile(seq_input=inseq, inset=subset)
    logger.info('Done')

    logger.info('Parsing structure file '+instr)
    strset, fileset = cin.parsestrfile(instr)
    logger.info('Done')

    # MAIN OPERATION / PRINT OUT RESULTS WITHIN
    logger.info('Renumbering structure(s)...')
    for pdbid, structure in strset.items():
        found = False
        for seqname in seqset:
            if ((seqname in pdbid) or
                    (len(seqset) == 1 and len(strset) == 1)):
                finalid = seqname
                try:
                    newstructure = cop.renumber_pdb(seqset[seqname], structure)
                except (AttributeError, IndexError, ValueError) as e:
                    logger.warning('Something has gone wrong during renumbering:\n{}'.format(e))
                    if args.force_alignment:
                        logger.info('Attempting Needleman-Wunsch...')
                        newstructure = cop.renumber_pdb_needleman(seqset[seqname], structure)
                    else:
                        logger.critical('Unable to renumber the structure, exiting now. '
                                        'Try again with -f option to force the alignment.')
                        return

                fout = finalid+infixlbl+os.path.splitext(instr)[1]
                outstr = outpathgen(outdir, subdir=finalid,
                                    filename=fout, mksubdir=True)
                newstructure.write_minimal_pdb(outstr)
                found = True
        if found is False:
            logger.warning("Identifier '"+pdbid+"' not found in sequence input.")

    # FINISH
    logger.info('Done' + os.linesep)

    return


if __name__ == "__main__":
    import sys
    import traceback

    try:
        main()
        logger.info(ccl._ok())
        sys.exit(0)
    except Exception as e:
        if not isinstance(e, SystemExit):
            msg = "".join(traceback.format_exception(*sys.exc_info()))
            logger.critical(msg)
        sys.exit(1)
