"""==========
This script will take a sequence (fasta) file as an input and produce
several different fasta files, one per sequence/chain.
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
    """Create a parser for the command line arguments used in crops-splitseqs"""

    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__+os.linesep+__doc__)

    parser.add_argument("input_seqpath", nargs=1, metavar="Sequence_filepath",
                        help="Input sequence filepath.")

    parser.add_argument("-i", "--individual", action='store_true', default=False,
                          help="One separated output fasta file per each sequence.")

    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence.")
    parser.add_argument('--version', action='version', version='%(prog)s '+ __version__)

    return parser

def main():
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = ccl.crops_logger(level="info")
    logger.info(ccl.welcome())

    inseq = check_path(args.input_seqpath[0], 'file')

    if args.outdir is None:
        outdir = check_path(os.path.dirname(inseq), 'dir')
    else:
        outdir = check_path(os.path.join(args.outdir[0], ''), 'dir')

    logger.info('Parsing sequence file '+inseq)
    seqset = cin.parseseqfile(inseq)
    logger.info('Done')

    for key, S in seqset.items():
        for key2, monomer in S.imer.items():
            if args.individual is True:
                fout = (key + '_' + key2 +
                        os.path.splitext(os.path.basename(inseq))[1])
                outseq = outpathgen(outdir, filename=fout)
                monomer.dump(outseq)

    logger.info('Done' + os.linesep)

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