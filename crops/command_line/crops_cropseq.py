"""This script will remove a number of residues from a sequence file in agreement to the intervals and other details supplied."""

from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

import argparse
import os

from crops.iomod import check_path
from crops.iomod import outpathgen
from crops.iomod import parsers as cin
from crops.iomod import taggers as ctg
from crops.core import ops as cop
from crops import command_line as ccl

import time
import copy

logger = None


def create_argument_parser():
    """Create a parser for the command line arguments used in crops-cropseq."""
    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__+os.linesep+__doc__)
    parser.add_argument("input_seqpath", nargs=1, metavar="Sequence_filepath",
                        help="Input sequence filepath.")
    parser.add_argument("input_database", nargs=1, metavar="Intervals_database",
                        help="Input intervals database filepath.")

    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence.")

    parser.add_argument("-p", "--preselect", nargs='+', metavar="Oligoseq_ids", default=None,
                        help="From all the sequences in the input sequence file, just print out this preselected subset.")

    outfiles = parser.add_mutually_exclusive_group(required=False)
    outfiles.add_argument("-s", "--sort", nargs=1, metavar="Sort_type",
                          help="Sort output sequences in descending order by criteria provided - 'ncrops' or 'percent'. Add 'IN' ('ncropsIN', 'percentIN') to ignore numbers from terminals. Only for multiple ID fasta inputs.")
    outfiles.add_argument("-i", "--individual", action='store_true', default=False,
                          help="One separated fasta file per each sequence.")

    sections = parser.add_mutually_exclusive_group(required=False)
    sections.add_argument("-t", "--terminals", action='store_true', default=False,
                          help="Ignore interval discontinuities and only crop the ends off.")
    sections.add_argument("-u", "--uniprot_threshold", nargs=2, metavar=("Uniprot_ratio_threshold", "Sequence_database"),
                          help='Act if SIFTS database is used as intervals source AND %% residues from single Uniprot sequence is above threshold. Threshold: [MIN,MAX)=[0,100). Database path: uniclust##_yyyy_mm_consensus.fasta-path or server-only. The latter requires internet connexion.')
    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser


def main():
    """Remove a number of residues from a sequence file in agreement to the intervals and other details supplied.

    :raises ValueError: For wrong argument values.

    """
    # INITIALISE AND PARSE ARGUMENTS FROM COMMAND LINE
    starttime = time.time()
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = ccl.crops_logger(level="info")
    logger.info(ccl._welcome())

    inseq = check_path(args.input_seqpath[0], 'file')
    indb = check_path(args.input_database[0], 'file')
    if args.uniprot_threshold is not None:
        if args.uniprot_threshold[1] != 'server-only':
            insprot = check_path(args.uniprot_threshold[1])
        else:
            insprot = 'server-only'
    else:
        insprot = None

    if args.uniprot_threshold is not None:
        minlen = float(args.uniprot_threshold[0])
        if minlen < 0.0 or minlen > 100.0:
            logger.critical('The UniProt threshold is a percentage and, therefore, it must fulfil 0 < threshold < 100.')
            raise ValueError
    else:
        minlen = 0.0
    targetlbl = ctg.target_format(indb, terms=args.terminals, th=minlen)
    infixlbl = ctg.infix_gen(indb, terms=args.terminals)

    if args.outdir is None:
        outdir = check_path(os.path.dirname(inseq), 'dir')
    else:
        outdir = check_path(os.path.join(args.outdir[0], ''), 'dir')

    if args.sort is not None:
        if (args.sort[0].lower() != 'ncrops' and
                args.sort[0].lower() != 'percent' and
                args.sort[0].lower() != 'ncropsin' and
                args.sort[0].lower() != 'percentin'):
            logger.critical("Arguments for sorting option can only be "
                            "either 'ncrops' or 'percent'.")
            raise ValueError
        else:
            sorter = args.sort[0].lower()

    # PARSE INPUT FILES
    logger.info('Parsing sequence file ' + inseq)
    if args.preselect is not None:
        subset = set(args.preselect)
    else:
        subset = None
    seqset = cin.parseseqfile(seq_input=inseq, inset=subset)
    logger.info('Done')

    logger.info('Parsing interval database file '+indb)
    if len(seqset) > 0:
        intervals = cin.import_db(indb, pdb_in=seqset)
    else:
        logger.critical('No chains were imported from sequence file.')
        raise ValueError
    logger.info('Done'+os.linesep)

    if insprot is not None and minlen > 0.0:
        logger.info('Parsing uniprot sequence file: '+insprot)
        uniprotset = set()
        for seqncid, seqnc in seqset.items():
            chains = seqnc.chainlist()
            for monomerid in chains:
                monomer = seqnc.imer[seqnc.whatseq(monomerid)]
                if 'uniprot' in intervals[seqncid][monomerid].tags:
                    for key in intervals[seqncid][monomerid].tags['uniprot']:
                        if key.upper() not in uniprotset:
                            uniprotset.add(key.upper())

        upserver = True if insprot == 'server-only' else False
        uniprotset = cin.parseseqfile(insprot, inset=uniprotset, use_UPserver=upserver)
        logger.info('Done'+os.linesep)

    # MAIN OPERATION
    logger.info('Cropping sequence(s)...')
    if len(seqset) > 1 and args.sort is not None:
        sorted_outseq = {}

    for key, S in seqset.items():
        for key2, monomer in S.imer.items():
            if key in intervals:
                if any(x in intervals[key] for x in S.imer[key2].chains):
                    for x in S.imer[key2].chains:
                        if (x in intervals[key]):
                            key3 = x
                            break
                    if insprot is not None and minlen > 0.0:
                        newinterval = intervals[key][key3].deepcopy()
                        newinterval.tags['description'] += ' - Uniprot threshold'
                        newinterval.subint = []
                        unilbl = ' uniprot chains included: '
                        for unicode, uniintervals in intervals[key][key3].tags['uniprot'].items():
                            val = 100*uniintervals.n_elements()
                            val = val / uniprotset[unicode].imer['1'].length()
                            if val >= minlen:
                                newinterval = newinterval.union(intervals[key][key3].intersection(uniintervals))
                                unilbl += unicode
                        monomer = cop.crop_seq(monomer, newinterval,
                                               targetlbl+unilbl,
                                               terms=args.terminals)
                    else:
                        monomer = cop.crop_seq(monomer, intervals[key][key3],
                                               targetlbl, terms=args.terminals)
                else:
                    monomer.infostring += (targetlbl + ' (reference not found in database)')
                    monomer.cropmap = {}
                    for n in range(1, monomer.length()+1):
                        monomer.cropmap[n] = n
                        monomer.cropbackmap = copy.deepcopy(monomer.cropmap)

            monomer.update_cropsheader()
            if args.individual is True:
                fout = (key + '_' + key2 + infixlbl["croprenum"] +
                        os.path.splitext(os.path.basename(inseq))[1])
                outseq = outpathgen(outdir, subdir=key, filename=fout,
                                    mksubdir=True)
                if monomer.cropmap is not None:
                    fout = (key + '_' + key2 + infixlbl["croprenum"] +
                            '.cropmap')
                    outmap = outpathgen(outdir, subdir=key, filename=fout,
                                        mksubdir=True)
            elif len(seqset) == 1 or args.sort is None:
                if len(seqset) > 1:
                    fout = (os.path.splitext(os.path.basename(inseq))[0] +
                            infixlbl["croprenum"] +
                            os.path.splitext(os.path.basename(inseq))[1])
                    outseq = outpathgen(outdir, filename=fout)
                    if monomer.cropmap is not None:
                        fout = (os.path.splitext(os.path.basename(inseq))[0] +
                                infixlbl["croprenum"] + '.cropmap')
                        outmap = outpathgen(outdir, filename=fout)
                else:
                    fout = (key + infixlbl["croprenum"] +
                            os.path.splitext(os.path.basename(inseq))[1])
                    outseq = outpathgen(outdir, subdir=key, filename=fout,
                                        mksubdir=True)
                    if monomer.cropmap is not None:
                        fout = key + infixlbl["croprenum"] + '.cropmap'
                        outmap = outpathgen(outdir, subdir=key,
                                            filename=fout, mksubdir=True)

            if len(seqset) > 1 and args.sort is not None:
                sorted_outseq[monomer.oligomer_id + '_' +
                              monomer.name] = monomer.deepcopy()
            else:
                monomer.dump(outseq)
                if monomer.cropmap is not None:
                    monomer.dumpmap(outmap)

    # PRINT OUT THE RESULTS AND FINISH
    croptime = time.time()
    logger.debug('Crop time = ' + str(croptime-starttime) + ' s')
    logger.info('Done'+os.linesep)

    if len(seqset) > 1 and args.sort is not None:
        logger.info('Sorting sequence(s)...')
        fout = (os.path.splitext(os.path.basename(inseq))[0] +
                infixlbl["croprenum"] + ".sorted_" +
                sorter + os.path.splitext(os.path.basename(inseq))[1])
        outseq = outpathgen(outdir, filename=fout)
        if sorter == 'ncrops':
            sorted_outseq2 = sorted(sorted_outseq.items(),
                                    key=lambda x: x[1].ncrops(),
                                    reverse=True)
        elif sorter == 'percent':
            sorted_outseq2 = sorted(sorted_outseq.items(),
                                    key=lambda x: x[1].ncrops()/x[1].full_length(),
                                    reverse=True)
        elif sorter == 'ncropsin':
            sorted_outseq2 = sorted(sorted_outseq.items(),
                                    key=lambda x: x[1].ncrops(offmidseq=True),
                                    reverse=True)
        elif sorter == 'percentin':
            sorted_outseq2 = sorted(sorted_outseq.items(),
                                    key=lambda x: x[1].ncrops(offmidseq=True)/x[1].full_length(),
                                    reverse=True)
        del sorted_outseq
        for monomer in sorted_outseq2:
            monomer.dump(outseq)
        logger.debug('Sort time = '+str(time.time()-croptime)+ ' s')
        logger.info('Done'+os.linesep)

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
