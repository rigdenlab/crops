"""This script will remove a number of residues from sequence and structure files in agreement to the intervals and other details supplied."""

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

logger = None


def create_argument_parser():
    """Create a parser for the command line arguments used in crops-cropstr."""
    parser = argparse.ArgumentParser(prog=__prog__, formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=__description__+' ('+__prog__+')  v.'+__version__+os.linesep+__doc__)
    parser.add_argument("input_seqpath", nargs=1, metavar="Sequence_filepath",
                        help="Input sequence filepath.")
    parser.add_argument("input_strpath", nargs=1, metavar="Structure_filepath",
                        help="Input structure filepath or dir. If a directory is inserted, it will act on all structure files in such directory.")
    parser.add_argument("input_database", nargs=1, metavar="Intervals_database",
                        help="Input intervals database filepath.")

    parser.add_argument("-o", "--outdir", nargs=1, metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence.")

    parser.add_argument("-p", "--preselect", nargs='+', metavar="Oligoseq_ids", default=None,
                        help="From all the sequences in the input sequence file, just print out this preselected subset.")
    parser.add_argument("-i", "--individual", action='store_true', default=False,
                        help="One separated output fasta file per each sequence.")

    sections = parser.add_mutually_exclusive_group(required=False)
    sections.add_argument("-t", "--terminals", action='store_true', default=False,
                          help="Ignore interval discontinuities and only crop the ends off.")
    sections.add_argument("-u", "--uniprot_threshold", nargs=2, metavar=("Uniprot_ratio_threshold", "Sequence_database"),
                          help='Act if SIFTS database is used as intervals source AND %% residues from single Uniprot sequence is above threshold. Threshold: [MIN,MAX)=[0,100). Database path: uniclust##_yyyy_mm_consensus.fasta-path or server-only. The latter requires internet connexion.')

    parser.add_argument('--version', action='version', version='%(prog)s ' + __version__)

    return parser


def main():
    """Remove a number of residues from sequence and structure files in agreement to the intervals and other details supplied.

    :raises ValueError: For wrong argument values.

    """
    # INITIALISE AND PARSE ARGUMENTS FROM COMMAND LINE
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = ccl.crops_logger(level="info")
    logger.info(ccl._welcome())

    inseq = check_path(args.input_seqpath[0], 'file')
    indb = check_path(args.input_database[0], 'file')
    instr = check_path(args.input_strpath[0])
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

    # PARSE INPUT FILES
    logger.info('Parsing sequence file ' + inseq)
    if args.preselect is not None:
        subset = set(args.preselect)
    else:
        subset = None
    seqset = cin.parseseqfile(seq_input=inseq, inset=subset)
    logger.info('Done')

    logger.info('Parsing structure file ' + instr)
    strset, fileset = cin.parsestrfile(instr)
    logger.info('Done')

    logger.info('Parsing interval database file ' + indb)
    if len(seqset) > 0:
        intervals = cin.import_db(indb, pdb_in=seqset)
    else:
        logger.critical('No chains were imported from sequence file.')
        raise ValueError
    logger.info('Done'+os.linesep)

    if insprot is not None and minlen > 0.0:
        logger.info('Parsing uniprot sequence file: ' + insprot)
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
        uniprotset = cin.parseseqfile(seq_input=insprot, inset=uniprotset, use_UPserver=upserver)
        logger.info('Done'+os.linesep)

    # MAIN OPERATION / PRINT OUT RESULTS WITHIN
    gseqset = {}
    strset2 = {}
    logger.info('Renumbering structure(s)...')
    for key, structure in strset.items():
        found = False
        for seqname in seqset:
            if ((seqname in key) or
                    (len(seqset) == 1 and len(strset) == 1)):
                finalid = seqname
                newstructure, gseqset[seqname] = cop.renumber_pdb(seqset[seqname],
                                                                  structure,
                                                                  seqback=True)
                fout = finalid + infixlbl["renumber"] + os.path.splitext(instr)[1]
                outstr = outpathgen(outdir, subdir=finalid,
                                    filename=fout, mksubdir=True)
                newstructure.write_minimal_pdb(outstr)
                strset2[finalid] = structure
                found = True
        if found is False:
            logger.warning("Identifier '"+key+"' not found in sequence input.")

    logger.info('Done'+os.linesep)
    logger.info('Cropping renumbered structure(s)...')
    outseq = os.path.join(outdir, os.path.splitext(os.path.basename(inseq))[0] +
                          infixlbl["croprenum"] +
                          os.path.splitext(os.path.basename(inseq))[1])
    for key, S in gseqset.items():
        newS = S.deepcopy()
        if key in intervals:
            if insprot is not None and minlen > 0.0:
                newinterval = {}
            for key2, monomer in S.imer.items():
                cropped_seq = False
                for key3 in monomer.chains:
                    if key3 in intervals[key]:
                        if insprot is not None and minlen > 0.0:
                            newinterval[key3] = intervals[key][key3].deepcopy()
                            newinterval[key3].tags['description'] += ' - Uniprot threshold'
                            newinterval[key3].subint = []
                            unilbl = ' uniprot chains included: '
                            for unicode, uniintervals in intervals[key][key3].tags['uniprot'].items():
                                uniseq = uniprotset[unicode].imer['1']
                                if 100*uniintervals.n_elements()/uniseq.length() >= minlen:
                                    newinterval[key3] = newinterval[key3].union(intervals[key][key3].intersection(uniintervals))
                                    unilbl += unicode + '|'
                            if cropped_seq is False:
                                monomer = cop.crop_seq(monomer,
                                                       newinterval[key3],
                                                       targetlbl+unilbl,
                                                       terms=args.terminals)
                                cropped_seq = True
                        else:
                            if cropped_seq is False:
                                monomer = cop.crop_seq(monomer,
                                                       intervals[key][key3],
                                                       targetlbl,
                                                       terms=args.terminals)
                                cropped_seq = True
                        if newS.imer[key2] != monomer:
                            newS.imer[key2] = monomer.deepcopy()
                    else:
                        logger.warning('Chain-name ' + key + '_' + str(key3) +
                                       ' not found in database. Cropping not performed.')

                monomer.update_cropsheader()

                hf = '_' + key2 if args.individual is True else ''
                ifx = infixlbl["croprenum"] if cropped_seq is True else ''
                fout = (key + hf + ifx +
                        os.path.splitext(os.path.basename(inseq))[1])
                outseq = outpathgen(outdir, subdir=key, filename=fout,
                                    mksubdir=True)
                monomer.dump(outseq)
                if monomer.cropmap is not None:
                    fout = key + hf + infixlbl["croprenum"] + '.cropmap'
                    outmap = outpathgen(outdir, subdir=key,
                                        filename=fout)
                    monomer.dumpmap(outmap)

            cropped_str = cop.crop_pdb(strset2[key], newS, original_id=True)
            fout = key + infixlbl["crop"] + os.path.splitext(instr)[1]
            outstr = outpathgen(outdir, subdir=key,
                                filename=fout, mksubdir=True)
            cropped_str.write_minimal_pdb(outstr)

            cropped_str2 = cop.crop_pdb(strset2[key], newS, original_id=False)

            fout = key + infixlbl["croprenum"] + os.path.splitext(instr)[1]
            outstr = outpathgen(outdir, subdir=key,
                                filename=fout, mksubdir=True)
            cropped_str2.write_minimal_pdb(outstr)
        else:
            logger.warning('PDB-ID ' + key.upper() +
                           ' not found in database. Cropping not performed.')
            for key2, monomer in newS.imer.items():
                hf = '_' + key2 if args.individual is True else ''
                fout = key + hf + os.path.splitext(os.path.basename(inseq))[1]
                outseq = outpathgen(outdir, subdir=key,
                                    filename=fout, mksubdir=True)
                monomer.dump(outseq)

    # FINISH
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
