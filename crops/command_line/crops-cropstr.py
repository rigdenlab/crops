"""==========
This script will remove a number of residues from a sequence file
in agreement to the intervals and other details supplied.

"""

from crops.about import __prog__, __description__, __author__, __date__, __version__

import argparse
import os

from crops.io import check_path
from crops.io import outpathgen
from crops.io import parsers as cin
from crops.io import taggers as ctg
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
    parser.add_argument("input_database",nargs=1, metavar="Intervals_database",
                        help="Input intervals database filepath.")

    parser.add_argument("-o","--outdir",nargs=1,metavar="Output_Directory",
                        help="Set output directory path. If not supplied, default is the one containing the input sequence.")

    sections=parser.add_mutually_exclusive_group(required=False)
    sections.add_argument("-t","--terminals",action='store_true',default=False,
                          help="Ignore interval discontinuities and only crop the ends off.")
    sections.add_argument("-u","--uniprot_threshold", nargs=2, metavar=("Uniprot_ratio_threshold","Sequence_database"),
                          help='Act if SIFTS database is used as intervals source AND %% residues from single Uniprot sequence is above threshold. [MIN,MAX)=[0,100) uniprot_sprot.fasta-path')

    parser.add_argument('--version', action='version', version='%(prog)s '+ __version__)

    return parser

def main():
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = ccl.crops_logger(level="info")
    logger.info(ccl.welcome())

    inseq=check_path(args.input_seqpath[0],'file')
    indb=check_path(args.input_database[0],'file')
    instr=check_path(args.input_strpath[0])
    insprot=check_path(args.uniprot_threshold[1]) if args.uniprot_threshold is not None else None

    minlen=float(args.uniprot_threshold[0]) if args.uniprot_threshold is not None else 0.0
    targetlbl=ctg.target_format(indb,terms=args.terminals, th=minlen)
    infixlbl=ctg.infix_gen(indb,terms=args.terminals)

    if args.outdir is None:
        outdir=check_path(os.path.dirname(inseq),'dir')
    else:
        outdir=check_path(os.path.join(args.outdir[0],''),'dir')
    ###########################################
    logger.info('Parsing sequence file '+inseq)
    seqset=cin.parseseqfile(inseq)
    logger.info('Done')

    logger.info('Parsing structure file '+instr)
    strset, fileset=cin.parsestrfile(instr)
    logger.info('Done')

    logger.info('Parsing interval database file '+indb)
    if len(seqset)>0:
        intervals=cin.import_db(indb,pdb_in=seqset)
    else:
        raise ValueError('No chains were imported from sequence file.')
    logger.info('Done\n')

    if insprot is not None and minlen>0.0:
        logger.info('Parsing uniprot sequence file '+insprot)
        uniprotset={}
        for seqncid, seqnc in seqset.items():
            for monomerid, monomer in seqnc.imer.items():
                if 'uniprot' in intervals[seqncid][monomerid].tags:
                    for key in intervals[seqncid][monomerid].tags['uniprot']:
                        if key.upper() not in uniprotset:
                            uniprotset[key.upper()]=None

        uniprotset=cin.parseseqfile(insprot, uniprot=uniprotset)['uniprot']
        logger.info('Done\n')

    ###########################################
    gseqset={}
    logger.info('Renumbering structure(s)...')
    for key, structure in strset.items():
        if key in seqset:
            newstructure,gseqset[key]=cop.renumber_pdb(seqset[key],structure,seqback=True)
            outstr=outpathgen(outdir,subdir=key,filename=key+infixlbl["renumber"]+os.path.splitext(instr)[1],mksubdir=True)
            #newstructure.write_pdb(outstr)
            newstructure.write_minimal_pdb(outstr)
    logger.info('Done\n')
    logger.info('Cropping renumbered structure(s)...')
    outseq=os.path.join(outdir,os.path.splitext(os.path.basename(inseq))[0]+infixlbl["crop"]+os.path.splitext(os.path.basename(inseq))[1])
    for key, S in gseqset.items():
        newS=S.deepcopy()
        if key in intervals:
            if insprot is not None and minlen>0.0:
                newinterval={}
            for key2,monomer in S.imer.items():
                if key2 in intervals[key]:
                    if insprot is not None and minlen>0.0:
                        newinterval[key2]=intervals[key][key2].deepcopy()
                        newinterval[key2].tags['description']+=' - Uniprot threshold'
                        newinterval[key2].subint=[]
                        unilbl=' uniprot chains included: '
                        for unicode,uniintervals in intervals[key][key2].tags['uniprot'].items():
                            if 100*uniintervals.n_elements()/uniprotset.imer[unicode].length()>=minlen:
                                newinterval[key2]=newinterval[key2].union(intervals[key][key2].intersection(uniintervals))
                                unilbl+=unicode +'|'
                        monomer=cop.crop_seq(monomer,newinterval[key2],targetlbl+unilbl,terms=args.terminals)
                    else:
                        monomer=cop.crop_seq(monomer,intervals[key][key2],targetlbl,terms=args.terminals)
                    newS.imer[key2]=monomer.deepcopy()
                else:
                    logger.warning('Chain-name '+key+'_'+str(key2)+' not found in database. Cropping not performed.')
                outseq=outpathgen(outdir,subdir=key,filename=key+infixlbl["croprenum"]+os.path.splitext(os.path.basename(inseq))[1])
                monomer.dump(outseq)
            if insprot is not None and minlen>0.0:
                cropped_str=cop.crop_pdb(strset[key],newS,original_id=True)
            else:
                cropped_str=cop.crop_pdb(strset[key],newS,original_id=True)
            outstr=outpathgen(outdir,subdir=key,filename=key+infixlbl["crop"]+os.path.splitext(instr)[1],mksubdir=True)
            #cropped_str.write_pdb(outstr)
            cropped_str.write_minimal_pdb(outstr)
            if insprot is not None and minlen>0.0:
                cropped_str2=cop.crop_pdb(strset[key],newS,original_id=False)
            else:
                cropped_str2=cop.crop_pdb(strset[key],newS,original_id=False)
            outstr=outpathgen(outdir,subdir=key,filename=key+infixlbl["croprenum"]+os.path.splitext(instr)[1],mksubdir=True)
            #cropped_str2.write_pdb(outstr)
            cropped_str2.write_minimal_pdb(outstr)
        else:
            logger.warning('PDB-ID '+key.upper()+' not found in database. Cropping not performed.')
            for key2,monomer in newS.imer.items():
                outseq=outpathgen(outdir,subdir=key,filename=key+os.path.splitext(os.path.basename(inseq))[1])
                monomer.dump(outseq)
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
