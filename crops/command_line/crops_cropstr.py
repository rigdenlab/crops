# -*- coding: utf-8 -*-

"""==========
This script will remove a number of residues from a sequence file
in agreement to the intervals and other details supplied.

"""

from crops.about import __prog__, __description__, __author__, __date__, __version__

import argparse
import os
#import copy
#import gemmi
from warnings import warn

from crops.core import cio
from crops.core import ops as cop
#from core import seq as csq

def main():
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


    args = parser.parse_args()

    inseq=cio.check_path(args.input_seqpath[0],'file')
    indb=cio.check_path(args.input_database[0],'file')
    instr=cio.check_path(args.input_strpath[0])
    insprot=cio.check_path(args.uniprot_threshold[1]) if args.uniprot_threshold is not None else None

    minlen=float(args.uniprot_threshold[0]) if args.uniprot_threshold is not None else 0.0
    targetlbl=cio.target_format(indb,terms=args.terminals, th=minlen)
    infixlbl=cio.infix_gen(indb,terms=args.terminals)

    if args.outdir is None:
        outdir=cio.check_path(os.path.dirname(inseq),'dir')
    else:
        outdir=cio.check_path(os.path.dirname(args.outdir[0]),'dir')
    ###########################################
    seqset=cio.parseseqfile(inseq)
    strset, fileset=cio.parsestrfile(instr)

    if len(seqset)==1:
        for key in seqset:
            pdbid=key
        intervals=cio.import_db(indb,pdbid)
    elif len(seqset)>1:
        intervals=cio.import_db(indb)
    else:
        raise ValueError('No chains were imported from sequence file.')

    if insprot is not None and minlen>0.0:
        uplist=[]
        for seqncid, seqnc in seqset.items():
            for monomerid, monomer in seqnc.imer.items():
                if 'uniprot' in intervals[seqncid][monomerid].tags:
                    for key in intervals[seqncid][monomerid].tags['uniprot']:
                        if key not in uplist:
                            uplist.append(key)

        uniprotset=cio.parseseqfile(insprot, uniprot=uplist)['uniprot']
    ###########################################
    gseqset={}
    for pdbid, structure in strset.items():
        if pdbid in seqset:
            newstructure,gseqset[pdbid]=cop.renumberpdb(seqset[pdbid],structure,seqback=True)
            outstr=cio.outpath(outdir,subdir=pdbid,filename=pdbid+infixlbl+os.path.splitext(instr)[1],mksubdir=True)
            newstructure.write_pdb(outstr)

    outseq=os.path.join(outdir,os.path.splitext(os.path.basename(inseq))[0]+infixlbl["crop"]+os.path.splitext(os.path.basename(inseq))[1])
    for pdbid, S in gseqset.items():
        if pdbid in intervals:
            for chainid,monomer in S.imer.items():
                if chainid in intervals[pdbid]:
                    if insprot is not None:
                        newinterval=intervals[pdbid][chainid].deepcopy(newdescription=intervals[pdbid][chainid].tags['description'] + ' - Uniprot threshold')
                        newinterval.subint=[]
                        unilbl=' uniprot chains included: '
                        for unicode,uniends in intervals[pdbid][chainid].tags['uniprot'].items():
                            if 100*(uniends[1]-uniends[0]+1)/uniprotset[unicode].length(unicode)>=minlen:
                                newinterval=newinterval.union(intervals[pdbid][chainid].intersection(uniends))
                                unilbl+=unicode +'|'
                        monomer=cop.crop_seq(monomer,newinterval,targetlbl+unilbl,terms=args.terminals)
                    else:
                        monomer=cop.crop_seq(monomer,intervals[pdbid][chainid],targetlbl,terms=args.terminals)
                else:
                    warn('Chain name '+pdbid+'_'+str(chainid)+' not found in database. Cropping not performed.')
                monomer=cop.crop_seq(monomer,intervals[pdbid][chainid],targetlbl,terms=args.terminals)
                outseq=cio.outpath(outdir,subdir=pdbid,filename=pdbid+infixlbl["crop"]+os.path.splitext(os.path.basename(inseq))[1])
                monomer.dump(outseq)
            cropped_str=cop.croppdb(strset[pdbid],S,intervals[pdbid],args.terminals)
            outstr=cio.outpath(outdir,subdir=pdbid,filename=pdbid+infixlbl["crop"]+os.path.splitext(instr)[1],mksubdir=True)
            cropped_str.write_pdb(outstr)
        else:
            warn('PDB ID '+pdbid+' not found in database. Cropping not performed.')
            for chainid,monomer in S.imer.items():
                outseq=cio.outpath(outdir,subdir=pdbid,filename=pdbid+os.path.splitext(os.path.basename(inseq))[1])
                monomer.dump(outseq)

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