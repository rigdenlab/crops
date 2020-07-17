# -*- coding: utf-8 -*-

from crops.about import __prog__, __description__, __author__, __date__, __version__

from crops.core.rescodes import ressymbol
#from .sequence import Sequence
from crops.core.sequence import monomer_sequence
#from .intervals import intinterval

def renumberpdb(INSEQ,INSTR,seqback=False):
    """Returns modified :class:`gemmi.Structure` with new residue numbers.

    :param INSEQ: Input sequence.
    :type INSEQ: :class:`~crops.core.sequence.Sequence`
    :param INSTR: Gemmi structure.
    :type INSTR: :class:`gemmi.Structure`
    :param seqback: If True, it additionally returns the Sequence with the gaps found in the structure, defaults to False.
    :type seqback: bool, optional
    :return INSTR: Renumbered structure.
    :return INSEQ: Sequence with extra information about gaps, only if seqback==True.

    """

    n_chains = 0
    n_resmax = 0
    for model in INSTR:
        n_chains += len(model)
        for chain in model:
            if len(chain) > n_resmax:
                n_resmax = len(chain)
    pos = [[0 for j in range(n_resmax)] for i in range(n_chains)]
    n_chains = 0
    #NUMBER OF CHAINS PER MODEL ->> DO
    if seqback:
        for monomer in INSEQ.imer.values():
            monomer.seqs['gapseq']=[]

    for model in INSTR:
        for chain in model:
            original_seq=INSEQ.imer[chain.name].seqs['mainseq']
            solved = False
            for shift in range(int(len(chain)/2)):
                cnt=0
                gap=0
                score=0
                nligands=0
                newseq=''
                newseq += '-'*shift
                for residue in chain:
                    if residue == chain[0]:
                        if ressymbol(residue.name) == original_seq[shift]:
                            score += 1
                            pos[n_chains][cnt]=1+shift
                            newseq += ressymbol(residue.name)
                    elif ressymbol(residue.name)==0:
                        nligands+=1
                        pos[n_chains][cnt]=-nligands
                    else:
                        if (chain[cnt].seqid.num-chain[cnt-1].seqid.num > 1):
                            gap += (chain[cnt].seqid.num-chain[cnt-1].seqid.num-1)
                            newseq += '-'*(chain[cnt].seqid.num-chain[cnt-1].seqid.num-1)
                        pos[n_chains][cnt]=cnt+1+gap+shift
                        if ressymbol(residue.name) == original_seq[cnt+gap+shift]:
                            score += 1
                            newseq += ressymbol(residue.name)
                        if residue==chain[-1]:
                            if cnt+gap+shift+1 < len(original_seq):
                                newseq += '-'*(len(original_seq)-(cnt+gap+shift+1))
                    cnt += 1
                if score == len(chain)-nligands:
                    solved = True
                    break
            if solved:
                cnt=0
                for residue in chain:
                    residue.seqid.num = pos[n_chains][cnt]
                    cnt += 1
            if seqback:
                INSEQ.imer[chain.name].seqs['gapseq'].append(newseq)
            n_chains += 1
            solved = False

    if seqback:
        return INSTR, INSEQ
    else:
        return INSTR

def crop_seq(INSEQ, segments, cut_type, terms=False):  #INPUTS MUST BE SINGLE MONOMERS
    """Returns modified :class:`~crops.core.sequence.monomer_sequence` without specified elements.

    :param INSEQ: Input sequence.
    :type INSEQ: :class:`~crops.core.sequence.monomer_sequence`
    :param segments: Input preserving interval.
    :type segments: :class:`~crops.core.intervals.intinterval`
    :param cut_type: Additional header information.
    :type cut_type: str
    :param terms: If True, only terminal ends are removed, defaults to False.
    :type terms: bool, optional
    :raises ValueError: If intervals given lie out of the sequence.
    :return newchain: Cropped sequence.
    :rtype newchain: :class:`~crops.core.sequence.monomer_sequence`

    """
    if segments.subint[-1][-1] > INSEQ.length():
        raise ValueError('One or many of the segment end values is outside the original sequence.')

    newchain=monomer_sequence(chid=INSEQ.info['chain_id'],header=INSEQ.info['header'])
    newchain.seqs['fullseq']=INSEQ.seqs['mainseq']
    newchain.seqs['cropseq']=''
    if 'gapseq' in INSEQ.seqs:
        newchain.seqs['gapseq']=['' for i in range(len(INSEQ.seqs['gapseq']))]
        newchain.seqs['cropgapseq']=['' for i in range(len(INSEQ.seqs['gapseq']))]
    cropint=segments.deepcopy() if not terms else segments.union(segments.terminals())

    for res in range(INSEQ.length()):
        if cropint.contains(res+1):
            newchain.seqs['mainseq'] += INSEQ.seqs['mainseq'][res]
            if 'gapseq' in INSEQ.seqs:
                for n in range(len(INSEQ.seqs['gapseq'])):
                    newchain.seqs['gapseq'][n] += INSEQ.seqs['gapseq'][n][res]
                    newchain.seqs['cropgapseq'][n] += INSEQ.seqs['gapseq'][n][res]
            newchain.seqs['cropseq'] += INSEQ.seqs['main'][res]
        else:
            if 'gapseq' in INSEQ.seqs:
                for n in range(len(INSEQ.seqs['gapseq'])):
                    newchain.seqs['cropgapseq'][n] += '*'
            newchain.seqs['cropseq'] += '*'

    if newchain.length()<len(newchain.seqs['cropseq']):
        newchain.info['header'] += cut_type

    return newchain

def croppdb(INSTR, INSEQ, segments, terms=False):
    """Returns modified :class:`gemmi.Structure` without specified elements.

    :param INSTR: Gemmi structure.
    :type INSTR: :class:`gemmi.Structure`
    :param INSEQ: Input sequence.
    :type INSEQ: :class:`~crops.core.sequence.Sequence`
    :param segments: Input preserving interval.
    :type segments: :class:`~crops.core.intervals.intinterval`
    :param terms: If True, only terminal ends are removed, defaults to False.
    :type terms: bool, optional
    :return INSTR: DESCRIPTION
    :rtype INSTR: :class:`gemmi.Structure`

    """

    n_chains = 0
    n_resmax = 0

    for model in INSTR:
        n_chains += len(model)
        for chain in model:
            if len(chain) > n_resmax:
                n_resmax = len(chain)
        delres = [[False for j in range(n_resmax)] for i in range(n_chains)]
        n_chains = 0

    for model in INSTR:
        for chain in model:
            if chain.name in segments:
                if not terms:
                    cropint=segments[chain.name].deepcopy()
                else:
                    cropint=segments[chain.name].union(segments[chain.name].terminals())
                original_seq=INSEQ.imer[chain.name].seqs['mainseq']

                r_bio=0
                pos_chainlist=0
                for r_original in range(len(original_seq)):
                    if cropint.contains(r_original+1):
                        r_bio+=1
                        if chain[pos_chainlist].seqid.num == r_original+1:
                            chain[pos_chainlist].seqid.num=r_bio
                            pos_chainlist += 1
                    else:
                        if chain[pos_chainlist].seqid.num == r_original+1:
                            delres[n_chains][pos_chainlist] = True
                            pos_chainlist += 1

        n_chains = 0
        for model in INSTR:
            for chain in model:
                for res in reversed(range(len(chain))):
                    if delres[n_chains][res]:
                        del chain[res]
                n_chains += 1

    return INSTR

