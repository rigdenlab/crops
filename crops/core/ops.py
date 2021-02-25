from crops.about import __prog__, __description__, __author__, __date__, __version__

from crops.libs import ressymbol
from crops.elements.sequence import monomer_sequence

import copy
import logging

def renumber_pdb(inseq,instr,seqback=False):
    """Returns modified :class:`gemmi.Structure` with new residue numbers.

    :param inseq: Input sequence.
    :type inseq: :class:`~crops.elements.sequence.Sequence`
    :param instr: Gemmi structure.
    :type instr: :class:`gemmi.Structure`
    :param seqback: If True, it additionally returns the Sequence with the gaps found in the structure, defaults to False.
    :type seqback: bool, optional
    :return instr: Renumbered structure.
    :rtype instr: :class:`gemmi.Structure`
    :return inseq: Sequence with extra information about gaps, only if seqback==True.
    :rtype inseq: :class:`~crops.elements.sequence.Sequence`

    """

    n_chains = 0
    n_resmax = 0
    for model in instr:
        n_chains += len(model)
        for chain in model:
            if len(chain) > n_resmax:
                n_resmax = len(chain)
    pos = [[0 for j in range(n_resmax)] for i in range(n_chains)]
    n_chains = 0
    if seqback:
        for monkey in inseq.imer:
            inseq.imer[monkey].seqs['gapseq']=[]

    for model in instr:
        for chain in model:
            solved = False
            if chain.name in inseq.imer:
                original_seq=inseq.imer[chain.name].seqs['mainseq']
                for shift in range(int(len(chain)/2)):
                    cnt=0
                    gap=0
                    score=0
                    nligands=0
                    newseq = '-'*shift
                    for residue in chain:
                        if str(residue.entity_type)!='EntityType.Polymer':
                            nligands+=1
                            pos[n_chains][cnt]=-nligands
                        else:
                            if residue != chain[0] and chain[cnt].seqid.num-chain[cnt-1].seqid.num > 1:
                                gap += (chain[cnt].seqid.num-chain[cnt-1].seqid.num-1)
                                newseq += '-'*(chain[cnt].seqid.num-chain[cnt-1].seqid.num-1)
                            pos[n_chains][cnt]=cnt+1+gap+shift
                            if (ressymbol(residue.name,pick=original_seq[cnt+gap+shift]) == original_seq[cnt+gap+shift]
                                or ressymbol(residue.name,pick=original_seq[cnt+gap+shift]) == 0): #second condition ignores possible gaps in crops.rescodes database
                                score += 1
                                newseq += original_seq[cnt+gap+shift]
                            else:
                                pdbseq=instr.get_entity(chain.name).full_sequence
                                if cnt>len(pdbseq):
                                    nligands+=1
                                    pos[n_chains][cnt]=-nligands
                                else:
                                    res_pdbid=pdbseq[cnt].split(",")
                                    if len(res_pdbid)>1:
                                        for n in range(1,len(res_pdbid)):
                                            if ressymbol(res_pdbid[n],pick=original_seq[cnt+gap+shift]) == original_seq[cnt+gap+shift]:
                                                score += 1
                                                newseq += original_seq[cnt+gap+shift]
                                                break
                        cnt += 1
                    if score == len(chain)-nligands:
                        solved = True
                        if len(newseq)<len(original_seq):
                            newseq += '-'*(len(original_seq)-len(newseq))
                        break
                if solved == False:
                    raise ValueError('The .fasta sequence and the structure given do not match.')
            else:
                logging.warning('.pdb chain '+str(chain.name)+' not found in .fasta file. All elements considered ligands.')
                ligandwarn=False
                nligands=0
                for residue in chain:
                    if ressymbol(residue.name)!=0:
                        ligandwarn=True
                    nligands+=1
                    pos[n_chains][nligands-1]=-nligands
                if ligandwarn==True:
                    logging.warning('Some of the ligands contain Aminoacid or Nucleotide residues. Please check that they actually are ligands.')
                solved=True

            if solved:
                cnt=0
                for residue in chain:
                    if pos[n_chains][cnt]<0:
                        residue.seqid.num = len(newseq)-pos[n_chains][cnt]
                    else:
                        residue.seqid.num = pos[n_chains][cnt]
                    cnt += 1
            if seqback:
                if chain.name in inseq.imer:
                    inseq.imer[chain.name].seqs['gapseq'].append(newseq)
            n_chains += 1
            solved = False
    if seqback:
        return instr,inseq
    else:
        return instr

def crop_seq(inseq, segments, cut_type, terms=False):  #INPUTS MUST BE SINGLE MONOMERS
    """Returns modified :class:`~crops.elements.sequence.monomer_sequence` without specified elements.

    :param inseq: Input sequence.
    :type inseq: :class:`~crops.elements.sequence.monomer_sequence`
    :param segments: Input preservation interval.
    :type segments: :class:`~crops.elements.intervals.intinterval`
    :param cut_type: Additional header information.
    :type cut_type: str
    :param terms: If True, only terminal ends are removed, defaults to False.
    :type terms: bool, optional
    :raises ValueError: If intervals given lie out of the sequence.
    :return: Cropped sequence.
    :rtype: :class:`~crops.elements.sequence.monomer_sequence`

    """
    if len(segments.subint)>0:
        if segments.subint[-1][-1] > inseq.length():
            logging.debug('On '+str(inseq))
            logging.debug('with '+str(segments))
            raise ValueError('One or many of the segment end values is outside the original sequence.')

    newchain=monomer_sequence(chid=inseq.info['chain_id'],header=inseq.info['header'])
    newchain.seqs['fullseq']=inseq.seqs['mainseq']
    newchain.seqs['cropseq']=''

    if 'gapseq' in inseq.seqs:
        newchain.seqs['gapseq']=['']*len(inseq.seqs['gapseq'])
        newchain.seqs['cropgapseq']=['']*len(inseq.seqs['gapseq'])

    cropint=segments.deepcopy() if not terms else segments.union(segments.terminals())
    cropmap={}
    for res in range(inseq.length()):
        cropmap[res+1]=None

    for res in range(inseq.length()):
        if cropint.contains(res+1):
            newchain.seqs['mainseq'] += inseq.seqs['mainseq'][res]
            cropmap[res+1]=len(newchain.seqs['mainseq'])
            if 'gapseq' in inseq.seqs:
                for n in range(len(inseq.seqs['gapseq'])):
                    newchain.seqs['gapseq'][n] += inseq.seqs['gapseq'][n][res]
                    newchain.seqs['cropgapseq'][n] += inseq.seqs['gapseq'][n][res]
            newchain.seqs['cropseq'] += inseq.seqs['mainseq'][res]
        else:
            if 'gapseq' in inseq.seqs:
                for n in range(len(inseq.seqs['gapseq'])):
                    newchain.seqs['cropgapseq'][n] += '*'
            newchain.seqs['cropseq'] += '*'

    if newchain.length()<len(newchain.seqs['cropseq']):
        newchain.info['header'] += cut_type
    newchain.info['cropmap']=copy.deepcopy(cropmap)

    return newchain

def crop_pdb(instr, inseq, original_id=True):
    """Returns modified :class:`gemmi.Structure` without specified elements.

    :param instr: Gemmi structure.
    :type instr: :class:`gemmi.Structure`
    :param inseq: Input previously-cropped-sequence.
    :type inseq: :class:`~crops.elements.sequence.Sequence`
    :param original_id: If True, it will keep residue ids alligned to original sequence, defaults to True
    :type original_id: bool, optional
    :return: Cropped structure.
    :rtype: :class:`gemmi.Structure`

    """
    for model in instr:
        for chain in model:
            if chain.name in inseq.imer:
                if 'cropmap' in inseq.imer[chain.name].info:
                    for res in reversed(range(len(chain))):
                        if chain[res].seqid.num in inseq.imer[chain.name].info['cropmap']:
                            if inseq.imer[chain.name].info['cropmap'][chain[res].seqid.num] is not None:
                                if not original_id:
                                    chain[res].seqid.num=inseq.imer[chain.name].info['cropmap'][chain[res].seqid.num]
                            else:
                                del chain[res]

    return instr

