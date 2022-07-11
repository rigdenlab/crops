from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

from crops.libs import ressymbol
from crops.elements.sequences import sequence

from Bio import Align
from Bio.Align import substitution_matrices

import copy
import gemmi
import logging


def get_sequence_alignment(sequence_1, sequence_2, mode='global', open_gap_score=-11, extend_gap_score=-2):
    """Perform a sequence alignment using Needleman-Wunsch algorithm.

    :param sequence_1: First input sequence.
    :type sequence_1: str
    :param sequence_2: Second input sequence.
    :type sequence_2: str
    :param mode: Alignment mode, defaults to global.
    :type mode: str, optional
    :param open_gap_score: Opening gap penalty, defaults to -11.
    :type open_gap_score: int
    :param extend_gap_score: Extension gap penalty, defaults to -2.
    :type extend_gap_score: int
    :return alignment_dict: Dictionary with the residue mapping between both input sequences.
    :rtype alignment_dict: dict

    """

    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = open_gap_score
    aligner.extend_gap_score = extend_gap_score
    try:
        alignments = list(aligner.align(sequence_1, sequence_2))
    except ValueError as e:
        logging.warning('Needleman-Wunsch alignment failed due to wrong alphabet:\n{}'.format(e))
        return None
    alignments.sort(key=lambda x: x.score, reverse=True)
    aligned_indices = alignments[0].aligned
    alignment_dict = {}

    for query_chunk, target_chunk in zip(*aligned_indices):
        for query_index, target_index in zip(range(*query_chunk), range(*target_chunk)):
            alignment_dict[target_index] = query_index

    return alignment_dict

def renumber_pdb_needleman(inseq, instr, seqback=False):
    """Returns modified :class:`gemmi.Structure` with new residue numbers. It uses Needleman-Wunsch
    algorithm to perform the sequence alignment

    :param inseq: Input sequence.
    :type inseq: :class:`~crops.elements.sequences.oligoseq`
    :param instr: Gemmi structure.
    :type instr: :class:`gemmi.Structure`
    :param seqback: If True, it additionally returns the :class:`~crops.elements.sequences.oligoseq` with the gaps found in the structure, defaults to False.
    :type seqback: bool, optional
    :return instr: Renumbered structure.
    :rtype instr: :class:`gemmi.Structure`
    :return inseq: Sequence with extra information about gaps, only if seqback==True.
    :rtype inseq: :class:`~crops.elements.sequences.oligoseq`

    """

    renumbered_structure = gemmi.Structure()
    renumbered_structure.name = instr.name

    if seqback:
        for monkey in inseq.imer:
            inseq.imer[monkey].seqs['gapseq'] = []

    for model in instr:
        renumbered_model = gemmi.Model(model.name)
        for chain in model:
            renumbered_chain = gemmi.Chain(chain.name)
            nseq = inseq.whatseq(chain.name)
            original_seq = inseq.imer[nseq].seqs['mainseq']
            model_seq = ''.join([ressymbol(x.name) for x in chain])
            aligned_dict = get_sequence_alignment(original_seq, model_seq)
            if aligned_dict is None:
                logging.warning('Alignment failed for chain {}, it will be excluded'.format(chain.name))
                continue
            for index, residue in enumerate(chain):
                _residue = residue.clone()
                _residue.seqid.num = aligned_dict[index] + 1
                renumbered_chain.add_residue(_residue)
            renumbered_model.add_chain(renumbered_chain)
            if seqback:
                res_set = {residue.seqid.num for residue in renumbered_chain}
                newseq = ''
                for n in range(len(original_seq)):
                    if (n+1) in res_set:
                        newseq += original_seq[n]
                    else:
                        newseq += '-'
                inseq.imer[nseq].seqs['gapseq'].append(newseq)
        renumbered_structure.add_model(renumbered_model)

    if seqback:
        return renumbered_structure, inseq
    else:
        return renumbered_structure

def renumber_pdb(inseq, instr, seqback=False):
    """Returns modified :class:`gemmi.Structure` with new residue numbers.

    :param inseq: Input sequence.
    :type inseq: :class:`~crops.elements.sequences.oligoseq`
    :param instr: Gemmi structure.
    :type instr: :class:`gemmi.Structure`
    :param seqback: If True, it additionally returns the :class:`~crops.elements.sequences.oligoseq` with the gaps found in the structure, defaults to False.
    :type seqback: bool, optional
    :return instr: Renumbered structure.
    :rtype instr: :class:`gemmi.Structure`
    :return inseq: Sequence with extra information about gaps, only if seqback==True.
    :rtype inseq: :class:`~crops.elements.sequences.oligoseq`

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
            inseq.imer[monkey].seqs['gapseq'] = []

    for model in instr:
        for chain in model:
            solved = False
            nseq = inseq.whatseq(chain.name)
            if nseq in inseq.imer:
                original_seq = inseq.imer[nseq].seqs['mainseq']
                for shift in range(int(len(chain)/2)):
                    cnt = 0
                    gap = 0
                    score = 0
                    nligands = 0
                    newseq = '-'*shift
                    for residue in chain:
                        if str(residue.entity_type) != 'EntityType.Polymer':
                            nligands += 1
                            pos[n_chains][cnt] = -nligands
                        else:
                            if residue != chain[0] and chain[cnt].seqid.num-chain[cnt-1].seqid.num > 1:
                                gap += (chain[cnt].seqid.num-chain[cnt-1].seqid.num-1)
                                newseq += '-'*(chain[cnt].seqid.num-chain[cnt-1].seqid.num-1)
                            pos[n_chains][cnt] = cnt+1+gap+shift
                            if (ressymbol(residue.name,pick=original_seq[cnt+gap+shift]) == original_seq[cnt+gap+shift]
                                or ressymbol(residue.name,pick=original_seq[cnt+gap+shift]) == 0): #second condition ignores possible gaps in crops.rescodes database
                                score += 1
                                newseq += original_seq[cnt+gap+shift]
                            else:
                                pdbseq = instr.get_entity(chain.name).full_sequence
                                if cnt > len(pdbseq):
                                    nligands += 1
                                    pos[n_chains][cnt] = -nligands
                                else:
                                    res_pdbid = pdbseq[cnt].split(",")
                                    if len(res_pdbid) > 1:
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
                    nerr = ("The .fasta sequence and the original "
                            "structure's numbering do not match. The number "
                            "of gaps in the structure must be consistent "
                            "with the number of residues in the sequence "
                            "in order for CROPS to make the renumbering.")
                    logging.critical(nerr)
                    raise ValueError
            else:
                nwarn = ('.pdb chain '+str(chain.name) + ' not found in '
                         '.fasta file. All elements considered ligands.')
                logging.warning(nwarn)
                ligandwarn = False
                nligands = 0
                for residue in chain:
                    if ressymbol(residue.name) != 0:
                        ligandwarn = True
                    nligands += 1
                    pos[n_chains][nligands-1] = -nligands
                if ligandwarn == True:
                    nwarn = ('Some of the ligands contain Aminoacid or '
                             'Nucleotide residues. Please check whether '
                             'they are actually ligands.')
                    logging.warning(nwarn)
                solved = True

            if solved:
                cnt = 0
                for residue in chain:
                    if pos[n_chains][cnt] < 0:
                        residue.seqid.num = len(newseq)-pos[n_chains][cnt]
                    else:
                        residue.seqid.num = pos[n_chains][cnt]
                    cnt += 1
            if seqback:
                if nseq in inseq.imer:
                    inseq.imer[nseq].seqs['gapseq'].append(newseq)
            n_chains += 1
            solved = False
    if seqback:
        return instr, inseq
    else:
        return instr


def crop_seq(inseq, segments, cut_type, terms=False):
    """Returns modified :class:`~crops.elements.sequences.sequence` without specified elements.

    :param inseq: Input sequence.
    :type inseq: :class:`~crops.elements.sequences.sequence`
    :param segments: Input preservation interval.
    :type segments: :class:`~crops.elements.intervals.intinterval`
    :param cut_type: Additional header information.
    :type cut_type: str
    :param terms: If True, only terminal ends are removed, defaults to False.
    :type terms: bool, optional
    :raises ValueError: If intervals given lie out of the sequence.
    :return: Cropped sequence.
    :rtype: :class:`~crops.elements.sequences.sequence`

    """
    if len(segments.subint) > 0:
        if segments.subint[-1][-1] > inseq.length():
            logging.debug('On '+str(inseq))
            logging.debug('with '+str(segments))
            logging.critical('One or many of the segment end values is '
                             'outside the original sequence.')
            raise ValueError

    newchain = inseq.deepcopy()
    newchain.seqs['mainseq'] = ''
    newchain.seqs['fullseq'] = inseq.seqs['mainseq']
    newchain.seqs['cropseq'] = ''

    if 'gapseq' in inseq.seqs:
        newchain.seqs['gapseq'] = ['']*len(inseq.seqs['gapseq'])
        newchain.seqs['cropgapseq'] = ['']*len(inseq.seqs['gapseq'])

    cropint = segments.deepcopy() if not terms else segments.union(segments.terminals())
    cropmap = {}
    cropmap['map'] = {}
    cropmap['backmap'] = {}
    for res in range(inseq.length()):
        cropmap['map'][res+1] = None

    for res in range(inseq.length()):
        if cropint.contains(res+1):
            newchain.seqs['mainseq'] += inseq.seqs['mainseq'][res]
            cropmap['map'][res+1] = len(newchain.seqs['mainseq'])
            cropmap['backmap'][len(newchain.seqs['mainseq'])] = res+1
            if 'gapseq' in inseq.seqs:
                for n in range(len(inseq.seqs['gapseq'])):
                    newchain.seqs['gapseq'][n] += inseq.seqs['gapseq'][n][res]
                    newchain.seqs['cropgapseq'][n] += inseq.seqs['gapseq'][n][res]
            newchain.seqs['cropseq'] += inseq.seqs['mainseq'][res]
        else:
            if 'gapseq' in inseq.seqs:
                for n in range(len(inseq.seqs['gapseq'])):
                    newchain.seqs['cropgapseq'][n] += '+'
            newchain.seqs['cropseq'] += '+'

    newchain.infostring += cut_type

    newchain.cropmap = copy.deepcopy(cropmap['map'])
    newchain.cropbackmap = copy.deepcopy(cropmap['backmap'])

    return newchain

def crop_pdb(instr, inseq, original_id=True):
    """Returns modified :class:`gemmi.Structure` without specified elements.

    :param instr: Gemmi structure.
    :type instr: :class:`gemmi.Structure`
    :param inseq: Input previously-cropped-sequence.
    :type inseq: :class:`~crops.elements.sequences.oligoseq`
    :param original_id: If True, it will keep residue ids alligned to original sequence, defaults to True
    :type original_id: bool, optional
    :return: Cropped structure.
    :rtype: :class:`gemmi.Structure`

    """
    for model in instr:
        for chain in model:
            if chain.name in inseq.chainlist():
                chseq = inseq.whatseq(chain.name)
                if inseq.imer[chseq].cropmap is not None:
                    for res in reversed(range(len(chain))):
                        if chain[res].seqid.num in inseq.imer[chseq].cropmap:
                            if inseq.imer[chseq].cropmap[chain[res].seqid.num] is not None:
                                if not original_id:
                                    chain[res].seqid.num = inseq.imer[chseq].cropmap[chain[res].seqid.num]
                            else:
                                del chain[res]

    return instr
