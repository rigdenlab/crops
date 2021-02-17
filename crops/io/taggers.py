from crops.about import __prog__, __description__, __author__, __date__, __version__

import os

def target_format(inpath,terms=False,th=0):
    """Returns extra information for .fasta headers.

    :param inpath: Path to interval database used.
    :type inpath: str
    :param terms: Are only terminal ends discarded?, defaults to False.
    :type terms: bool, optional
    :param th: Uniprot threshold, defaults to 0.
    :type th: int, float, optional
    :return: Extra information for .fasta headers
    :rtype: str

    """

    if os.path.basename(inpath)=='pdb_chain_uniprot.csv':
        outcome=' | CROPS | UNIPROT via SIFTS'
        if th>0:
            outcome += ' - UNIPROT CHAIN INCLUDED THRESHOLD = '+str(th)
    else:
        outcome=' | CROPS | CUSTOM'
    if terms and th==0:
        outcome += ' - ONLY TERMINALS REMOVED'

    return outcome

def infix_gen(inpath,terms=False):
    """Returns filename tag for outputs.

    :param inpath: Path to interval database used.
    :type inpath: str
    :param terms: Are terminal ends the only segments to be discarded?, defaults to False.
    :type terms: bool, optional
    :return: Filename tag.
    :rtype: str

    """
    if os.path.basename(inpath)=='pdb_chain_uniprot.csv':
        cut=".to_uniprot"
    else:
        cut=".custom"

    if terms:
        cut=".custom"

    infix_out={
        "croprenum" : ".crops"+cut,
        "crop" : ".crops.oldids"+cut,
        "renumber" : ".crops.seq"}

    return infix_out

def retrieve_id(seqheader,extrainfo=False):
    """Extracts sequence IDs from a standard .fasta header.

    :param seqheader: Standard .fasta header, starting with ">".
    :type seqheader: str
    :param extrainfo: If True, extra information string is returned instead of sequence IDs, defaults to False.
    :type extrainfo: bool, optional
    :raises ValueError: If seqheader is not a string.
    :return: A list with the two sequence identifiers (e.g. [pdb ID, chain ID]) or a single string if extrainfo==True.
    :rtype: list [str], str

    """

    if not isinstance(seqheader,str):
        raise ValueError('Argument is not a str')
    namechar=False
    idchar=False
    nameseq=['',[]]
    newchid=''
    if seqheader.startswith('>sp|'):
        for i in range(4,len(seqheader)):
            if seqheader[i]=='|':
                if extrainfo:
                    return seqheader[i:]
                break
            else:
                nameseq[0]+=seqheader[i]
        nameseq[1]=[nameseq[0]]
        return nameseq

    for j in range(len(seqheader)):
        if seqheader[j]==">":
            idchar=True
        elif seqheader[j]==":" or seqheader[j]=="_":
            idchar=False
            namechar=True
        elif seqheader[j]=="|":
            if seqheader[j+1:j+6]=='Chain' or seqheader[j+1:j+6]=='chain':
                k=0 if seqheader[j+6]==' ' else 1
                newchid=''
                for jj in range(j+6+k+1,len(seqheader)):
                    if seqheader[jj]==',':
                        nameseq[1].append(newchid)
                        newchid=''
                    elif seqheader[jj]=="|" or seqheader[jj]==":" or jj==len(seqheader)-1:
                        if extrainfo:
                            if jj==len(seqheader)-1:
                                return ''
                            else:
                                return seqheader[jj:]
                        if jj==len(seqheader)-1:
                            newchid+=seqheader[jj]
                        nameseq[1].append(newchid)
                        newchid=''
                        return nameseq
                    else:
                        newchid+=seqheader[jj]
            else:
                if extrainfo:
                    return seqheader[j:]
                nameseq[1].append(newchid)
                return nameseq
        elif seqheader[j]==" ":
            if extrainfo:
                return seqheader[j:]
            if namechar:
                nameseq[1].append(newchid)
                return nameseq
        else:
            if namechar:
                newchid += seqheader[j]
            elif idchar:
                nameseq[0] += seqheader[j].lower()


