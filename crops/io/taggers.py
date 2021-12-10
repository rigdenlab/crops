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
        "cropseq" : ".crops"+cut,
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
    nameseq=["", [], None] #PDB IF, chain IDs for seqgroup,
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

    if seqheader.startswith('>pdb|'):
        for i in range(5,len(seqheader)):
            if seqheader[i]=='|':
                for j in range(i+1,len(seqheader)):
                    if seqheader[j]!=' ' and seqheader[i]!='|':
                        newchid += seqheader[j]
                    elif seqheader[j]==' ':
                        if newchid != '':
                            nameseq[1].append(newchid)
                            newchid = ''
                    else:
                        if extrainfo:
                            return seqheader[i:]
                        break
            else:
                nameseq[0]+=seqheader[i]
        return nameseq

    for j in range(len(seqheader)):
        if seqheader[j]==">":
            idchar=True
        elif seqheader[j]==":" or seqheader[j]=="_":
            idchar=False
            namechar=True
        elif seqheader[j]==" ":
            pass
        elif seqheader[j]=="[":
            if seqheader[j:j+5]=="[auth" or seqheader[j:j+6]=="[ auth":
                newchid=''
        elif ((seqheader[j]=="a" and seqheader[j:j+4]=='auth') or
              (seqheader[j]=="u" and seqheader[j-1:j+3]=='auth') or
              (seqheader[j]=="t" and seqheader[j-2:j+2]=='auth') or
              (seqheader[j]=="h" and seqheader[j-3:j+1]=='auth')):
            pass
        elif seqheader[j]=="]":
            pass
        elif seqheader[j]=="|":
            if seqheader[j+1:j+6]=='Chain' or seqheader[j+1:j+6]=='chain':
                if newchid != '':
                    nameseq[2] = newchid
                k=0 if seqheader[j+6]==' ' else 1
                newchid=''
                for jj in range(j+6+k+1,len(seqheader)):
                    if seqheader[jj]==',':
                        nameseq[1].append(newchid)
                        newchid=''
                    elif seqheader[jj]==" ":
                        pass
                    elif seqheader[jj]=="[":
                        if seqheader[jj:jj+5]=="[auth" or seqheader[jj:jj+6]=="[ auth":
                            newchid=''
                    elif ((seqheader[jj]=="a" and seqheader[jj:jj+4]=='auth') or
                          (seqheader[jj]=="u" and seqheader[jj-1:jj+3]=='auth') or
                          (seqheader[jj]=="t" and seqheader[jj-2:jj+2]=='auth') or
                          (seqheader[jj]=="h" and seqheader[jj-3:jj+1]=='auth')):
                        pass
                    elif seqheader[jj]=="]":
                        pass
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
        elif seqheader[j]==" " and seqheader[j-1]!="_":
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


