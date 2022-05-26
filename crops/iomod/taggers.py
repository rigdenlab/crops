from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

import os
import logging

def target_format(inpath, terms=False, th=0):
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

    if os.path.basename(inpath) == 'pdb_chain_uniprot.csv':
        outcome = '|CROPS (UniProt via SIFTS)'
        if th > 0:
            outcome += ' - UniProt chain included threshold = ' + str(th)
    else:
        outcome = '|CROPS (Custom database)'
    if terms is True and th == 0:
        outcome += ' - Only terminals removed'

    return outcome

def infix_gen(inpath, terms=False):
    """Returns filename tag for outputs.

    :param inpath: Path to interval database used.
    :type inpath: str
    :param terms: Are terminal ends the only segments to be discarded?, defaults to False.
    :type terms: bool, optional
    :return: Filename tag.
    :rtype: str

    """
    if os.path.basename(inpath) == 'pdb_chain_uniprot.csv':
        cut = ".to_uniprot"
    else:
        cut = ".custom"

    if terms:
        cut = ".custom"

    infix_out = {
        "croprenum": ".crops" + cut,
        "cropseq": ".crops" + cut,
        "crop": ".crops.oldids" + cut,
        "renumber": ".crops.seq"}

    return infix_out

def retrieve_id(seqheader):
    """Extracts sequence IDs and additional comments from a standard .fasta header.

    :param seqheader: Standard .fasta header, starting with ">".
    :type seqheader: str
    :raises ValueError: If seqheader is not a string.
    :return: A list with the two sequence identifiers (e.g. [pdb ID, chain ID]) or a single string if extrainfo==True.
    :rtype: list [str], str

    """

    if not isinstance(seqheader, str):
        logging.critical('Argument is not a string.')
        raise ValueError
    headerinfo = {}
    headerinfo['mainid'] = ""
    headerinfo['chains'] = None
    headerinfo['seqid'] = None
    headerinfo['source'] = None
    headerinfo['comments'] = None
    namechar = False
    idchar = False
    newchid = ''
    # UniProt Swiss-Prot/TrEMBL
    if seqheader.startswith('>sp|') or seqheader.startswith('>tr|'):
        headerinfo['seqid'] = '1'
        if seqheader.startswith('>sp|'):
            headerinfo['source'] = 'UniProtKB/SwissProt'
        elif seqheader.startswith('>tr|'):
            headerinfo['source'] = 'UniProtKB/TrEMBL'
        for i in range(4, len(seqheader)):
            if seqheader[i] == '|':
                headerinfo['comments'] = seqheader[i+1:]
                break
            elif seqheader[i] == ' ' and seqheader[i:i+9] == ' archived':
                headerinfo['source'] += ' (archived)'
                headerinfo['comments'] = seqheader[i+1:]
                break
            else:
                headerinfo['mainid'] += seqheader[i]
        if headerinfo['chains'] is None:
            headerinfo['chains'] = set()
        headerinfo['chains'].add(headerinfo['mainid'])

    # UniRef
    elif seqheader.startswith('>UniRef'):
        headerinfo['seqid'] = '1'
        for i in range(1, len(seqheader)):
            if seqheader[i] == '_':
                headerinfo['source'] = seqheader[1:i]
                chi = i+1
            elif seqheader[i] == ' ':
                headerinfo['comments'] = seqheader[i+1:]
                headerinfo['mainid'] = seqheader[chi:i]
                if headerinfo['chains'] is None:
                    headerinfo['chains'] = set()
                headerinfo['chains'].add(headerinfo['mainid'])
                break

    # UniParc
    elif seqheader.startswith('>UPI'):
        headerinfo['seqid'] = '1'
        for i in range(4, len(seqheader)):
            headerinfo['source'] = 'UniParc'
            if seqheader[i] == ' ':
                headerinfo['comments'] = seqheader[i+1:]
                headerinfo['mainid'] = seqheader[1:i]
                if headerinfo['chains'] is None:
                    headerinfo['chains'] = set()
                headerinfo['chains'].add(headerinfo['mainid'])
                break

    # UniClust
    elif seqheader.startswith('>uc'):
        newchid = ''
        headerinfo['source'] = 'UniClust'
        tag = 0
        for i in range(1, len(seqheader)):
            newchid += seqheader[i]
            if seqheader[i] == '-':
                if tag == 0:
                    headerinfo['source'] = ('UniClust'+seqheader[3:i] + '_20' +
                                            seqheader[i+1:i+3] + '_' +
                                            seqheader[i+3:i+5])
                    tag = 1
                elif tag == 1:
                    ii = i
                    newchid = ''
                    while True:
                        ii += 1
                        if seqheader[ii] == '|':
                            headerinfo['seqid'] = newchid
                            headerinfo['comments'] = seqheader[ii+1:]
                            tag = 2
                            break
                        newchid += seqheader[ii]
            elif seqheader[i:i+16] == '|Representative=':
                newchid = ''
                ii = i + 16
                while True:
                    newchid += seqheader[ii]
                    ii += 1
                    if seqheader[ii] == ' ':
                        headerinfo['mainid'] = newchid
                        break
            elif seqheader[i:i+8] == 'Members=':
                newchid = ''
                headerinfo['chains'] = set()
                for ii in range(i+8, len(seqheader)):
                    if seqheader[ii] == ',' or ii == len(seqheader)-1:
                        headerinfo['chains'].add(newchid)
                        newchid = ''
                    else:
                        newchid += seqheader[ii]
                break

    # PDBe
    elif seqheader.startswith('>pdb|'):
        headerinfo['source'] = 'PDBe'
        newchid = ''
        tag = 'mainid'
        for i in range(5, len(seqheader)):
            if seqheader[i] == '|':
                headerinfo[tag] = newchid
                tag = 'chains'
                newchid = ''
            elif seqheader[i] == '_':
                headerinfo[tag] = newchid
                tag = 'seqid'
                newchid = ''
            elif (seqheader[i] == ' ' or
                  seqheader[i] == ',' or
                  i == len(seqheader)-1):
                if i == len(seqheader)-1:
                    newchid += seqheader[i]
                if headerinfo[tag] is None:
                    headerinfo[tag] = set()
                headerinfo[tag].add(newchid)
                newchid = ''
            else:
                if tag == 'mainid':
                    newchid += seqheader[i].lower()
                else:
                    newchid += seqheader[i]

    else:
        # RCSB PDB, CROPS, MrBUMP, others
        if seqheader.startswith('>crops|'):
            seqheader = ">" + seqheader[7:]
            headerinfo['source'] = 'CROPS'
        else:
            headerinfo['source'] = 'RCSB PDB'
        for j in range(len(seqheader)):
            if seqheader[j] == ">":
                idchar = True
            elif seqheader[j] == ":" or seqheader[j] == "_":
                idchar = False
                namechar = True
            elif seqheader[j] == " ":
                if seqheader[j-2] == "_" and newchid != '':
                    try:
                        int(newchid)
                        headerinfo['seqid'] = newchid
                        newchid = ''
                    except Exception:
                        if headerinfo['chains'] is None:
                            headerinfo['chains'] = set()
                        headerinfo['chains'].add(newchid)
                        headerinfo['comments'] = seqheader[j+1:]
                        if seqheader[j+1:j+12] == "resolution:":
                            headerinfo['source'] = 'MrBUMP'
                        break
                else:
                    pass
            elif seqheader[j] == "[":
                if seqheader[j:j+5] == "[auth" or seqheader[j:j+6] == "[ auth":
                    newchid = ''
            elif ((seqheader[j] == "a" and seqheader[j:j+4] == 'auth') or
                  (seqheader[j] == "u" and seqheader[j-1:j+3] == 'auth') or
                  (seqheader[j] == "t" and seqheader[j-2:j+2] == 'auth') or
                  (seqheader[j] == "h" and seqheader[j-3:j+1] == 'auth')):
                pass
            elif seqheader[j] == "]":
                pass
            elif seqheader[j] == "|":
                if idchar is True:
                    idchar = False
                if seqheader[j+1:j+6].lower() == 'chain':
                    if newchid != '':
                        headerinfo['seqid'] = newchid
                    k = 0 if seqheader[j+6] == ' ' else 1
                    newchid = ''
                    for jj in range(j+6+k+1, len(seqheader)):
                        if seqheader[jj] == ',':
                            if headerinfo['chains'] is None:
                                headerinfo['chains'] = set()
                            headerinfo['chains'].add(newchid)
                            newchid = ''
                        elif seqheader[jj] == " ":
                            pass
                        elif seqheader[jj] == "[":
                            if (seqheader[jj:jj+5] == "[auth" or
                                    seqheader[jj:jj+6] == "[ auth"):
                                newchid = ''
                        elif ((seqheader[jj] == "a" and seqheader[jj:jj+4] == 'auth') or
                              (seqheader[jj] == "u" and seqheader[jj-1:jj+3] == 'auth') or
                              (seqheader[jj] == "t" and seqheader[jj-2:jj+2] == 'auth') or
                              (seqheader[jj] == "h" and seqheader[jj-3:jj+1] == 'auth')):
                            pass
                        elif seqheader[jj] == "]":
                            pass
                        elif (seqheader[jj] == "|" or seqheader[jj] == ":" or
                              jj == len(seqheader)-1):
                            if jj == len(seqheader)-1:
                                newchid += seqheader[jj]
                            else:
                                headerinfo['comments'] = seqheader[jj+1:]
                            if headerinfo['chains'] is None:
                                headerinfo['chains'] = set()
                            headerinfo['chains'].add(newchid)
                            newchid = ''
                            namechar = False
                            return headerinfo
                        else:
                            newchid += seqheader[jj]
                else:
                    if namechar is True:
                        headerinfo['comments'] = seqheader[j+1:]
                        if headerinfo['chains'] is None:
                            headerinfo['chains'] = set()
                        headerinfo['chains'].add(newchid)
                    else:
                        pass
                    break
            elif seqheader[j] == " " and seqheader[j-1] != "_":
                headerinfo['comments'] = seqheader[j+1:]
                if namechar is True:
                    if headerinfo['chains'] is None:
                        headerinfo['chains'] = set()
                    headerinfo['chains'].add(newchid)
                    break
            else:
                if namechar is True:
                    newchid += seqheader[j]
                elif idchar is True:
                    headerinfo['mainid'] += seqheader[j].lower()

    return headerinfo

def makeheader(mainid=None, seqid=None, chains=None,
               source=None, extrainfo=None, short=False):
    """Returns a fasta header of the format ">MainID_seqID|Chains chain list|extrainfo".

    :param mainid: PDB ID, or Uniprot ID, etc.
    :type mainid: str
    :param seqid: Sequence Identifier, usually a natural number: "1", "2", etc, defaults to None.
    :type seqid: str, optional
    :param chains: A set containing the chain IDs of monomers sharing the same sequence, defaults to None.
    :type chains: set [str], optional
    :param source: A set containing the chain IDs of monomers sharing the same sequence, defaults to None.
    :type source: str, optional
    :param extrainfo: Additional information to be included in the header, defaults to None.
    :type extrainfo: str, optional
    :param short: If selected, it returns a header of format '>main_seq|Chains A,B,C', defaults to False.
    :type extrainfo: bool, optional
    :raises ValueError: If any of mainid, seqid, extrainfo or elements of chains are not strings, or chains is not a set.
    :return: A fasta header.
    :rtype: str

    """
    if not isinstance(mainid, str):
        try:
            mainid = str(mainid)
        except Exception:
            logging.critical("Argument 'mainid' is not a string.")
            raise ValueError

    if short is True:
        newheader = '>' + mainid.upper()
    else:
        newheader = '>crops|' + mainid
    if seqid is not None:
        if not isinstance(seqid, str):
            try:
                seqid = str(seqid)
            except Exception:
                logging.critical("Argument 'seqid' is not a string.")
                raise ValueError
        newheader += '_'
        newheader += seqid
    newheader += '|'
    if chains is not None:
        if not isinstance(chains, set):
            logging.critical("Argument 'chains' is not a set of strings.")
            raise ValueError
        else:
            for element in chains:
                if not isinstance(element, str):
                    logging.critical("Argument 'chains' is not a set of strings.")
                    raise ValueError

        if len(chains) > 0:
            if len(chains) > 1:
                newheader += 'Chains '
            else:
                newheader += 'Chain '
            for element in sorted(chains):
                newheader += element
                newheader += ','
            newheader = newheader[:-1]
        if not short:
            newheader += '|'
    if not short:
        if source is None:
            source = 'Unknown'
        if not isinstance(source, str):
            try:
                source = str(source)
            except Exception:
                logging.critical("Argument 'source' is not a string.")
                raise ValueError
        newheader += 'Source: '
        newheader += source
        newheader += '|'
    if extrainfo is not None and not short:
        if not isinstance(extrainfo, str):
            try:
                extrainfo = str(extrainfo)
            except Exception:
                logging.critical('Argument extrainfo is not a list of strings.')
                raise ValueError
        newheader += extrainfo

    return newheader
