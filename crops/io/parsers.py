from crops import __prog__, __description__, __author__, __date__, __version__

import gemmi
import os
import csv
import urllib3
import copy

from crops.elements.sequences import oligoseq
from crops.elements.sequences import sequence
from crops.elements.sequences import guess_type
from crops.io.taggers import retrieve_id
from crops.elements.intervals import intinterval

def import_db(inpath, pdb_in=None):
    """Imports intervals database. Input must be a .csv file (filepath).
    If imported file is not 'pdb_chain_uniprot.csv' from SIFTS database,
    the columns must contain molecule ID, chain ID, lower element of subset,
    and higher element of subset, in this order.

    :param inpath: Path to interval database used.
    :type inpath: str
    :param pdb_in: Chain ID(s). If given, the imported values
        will be filtered to contain only IDs provided, defaults to None.
    :type pdb_in: str, set, dict, optional
    :raises TypeError: When pdb_in is given and is neither a string nor a dictionary.
    :return: dict [str, :class:`~crops.elements.intervals.intinterval`im]
    :rtype: A dictionary of :class:`~crops.elements.intervals.intinterval`.

    """
    database_out = {}
    if isinstance(pdb_in, str) is True:
        pdb_in_lower = set()
        pdb_in_lower.add(pdb_in.lower())
    elif isinstance(pdb_in, dict) is True or isinstance(pdb_in, set) is True:
        pdb_in_lower = set()
        for key in pdb_in:
            if isinstance(key, str) is False:
                raise TypeError('Argument should be either None, a string, ' +
                                'a set, or a dictionary with empty values.')
            pdb_in_lower.add(key.lower())
    elif pdb_in is None:
        pass
    else:
        raise TypeError('Argument should be either None, a string, ' +
                        'a set, or a dictionary with empty values.')

    if os.path.basename(inpath) == 'pdb_chain_uniprot.csv':
        mol = 0
        chain = 1
        up = 2
        leftend = 3
        rightend = 4
    else:
        mol = 0
        chain = 1
        leftend = 2
        rightend = 3
        up = None

    csv_chain_file = open(inpath)
    csv_chain = csv.reader(csv_chain_file)
    for entry in csv_chain:
        if entry[0][0] != "#" and entry[0] != "PDB":
            if pdb_in is None or entry[mol].lower() in pdb_in_lower:
                if entry[mol].lower() not in database_out:
                    database_out[entry[mol].lower()] = {}
                if entry[chain] not in database_out[entry[mol].lower()]:
                    database_out[entry[mol].lower()][entry[chain]] = intinterval(description=entry[mol].lower()+'_'+entry[chain])
                    if up is not None:
                        database_out[entry[mol].lower()][entry[chain]].tags['uniprot'] = {}
                database_out[entry[mol].lower()][entry[chain]] = \
                    database_out[entry[mol].lower()][entry[chain]].union(other=[int(entry[leftend]), int(entry[rightend])])
                if up is not None:
                    if entry[up].upper() not in database_out[entry[mol].lower()][entry[chain]].tags['uniprot']:
                        database_out[entry[mol].lower()][entry[chain]].tags['uniprot'][entry[up]] = \
                            intinterval(description=entry[up].upper())
                    database_out[entry[mol].lower()][entry[chain]].tags['uniprot'][entry[up]] = \
                        database_out[entry[mol].lower()][entry[chain]].tags['uniprot'][entry[up]].union([int(entry[leftend]), int(entry[rightend])])
    return database_out


def parsestrfile(str_inpath):
    """Returns dictionary containing :class:`~gemmi.Structure` objects and another one with the file names.

    :param str_inpath: Either a directory or file path.
    :type str_inpath: str
    :raises KeyError: More than one structure file containing same identifier.
    :return strdict: A dictionary containing imported :class:`~gemmi.Structure` objects.
    :rtype strdict: dict [str, :class:`~gemmi.Structure`]
    :return filedict: A dictionary containing file names.
    :rtype filedict: dict [str, str]

    """
    strdict = {}
    filedict = {}
    if os.path.isfile(str_inpath):

        structure = gemmi.read_structure(str_inpath)
        pdbid = structure.name.lower()
        strdict[pdbid] = structure
        filedict[pdbid] = os.path.basename(str_inpath)
    elif os.path.isdir(str_inpath):
        filelist = os.listdir(str_inpath)
        for file in filelist:
            if os.isfile(file):
                try:
                    structure = gemmi.read_structure(file)
                    pdbid = structure.name.lower()
                    if pdbid in strdict:
                        raise KeyError('Structure '+pdbid+' loaded more ' +
                                       'than once. Check files in directory ' +
                                       'and remove duplicates.')
                    strdict[pdbid] = structure
                    filedict[pdbid] = os.path.basename(str_inpath)
                except:
                    pass

    return strdict, filedict

def parseseqfile(inpath, uniprot=None):
    """Sequence file parser.

    :param inpath: Sequence file path.
    :type inpath: str
    :param uniprot: A dictionary or set of Uniprot codes, defaults to None.
    :type uniprot: str, dict [str, any], optional
    :return: A dictionary containing parsed :class:`~crops.elements.sequences.oligoseq` objects.
    :rtype: dict [str, :class:`~crops.elements.sequences.oligoseq`]

    """
    newseqs = {}
    newid = []
    head = ''
    chain = ''
    ignore = False

    if uniprot is not None:
        if (not isinstance(uniprot, str) and not isinstance(uniprot, dict) and
                not isinstance(uniprot, set)):
            raise TypeError('Input argument uniprot must be one of ' +
                            'a string, a set, or a dictionary.')
        elif isinstance(uniprot, str):
            unitemp = uniprot
            uniprot = set()
            uniprot.add(unitemp)
        for upcode in uniprot:
            if not isinstance(upcode, str):
                raise TypeError('Input argument uniprot must be one of ' +
                                'a string, a set, or a dictionary.')

    if inpath == 'server-only' and uniprot is None:
        raise TypeError("Input argument inpath cannot be 'server-only' " +
                        "when a set or dict of uniprot ids is not provided.")
    elif inpath == 'server-only' and uniprot is not None:
        pass
    else:
        with open(inpath, 'r') as f:
            indx = -1
            while True:
                line = f.readline().rstrip()
                if (not line or line.startswith(">")) and not ignore:
                    if uniprot is not None:
                        if indx >= 0:
                            if newid['mainid'].upper() not in newseqs:
                                newseqs[newid['mainid']] = oligoseq(oligomer_id=newid['mainid'])
                            aseq = sequence(seqid=newid['seqid'],
                                            oligomer=newid['mainid'].upper(),
                                            seq=chain, chains=newid['chains'],
                                            source=newid['source'],
                                            header=head, extrainfo=newid['comments'])
                            newseqs[newid['mainid']].add_sequence(aseq)
                            if len(newseqs) == len(uniprot):
                                break
                    else:
                        if indx >= 0:
                            if newid['mainid'].lower() not in newseqs:
                                newseqs[newid['mainid'].lower()] = oligoseq(oligomer_id=newid['mainid'].lower())
                            aseq = sequence(seqid=newid['seqid'],
                                            oligomer=newid['mainid'].upper(),
                                            seq=chain, chains=newid['chains'],
                                            source=newid['source'],
                                            header=head, extrainfo=newid['comments'])
                            newseqs[newid['mainid'].lower()].add_sequence(aseq)
                    if not line:
                        try:
                            line = f.readline().rstrip()
                            if not line:
                                break
                        except:
                            break
                if line.startswith(">"):
                    newid = retrieve_id(line)
                    head = line
                    indx += 1
                    chain = ''
                    if uniprot is not None:
                        ignore = False if newid['mainid'] in uniprot else True

                elif line.startswith("#") or line.startswith(' #'):
                    pass
                else:
                    if not ignore:
                        chain += str(line)

    if inpath == 'server-only' or uniprot is not None:
        for upcode in uniprot:
            if upcode.upper() not in newseqs:
                try:
                    for line in urllib3.urlopen('https://www.uniprot.org/uniprot/' +
                                                upcode.upper()+'.fasta'):
                        if line.startswith(">"):
                            chain = ''
                            head = line
                            newid = retrieve_id(line)
                        else:
                            chain += str(line)
                except:
                    if inpath == 'server-only':
                        msg = ('Uniprot sequence ' + upcode.upper() +
                               ' not found online. Check your internet connexion.')
                    else:
                        msg = ('Uniprot sequence ' + upcode.upper() +
                               ' not found in local file or online. Check your internet connexion.')
                    raise OSError(msg)
                if upcode.upper() not in newseqs:
                    newseqs[newid['mainid']] = oligoseq(oligomer_id=newid['mainid'])
                aseq = sequence(seqid=newid['seqid'],
                                oligomer=newid['mainid'].upper(),
                                seq=chain, chains=newid['chains'],
                                source=newid['source'],
                                header=head, extrainfo=newid['comments'])
                newseqs[newid['mainid']].add_sequence(aseq)

    return newseqs

def parsemapfile(inpath):
    """Cropmap file parser.

    :param inpath: Cropmap file path.
    :type inpath: str
    :return: A dictionary containing parsed mapping and backmapping coordinates.
    :rtype: dict [str, dict[str, dict[str, dict[int, int]]]]

    """
    mapdict = {}
    newid = []
    with open(inpath, 'r') as f:
        indx = -1
        while True:
            line = f.readline().rstrip()
            if (not line or line.startswith(">")):
                if indx >= 0:
                    if newid['mainid'].lower() not in mapdict:
                        mapdict[newid['mainid'].lower()] = {}
                    if newid['seqid'] not in mapdict[newid['mainid'].lower()]:
                        mapdict[newid['mainid'].lower()][newid['seqid']] = {}
                        mapdict[newid['mainid'].lower()][newid['seqid']]['cropmap'] = copy.deepcopy(forthmap)
                        mapdict[newid['mainid'].lower()][newid['seqid']]['cropbackmap'] = copy.deepcopy(backmap)
                if not line:
                    try:
                        line = f.readline().rstrip()
                        if not line:
                            break
                    except:
                        break

            if line.startswith(">"):
                newid = retrieve_id(line)
                indx += 1
                forthmap = {}
                backmap = {}
            elif line.startswith("#") or line.startswith(' #'):
                pass
            else:
                m = line.split('  ')
                if m[1] != '0':
                    forthmap[int(m[0])] = int(m[1])
                    backmap[int(m[1])] = int(m[0])
                else:
                    forthmap[int(m[0])] = None

    return mapdict
