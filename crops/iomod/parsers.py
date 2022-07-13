"""This is CROPS: Cropping and Renumbering Operations for PDB structure and Sequence files."""

from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

import gemmi
import os
import csv
from urllib import request as ur
import copy
import logging

from crops.elements.sequences import oligoseq
from crops.elements.sequences import sequence
from crops.elements.sequences import guess_type
from crops.iomod.taggers import retrieve_id
from crops.elements.intervals import intinterval


def parse_db(instream, pdbset=None):
    """Import intervals database from csv-formatted string.

    If imported file is not 'pdb_chain_uniprot.csv' from SIFTS database,
    the columns must contain, in this order, molecule ID, chain ID,
    lower element of subset, and higher element of subset, in this order.
    More than one row with the same molecule and chain IDs are used to indicate
    a discontinuous interval with more than one subset.

    :param instream: Interval database, csv-formatted string.
    :type instream: str
    :param pdbset: Molecule IDs to return, if None it returns them all, defaults to None.
    :type pdbset: str or set or dict, optional

    :raises TypeError: When pdbset is given and is not one of a string, set or dictionary.
        It will also raise this error when the database is not from SIFTS or a minimal file (4 elements per line).

    :return: A dictionary containing :class:`crops.elements.intervals.intinterval` objects.
    :rtype: dict [str, :class:`crops.elements.intervals.intinterval`]

    """
    database_out = {}
    if isinstance(pdbset, str) is True:
        pdb_in_upper = set()
        pdb_in_upper.add(pdbset.upper())
    elif isinstance(pdbset, dict) is True or isinstance(pdbset, set) is True:
        pdb_in_upper = set()
        for key in pdbset:
            if isinstance(key, str) is False:
                logging.critical('Argument pdbset should be either None, a string, '
                                 'a set, or a dictionary with empty values.')
                raise TypeError
            pdb_in_upper.add(key.upper())
    elif pdbset is None:
        pdb_in_upper = set()
    else:
        logging.critical('Argument pdbset should be either None, a string, '
                         'a set, or a dictionary with empty values.')
        raise TypeError

    csv_list = instream.splitlines()
    csv_chain = csv.reader(csv_list)
    if csv_list[-1].count(',') == 8:
        mol = 0
        chain = 1
        up = 2
        leftend = 3
        rightend = 4
    elif csv_list[-1].count(',') == 3:
        mol = 0
        chain = 1
        leftend = 2
        rightend = 3
        up = None
    else:
        logging.critical('Database format is neither SIFTS-like nor minimal.')
        raise TypeError

    csv_chain = csv.reader(instream.splitlines())
    for entry in csv_chain:
        if entry[0][0] != "#" and entry[0] != "PDB":
            molid = entry[mol].upper()
            if pdbset is None or molid in pdb_in_upper:
                if molid not in database_out:
                    database_out[molid] = {}
                if entry[chain] not in database_out[molid]:
                    database_out[molid][entry[chain]] = intinterval(description=molid+'_'+entry[chain])
                    if up is not None:
                        database_out[molid][entry[chain]].tags['uniprot'] = {}
                database_out[molid][entry[chain]] = \
                    database_out[molid][entry[chain]].union(other=[int(entry[leftend]), int(entry[rightend])])
                if up is not None:
                    upid = entry[up].upper()
                    if upid not in database_out[molid][entry[chain]].tags['uniprot']:
                        database_out[molid][entry[chain]].tags['uniprot'][upid] = \
                            intinterval(description=upid)
                    database_out[molid][entry[chain]].tags['uniprot'][upid] = \
                        database_out[molid][entry[chain]].tags['uniprot'][upid].union([int(entry[leftend]), int(entry[rightend])])

    return database_out


def import_db(inpath, pdb_in=None):
    """Import intervals database from a .csv file.

    If the imported file is not 'pdb_chain_uniprot.csv' from the SIFTS database,
    the .csv file must be formatted as follows: four columns containing molecule ID, chain ID, lower element of subset,
    and higher element of subset, in this order.

    :param inpath: Path to the interval database to be imported.
    :type inpath: str
    :param pdb_in: Molecule IDs to return, if None it returns them all, defaults to None.
    :type pdb_in: str or set or dict, optional

    :return: Parsed interval database.
    :rtype: dict [str,  :class:`crops.elements.intervals.intinterval`]

    """
    with open(inpath, 'r') as f:
        csv_stream = f.read()

    database_out = parse_db(csv_stream, pdbset=pdb_in)

    return database_out


def parsestr(instream):
    """Parse structure file from a string.

    :param instream: Imported-to-string structure file.
    :type instream: str

    :return: Parsed structure
    :rtype: :obj:`gemmi.Structure`

    """
    strout = gemmi.read_pdb_string(instream)

    return strout


def parsestrfile(str_input, intype='path'):
    """Parse structure file(s).

    :param str_input: Either a directory or file path or a structure in string format.
    :type str_input: str
    :param intype: One of 'path' or 'string', defaults to 'path'.
    :type intype: str, optional

    :raises KeyError: If more than one structure file contains the same identifier.
    :raises ValueError: If the argument 'intype' has an invalid value.

    :return strdict: A dictionary containing parsed structures.
    :rtype strdict: dict [str, :obj:`gemmi.Structure`]
    :return filedict: A dictionary containing the structure file name(s).
    :rtype filedict: dict [str, str]

    """
    strdict = {}
    filedict = {}
    if intype == 'string':
        structure = parsestr(str_input)
        pdbid = structure.name.upper()
        strdict[pdbid] = structure
        filedict[pdbid] = None
    elif intype == 'path':
        if os.path.isfile(str_input):
            with open(str_input, 'r') as f:
                strfile = f.read()
            structure = parsestr(strfile)
            if os.path.splitext(os.path.basename(str_input))[1].lower() == '.pdb':
                structure.name = os.path.splitext(os.path.basename(str_input))[0]
            else:
                structure.name = os.path.basename(str_input)
            pdbid = structure.name.upper()
            strdict[pdbid] = structure
            filedict[pdbid] = os.path.basename(str_input)
        elif os.path.isdir(str_input):
            filelist = os.listdir(str_input)
            for file in filelist:
                if os.isfile(file):
                    try:
                        with open(str_input, 'r') as f:
                            strfile = f.read()
                        structure = parsestr(strfile)
                        pdbid = structure.name.upper()
                        if pdbid in strdict:
                            logging.critical('Structure ' + pdbid + ' loaded more '
                                             'than once. Check files in directory '
                                             'and remove duplicates.')
                            raise KeyError
                        strdict[pdbid] = structure
                        filedict[pdbid] = os.path.basename(str_input)
                    except Exception:
                        logging.warning("There was some error while processing '" +
                                        pdbid + "'. Ignoring structure.")
                        pass
    else:
        logging.critical("Invalid value for argument 'intype'")
        raise ValueError

    return strdict, filedict


def parseseq(instream, inset=None):
    """Parse sequence(s).

    :param instream: Imported-to-string sequence file content (fasta format).
    :type instream: str
    :param inset: Sequence IDs to return, if None it returns them all, defaults to None.
    :type inset: set or dict or str, optional

    :raises TypeError: When inset a set [str]; or instream is not a string.

    :return: Parsed sequences.
    :rtype: dict [str, :class:`crops.elements.sequences.oligoseq`]

    """
    if isinstance(instream, str) is False:
        logging.critical('Input argument instream should be a string.')
        raise TypeError

    if inset is not None:
        if (not isinstance(inset, str) and not isinstance(inset, dict) and
                not isinstance(inset, set)):
            logging.critical('Input argument inset should be a set or, '
                             'alternatively a string or a dictionary.')
            raise TypeError
        elif isinstance(inset, str):
            temp = inset
            inset = set()
            inset.add(temp)
        upperset = set()
        for element in inset:
            if not isinstance(element, str):
                logging.critical('Elements in inseq should be strings.')
                raise TypeError
            upperset.add(element.upper())

    newseqs = {}
    newid = []
    head = ''
    chain = ''
    ignore = False
    ignore = False
    indx = -1
    inseqlines = instream.splitlines()
    inseqlines.append('')
    for raw in range(len(inseqlines)):
        line = inseqlines[raw].rstrip()
        if (not line or line.startswith(">")) and not ignore:
            if indx >= 0:
                if newid['mainid'] not in newseqs:
                    newseqs[newid['mainid']] = oligoseq(oligomer_id=newid['mainid'])
                aseq = sequence(seqid=newid['seqid'],
                                oligomer=newid['mainid'],
                                seq=chain, chains=newid['chains'],
                                source=newid['source'],
                                header=head, extrainfo=newid['comments'])
                newseqs[newid['mainid']].add_sequence(aseq)
        else:
            pass

        if line.startswith(">"):
            newid = retrieve_id(line)
            head = line
            indx += 1
            chain = ''
            if inset is not None:
                ignore = False if newid['mainid'] in upperset else True
        elif line.startswith("#") or line.startswith(' #'):
            continue
        else:
            if not ignore:
                chain += str(line)

    return newseqs


# TODO: ADD use_PDBserver argument
def parseseqfile(seq_input='server-only', inset=None, use_UPserver=False):
    """Parse sequence file containing one or more sequences.

    If 'server-only' is inserted instead of a local file name,

    :param seq_input: Sequence file path, defaults to 'server-only'.
    :type seq_input: str, optional
    :param inset: Sequence IDs to return, if None it returns them all, defaults to None.
    :type inset: set or dict or str, optional
    :param use_UPserver: Use UniProt server as a backup for those ids not found in `seq_input` (all of them if `seq_input` == 'server-only'), defaults to False.
    :type use_UPserver: bool, optional

    :raises TypeError: If `inset` is not a str or set [str] or dict [str, str].
    :raises ValueError: If seq_input`=='server-only' but `use_UPserver` is False or `inset` is None.

    :return: Parsed sequences.
    :rtype: dict [str, :class:`crops.elements.sequences.oligoseq`]

    """
    if ((seq_input == 'server-only' and use_UPserver is False) or
            (seq_input == 'server-only' and inset is None)):
        logging.critical("Input argument seq_input cannot be 'server-only' "
                         "if no inset is given and use server is False.")
        raise ValueError

    if inset is not None:
        if (not isinstance(inset, str) and not isinstance(inset, dict) and
                not isinstance(inset, set)):
            logging.critical('Input argument inset should be a set or, '
                             'alternatively a string or a dictionary.')
            raise TypeError
        elif isinstance(inset, str):
            temp = inset
            inset = set()
            inset.add(temp)
        upperset = set()
        for element in inset:
            if not isinstance(element, str):
                logging.critical('Elements in inseq should be strings.')
                raise TypeError
            upperset.add(element.upper())
    else:
        upperset = None

    if seq_input != 'server-only':
        with open(seq_input, 'r') as f:
            inseq = f.read()
        newseqs = parseseq(inseq, inset=upperset)
    else:
        newseqs = {}

    inseq = ''
    if use_UPserver:
        for upcode in upperset:
            if upcode not in newseqs:
                try:
                    download = ur.urlopen('https://www.uniprot.org/uniprot/' +
                                          upcode.upper() + '.fasta')
                    inseq += download + os.linesep
                except Exception:
                    if seq_input == 'server-only':
                        msg = ('Uniprot sequence ' + upcode +
                               ' not found online. Check your internet connexion.')
                    else:
                        msg = ('Uniprot sequence ' + upcode +
                               ' not found either in local file or online. '
                               'Check your internet connexion.')
                    logging.warning(msg)
                    continue
        if inseq != '':
            upseqs = parseseq(inseq)
            newseqs.update(upseqs)

    return newseqs


def parsemap(instream):
    """Parse cropmap from string.

    :param instream: Imported-to-string cropmap file content.
    :type instream: str

    :return: Mapping and backmapping coordinates.
    :rtype: dict [str, dict [str, dict [str, dict [int, int]]]]

    """
    mapdict = {}
    newid = []
    indx = -1
    inmaplines = instream.splitlines()
    inmaplines.append('')
    for raw in range(len(inmaplines)):
        line = inmaplines[raw].rstrip()
        if (not line or line.startswith(">")):
            if indx >= 0:
                if newid['mainid'] not in mapdict:
                    mapdict[newid['mainid']] = {}
                if newid['seqid'] not in mapdict[newid['mainid']]:
                    mapdict[newid['mainid']][newid['seqid']] = {}
                    mapdict[newid['mainid']][newid['seqid']]['cropmap'] = copy.deepcopy(forthmap)
                    mapdict[newid['mainid']][newid['seqid']]['cropbackmap'] = copy.deepcopy(backmap)
            if not line:
                try:
                    line = f.readline().rstrip()
                    if not line:
                        break
                except Exception:
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


def parsemapfile(input_map):
    """Cropmap file parser.

    :param input_map: Cropmap file path.
    :type input_map: str

    :return: Mapping and backmapping coordinates.
    :rtype: dict [str, dict[str, dict[str, dict[int, int]]]]

    """
    with open(input_map, 'r') as f:
        infile = f.read()

    mdict = parsemap(infile)

    return mdict
