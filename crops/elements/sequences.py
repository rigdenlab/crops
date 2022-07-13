"""This is CROPS: Cropping and Renumbering Operations for PDB structure and Sequence files."""

from crops import __prog__, __description__, __author__
from crops import __date__, __version__, __copyright__

import os
import io
import copy
import logging

from crops.iomod.taggers import retrieve_id
from crops.iomod.taggers import makeheader
from crops.libs.rescodes import reslist
from crops.libs.rescodes import nuclist
from crops.elements.intervals import intinterval


def guess_type(inseq):
    """Return the biological type of the sequence as guessed from residue types.

    :param inseq: Sequence to be evaluated.
    :type inseq: str

    :return: Sequence type ('Protein' or 'DNA' or 'RNA' or 'Unknown').
    :rtype: str

    """
    if not isinstance(inseq, str):
        logging.critical("Sequence 'inseq' should be a string.")
        raise TypeError

    outtype = None
    for char in inseq:
        if (char == 'C' or char == 'A' or char == 'G' or char == 'I' or
                char == 'X' or char == '-' or char == '+' or char == '*'):
            pass
        elif char not in nuclist.values():
            if char in reslist.values():
                outtype = 'Protein'
            else:
                outtype = 'Unknown'
        else:
            if char == 'T':
                if outtype == 'DNA' or outtype == 'Protein':
                    pass
                elif outtype is None:
                    outtype = 'DNA'
                elif outtype == 'RNA':
                    outtype = 'Protein'
            elif char == 'U':
                if outtype == 'RNA' or outtype == 'Protein':
                    pass
                elif outtype is None:
                    outtype = 'RNA'
                elif outtype == 'DNA':
                    outtype = 'Protein'

    if outtype is None:
        outtype = 'DNA or RNA'

    return outtype


class sequence:
    """A :class:`crops.elements.sequences.sequence` object representing a single chain sequence.

    The :class:`crops.elements.sequences.sequence` class represents a data structure to hold all
    sequence versions and other useful information characterising it.
    It contains functions to store, manipulate and organise sequence versions.

    :param seqid: Sequence identifier. Can be used alone or together with oligomer ID, defaults to None.
    :type seqid: str
    :param oligomer: Oligomer identifier. Sometimes as important as seqid, defaults to None.
    :type oligomer: str, optional
    :param seq: Sequence string, defaults to None.
    :type seq: str, optional
    :param chains: The names of chains having this sequence, defaults to None.
    :type chains: set [str], optional
    :param source: Source of the sequence, defaults to None
    :type source: str, optional
    :param header: Standard .fasta header, starting with ">", defaults to None.
    :type header: str, optional
    :param biotype: Type of molecule ('Protein', 'DNA', 'RNA'...), defaults to None.
    :type biotype: str, optional
    :param extrainfo: Other useful information about the sequence, defaults to None.
    :type extrainfo: str, optional

    :ivar name: Sequence identifier.
    :vartype name: str
    :ivar oligomer_id: Oligomer identifier.
    :vartype oligomer_id: str
    :ivar chains: The names of chains having this sequence.
    :vartype chains: set [str]
    :ivar seqs: The set of sequences, including default "mainseq".
    :vartype seqs: dict [str, str]
    :ivar source: Source of the sequence.
    :vartype source: str
    :ivar source_headers: A list of headers from input files.
    :vartype source_headers: list [str]
    :ivar crops_header: A new header containing the information from the object that will be used when printing sequence and cropmap.
    :vartype crops_header: str
    :ivar biotype: Type of molecule ('Protein', 'DNA', 'RNA'...).
    :vartype biotype: str
    :ivar infostring: Other useful information about the sequence.
    :vartype infostring: str
    :ivar cropmap: A dictionary mapping residue numbers from original sequence to cropped sequence.
    :vartype cropmap: dict [int, int]
    :ivar cropbackmap: A dictionary mapping residue numbers from cropped sequence to original sequence.
    :vartype cropbackmap: dict [int, int]
    :ivar msa: A free variable not used by CROPS itself.
    :vartype msa: Any
    :ivar cropmsa: A free variable not used by CROPS itself.
    :vartype cropmsa: Any
    :ivar intervals: The integer interval object containing the cropping information.
    :vartype intervals: :class:`crops.elements.intervals.intinterval`

    :raises `TypeError`: For wrong input formats.

    :example:

    >>> from crops.elements import sequences as ces
    >>> myseq = ces.sequence(seqid='1', oligomer = 'exampleID')
    >>> myseq.mainseq('GATTACA')
    >>> myseq.mainseq()
    'GATTACA'
    >>> myseq.chains = {'A', 'B'}
    >>> myseq.addseq('gapseq','GAT--C-')
    >>> myseq.addseq('cobra','TACATACA')
    >>> myseq.length()
    7
    >>> myseq.ngaps('gapseq')
    3
    >>> myseq.guess_biotype()
    'DNA'
    >>> print(myseq)
    Sequence object >EXAMPLEID_1|Chains A,B (seq=GATTACA, type=DNA, length=7)
    >>> myseq.source = 'Example'
    >>> myseq.addseq('cropseq', '+A+T++')
    >>> myseq.addseq('cropgapseq', '+A+-++')
    >>> myseq.full_length()
    7
    >>> myseq.mainseq('AT')
    'AT'
    >>> myseq.ncrops()
    4
    >>> myseq.update_cropsheader()
    >>> myseq.cropinfo()
    '#Residues cropped: 4 (1 not from terminals) ; % cropped: 66.67 (16.67 not from terminal segments)'
    >>> myseq.dump(out='string')
    '>crops|exampleID_1|Chains A,B|Source: Example|#Residues cropped: 4 (1 not from terminal segments) ; % cropped: 66.67 (16.67 not from terminal segments)\\nAT\\n'

    :example:

    >>> from crops.elements import sequences as ces
    >>> from crops.iomod import parsers as cip
    >>> myseq = cip.parseseqfile('7M6C.fasta')
    >>> myseq
    Sequence object: (>7M6C_1|Chain A, seq=MRTLWIMAVL[...]KPLCKKADPC, type=Undefined, length=138)
    >>> myseq.guess_biotype()
    'Protein'
    >>> myseq
    Sequence object: (>7M6C_1|Chain A, seq=MRTLWIMAVL[...]KPLCKKADPC, type=Protein, length=138)
    """

    _kind = 'Sequence'
    __slots__ = ['oligomer_id', 'name', 'chains', 'source', 'seqs', 'biotype',
                 'source_headers', 'crops_header', 'cropmap', 'cropbackmap',
                 'infostring', 'msa', 'cropmsa', 'intervals']

    def __init__(self, seqid=None, oligomer=None, seq=None, chains=None,
                 source=None, header=None, biotype=None, extrainfo=None):
        self.oligomer_id = None
        self.name = None
        self.chains = set()
        self.source = None
        self.source_headers = []
        self.crops_header = None
        self.seqs = {}
        self.biotype = None
        self.infostring = None
        self.cropmap = None
        self.cropbackmap = None
        self.msa = None
        self.cropmsa = None
        self.intervals = None

        if header is not None:
            if isinstance(header, str):
                self.source_headers.append(header)
                try:
                    header_info = retrieve_id(header)
                except Exception:
                    logging.warning('Header format not recognised. Information not extracted.')
                    header_info = None
            else:
                logging.critical("Argument 'header' should be a string.")
                raise TypeError
        else:
            header_info = None

        if seqid is not None:
            if isinstance(seqid, str):
                self.name = seqid
            elif isinstance(seqid, int):
                self.name = str(seqid)
            else:
                logging.critical("Sequence ID 'seqid' should be a string.")
                raise TypeError
        else:
            if header_info is not None:
                if 'seqid' in header_info:
                    self.name = header_info['seqid']
            else:
                self.name = '1'

        if seq is not None:
            if isinstance(seq, str):
                self.seqs['mainseq'] = seq
            else:
                logging.critical("Chain sequence 'seq' should be a string.")
                raise TypeError
        else:
            self.seqs['mainseq'] = ''

        if oligomer is not None:
            if isinstance(oligomer, str):
                self.oligomer_id = oligomer
            else:
                logging.critical("Oligomer ID 'oligomer' should be a string.")
                raise TypeError
        else:
            if header_info is not None:
                if 'mainid' in header_info:
                    self.oligomer_id = header_info['mainid']

        if chains is not None:
            if isinstance(chains, set):
                for ch in chains:
                    if isinstance(ch, str):
                        self.chains.add(ch)
                    else:
                        logging.critical("Chain IDs in 'chains' set should be strings.")
                        raise TypeError
            else:
                logging.critical("Argument 'chains' should be a set of strings.")
                raise TypeError
        else:
            if header_info is not None:
                if 'chains' in header_info:
                    self.chains = header_info['chains']

        if source is not None:
            if isinstance(source, str):
                self.source = source
            else:
                logging.critical("Argument 'source' should be a string.")
                raise TypeError
        else:
            if header_info is not None:
                if 'source' in header_info:
                    self.source = header_info['source']

        if biotype is not None:
            if biotype.lower() == 'guess':
                self.biotype = guess_type(seq)
            else:
                self.biotype = biotype
        else:
            self.biotype = None

        if extrainfo is not None:
            if isinstance(extrainfo, str):
                self.infostring = extrainfo
            else:
                logging.critical("Argument 'extrainfo' should be a string.")
                raise TypeError
        else:
            if header_info is not None:
                if 'comments' in header_info:
                    self.infostring = header_info['comments']
            else:
                self.infostring = ""

        if oligomer is None:
            self.crops_header = makeheader(mainid='NOID', seqid=self.name,
                                           chains=self.chains, source=self.source,
                                           extrainfo=self.infostring)
        else:
            self.crops_header = makeheader(mainid=self.oligomer_id, seqid=self.name,
                                           chains=self.chains, source=self.source,
                                           extrainfo=self.infostring)

    def __repr__(self):
        chtype = self.biotype if self.biotype is not None else 'Undefined'
        if 'mainseq' not in self.seqs:
            logging.critical("'mainseq' sequence not found.")
            raise ValueError
        if len(self.seqs['mainseq']) <= 20:
            showseq = self.seqs['mainseq']
        else:
            showseq = (self.seqs['mainseq'][:10]+'[...]' +
                       self.seqs['mainseq'][len(self.seqs['mainseq'])-10:])
        tempolig = self.oligomer_id if self.oligomer_id is not None else 'NOID'
        shortid = makeheader(mainid=tempolig, seqid=self.name,
                             chains=self.chains, short=True)
        string = (self._kind+" object "+shortid+" (seq="+str(showseq) +
                  ", type=" + chtype + ", length=" +
                  str(len(self.seqs['mainseq']))+")")
        return string

    def __iter__(self):
        return iter(self.seqs['mainseq'].values())

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    def addseq(self, newid, newseq):
        """Add sequence to `seqs` dictionary.

        :param newid: New sequence's identifier.
        :type newid: str
        :param newseq: New sequence.
        :type newseq: str

        :raises `TypeError`: If newid is not a string.
        :raises `KeyError`: If newseq is not a string.

        """
        if not isinstance(newid, str):
            logging.critical("New sequence ID 'newid' should be a string.")
            raise TypeError
        if not isinstance(newseq, str):
            logging.critical("New sequence string 'newseq' should be a string.")
            raise TypeError
        if newid in self.seqs:
            logging.critical("Key name 'newid' already exists.")
            raise KeyError

        self.seqs[newid] = newseq

    def delseq(self, delid=None, wipeall=False):
        """Delete sequence(s) from the `seqs` dictionary.

        :param delid: ID of sequence to be deleted, defaults to None.
        :type delid: str, optional
        :param wipeall: If True, all the sequences are deleted, defaults to False.
        :type wipeall: bool, optional

        :raises `TypeError`: If delid is not a string or wipeall is not a boolean.

        """
        if not isinstance(delid, str) and delid is not None:
            logging.critical("Sequence ID 'delid' should be a string.")
            raise TypeError
        if not isinstance(wipeall, bool):
            logging.critical("Boolean switch 'wipeall' is neither True nor False.")
            raise TypeError

        if wipeall:
            self.seqs = {}
            self.seqs['mainseq'] = ''
            return
        if delid is None:
            return

        if delid == 'mainseq':
            self.seqs['mainseq'] = ''
        else:
            self.seqs.pop(delid)

    def mainseq(self, add=None):
        """Return or modifies the main sequence.

        :param add: If given, the main sequence is replaced by 'add', defaults to None.
        :type add: str, optional

        :raises `TypeError`: If 'add' is given and is not a string.

        :return: The (new) main sequence.
        :rtype: str

        """
        if not isinstance(add, str) and add is not None:
            logging.critical("If included, sequence 'add' should be a string.")
            raise TypeError

        if add is not None:
            self.seqs['mainseq'] = add

        return self.seqs['mainseq']

    def guess_biotype(self):
        """Save the guessed biotype and return it.

        :return: Guessed biotype.
        :rtype: str

        """
        if self.seqs['mainseq'] is None:
            self.biotype = None
        else:
            self.biotype = guess_type(self.seqs['mainseq'])

        return self.biotype

    def dump(self, out, split=False, oneline=False):
        """Write header and main sequence to a file. If the file exists, output is appended.

        :param out: An output filepath (str), 'string', or an open file.
        :type out: str, file
        :param split: If True, identical sequences are dumped for every chain, defaults to False.
        :type split: bool, optional
        :param oneline: If True, sequences are not split in 80 residue-lines, defaults to False.
        :type oneline: bool, optional

        :raises `TypeError`: If `out` is neither a string nor an open file.
        :raises `KeyError`: If object contains no chains.

        :return: A string containing the output if and only if out=='string'.
        :rtype: str

        """
        if not isinstance(out, str) and not isinstance(out, io.IOBase):
            logging.critical("Argument 'out' should be a string or a file.")
            raise TypeError

        if (self.chains is None or
                (isinstance(self.chains, set) and len(self.chains) == 0)):
            logging.critical('No chains defined in sequence.')
            raise KeyError

        outheader = []

        if split:
            chset = []
            for ch in self.chains:
                chset.append({ch})
        else:
            chset = [self.chains]

        if self.oligomer_id is None:
            tag1 = 'NoID'
        else:
            tag1 = self.oligomer_id
        tag2 = self.infostring
        if self.ncrops() > 0:
            if self.infostring[-1] != "|":
                tag2 += '|'
            tag2 += self.cropinfo()

        for ch in chset:
            outheader.append(makeheader(mainid=tag1,
                                        seqid=self.name,
                                        chains=ch,
                                        source=self.source,
                                        extrainfo=tag2))

        if not oneline:
            lenseq = len(self.seqs['mainseq'])
            nlines = int((lenseq-1)/80)+1
        output = ''
        for header in outheader:
            if isinstance(out, io.IOBase) is True:
                out.write(header+'\n')
                if oneline:
                    out.write(self.seqs['mainseq']+'\n')
                else:
                    for n in range(nlines):
                        out.write(self.seqs['mainseq'][n*80:(n+1)*80]+'\n')
            else:
                output += header + os.linesep
                if oneline:
                    output += self.seqs['mainseq'] + os.linesep
                else:
                    for n in range(nlines):
                        output += self.seqs['mainseq'][n*80:(n+1)*80] + os.linesep
        if isinstance(out, io.IOBase) is False:
            if out.lower() == 'string':
                return output
            else:
                outpath = out
                op = 'a' if os.path.isfile(outpath) else 'w'
                with open(outpath, op) as out:
                    out.write(output)
        return

    def dumpmap(self, out, split=False):
        """Write header and cropmap to a file. If file exists, output is appended.

        :param out: An output filepath (str) or an open file.
        :type out: str, file
        :param backmap: If True, the output will be self.cropbackmap, defaults to False.
        :type backmap: bool, optional
        :param split: If True, identical maps are dumped for every chain, defaults to False.
        :type split: bool, optional

        :raises `TypeError`: If `out` is neither a string nor an open file.
        :raises `ValueError`: If one or both of `cropmap` and `cropbackmap` are empty.
        :raises `KeyError`: If object contains no chains.

        """
        if not isinstance(out, str) and not isinstance(out, io.IOBase):
            logging.critical("Argument 'out' should be a string or a file.")
            raise TypeError

        if self.cropmap is None:
            stringerr = "Cropmap not found in sequence."
            logging.critical(stringerr)
            raise ValueError

        if (self.chains is None or
                (isinstance(self.chains, set) and len(self.chains) == 0)):
            logging.critical('No chains defined in sequence.')
            raise KeyError

        outheader = []

        if split:
            chset = []
            for ch in self.chains:
                chset.append({ch})
        else:
            chset = [self.chains]

        if self.oligomer_id is None:
            tag1 = 'NoID'
        else:
            tag1 = self.oligomer_id
        tag2 = self.infostring
        if self.ncrops() > 0:
            if self.infostring[-1] != "|":
                tag2 += '|'
            tag2 += self.cropinfo()

        for ch in chset:
            outheader.append(makeheader(mainid=tag1,
                                        seqid=self.name,
                                        chains=ch,
                                        source=self.source,
                                        extrainfo=tag2))
        output = ''
        for header in outheader:
            if isinstance(out, io.IOBase):
                out.write(header+'\n')
                for key, value in self.cropmap.items():
                    if value is not None:
                        out.write(str(key)+'  '+str(value)+'\n')
                    else:
                        out.write(str(key)+'  0\n')
            else:
                output += header + os.linesep
                for key, value in self.cropmap.items():
                    if value is not None:
                        output += str(key) + '  ' + str(value) + os.linesep
                    else:
                        output += str(key) + '  0' + os.linesep

        if isinstance(out, io.IOBase) is False:
            if out.lower() == 'string':
                return output
            else:
                outpath = out
                op = 'a' if os.path.isfile(outpath) else 'w'
                with open(outpath, op) as out:
                    out.write(output)
        return

    def length(self):
        """Return the length of the main sequence.

        :return: Length of the main sequence.
        :rtype: int

        """
        return len(self.seqs['mainseq'])

    def full_length(self):
        """Return the length of the full sequence. If not found, the main sequence
        will be considered the full sequence, and will be saved as so.

        :return: Length of the full sequence.
        :rtype: int

        """
        if 'fullseq' not in self.seqs:
            self.seqs['fullseq'] = self.seqs['mainseq']

        return len(self.seqs['fullseq'])

    def ngaps(self, seqid='gapseq'):
        """Return the number of gaps ('-') in a sequence.

        :param seqid: The ID of the sequence containing the gaps, defaults to 'gapseq'.
        :type seqid: str, optional

        :raises `TypeError`: If seqid is not a string.

        :return: Number of gaps in `seqid`. If 'gapseq' is a list of several models, a list is returned. If `seqid` not found, 0 is returned.
        :rtype: int or list [int]

        """
        if not isinstance(seqid, str):
            logging.critical("Sequence ID 'seqid' should be a string.")
            raise TypeError
        if seqid in self.seqs:
            if isinstance(self.seqs[seqid], str):
                nseqid = [self.seqs[seqid]]
            else:
                nseqid = self.seqs[seqid]
            ng = []
            for altseq in nseqid:
                n = 0
                for char in altseq:
                    if char == '-':
                        n += 1
                ng.append(n)
            if len(ng) == 1:
                ng = ng[0]
        else:
            ng = 0

        return ng

    def ncrops(self, seqid='cropseq', offterminals=False, offmidseq=False):
        """Return the number of cropped elements ('+','*') in a sequence.

        :param seqid: The ID of the sequence containing the cropped elements, defaults to 'cropseq'.
        :type seqid: str, optional
        :param offterminals: Count those removed from terminal segments only, defaults to False.
        :type offterminals: bool, optional
        :param offmidseq: Count those removed NOT from terminal segments only, defaults to False.
        :type offmidseq: bool, optional

        :raises `TypeError`: If `seqid` is not a string, or `offterminals`, `offmidseq` are not boolean.

        :return: Number of cropped elements in `seqid` according to interval chosen. If `seqid` not found, 0 is returned.
        :rtype: int

        """
        if not isinstance(seqid, str):
            logging.critical("Sequence ID 'seqid' should be a string.")
            raise TypeError

        n = 0
        if seqid not in self.seqs:
            return n

        for char in self.seqs[seqid]:
            if char == '+' or char == '*':
                n += 1

        if ((offterminals is False and offmidseq is False) or
                (offterminals is True and offmidseq is True)):
            return n
        else:
            nterms = 0
            for char in self.seqs[seqid]:
                if char == '+' or char == '*':
                    nterms += 1
                else:
                    break
            for char in reversed(self.seqs[seqid]):
                if char == '+' or char == '*':
                    nterms += 1
                else:
                    break

        if offterminals is False and offmidseq is True:
            return n-nterms
        elif offterminals is True and offmidseq is False:
            return nterms

    def update_cropsheader(self):
        """Update `cropsheader`. Useful after updating any information from the sequence."""
        if self.oligomer_id is None:
            tag1 = 'NoID'
        else:
            tag1 = self.oligomer_id
        tag2 = self.infostring
        if self.ncrops() > 0:
            if tag2[-1] != "|":
                tag2 += '|'
            tag2 += self.cropinfo()

        self.crops_header = makeheader(mainid=tag1,
                                       seqid=self.name,
                                       chains=self.chains,
                                       source=self.source,
                                       extrainfo=tag2)

    def cropinfo(self):
        """Return a string containing statistics about the cropped residues.

        :return: Statistics on number of crops.
        :rtype: str

        """
        cropstr = ""
        if 'cropseq' in self.seqs:
            cropstr += ('#Residues cropped: ')
            if self.ncrops() == 0:
                cropstr += '0'
            else:
                cropstr += (str(self.ncrops()) + ' (' +
                            str(self.ncrops(offmidseq=True)) +
                            ' not from terminal segments) ' +
                            '; % cropped: ' +
                            str(round(100*self.ncrops()/len(self.seqs['cropseq']), 2)) +
                            ' (' + str(round(100*self.ncrops(offmidseq=True)/len(self.seqs['cropseq']), 2)) +
                            ' not from terminal segments)')
        else:
            pass

        return cropstr


class oligoseq:
    """An object grouping several :class:`crops.elements.sequences.sequence` objects pertaining to a common oligomer.

    :param oligomer_id: Oligomer identifier (e.g. PDB id), defaults to None.
    :type oligomer_id: str
    :param imer: Container of several :class:`crops.elements.sequences.sequence` objects making up the oligomer, defaults to empty dict.
    :type imer: dict [str, :class:`crops.elements.sequences.sequence`], optional

    :ivar id: Oligomer sequence identifier (e.g. PDB id).
    :vartype id: str
    :ivar imer: Container of several :class:`crops.elements.sequence.monomer_sequence` making up the oligomer.
    :vartype imer: dict [str, :class:`crops.elements.sequence.monomer_sequence`]

    :raises `TypeError`: If the input formats are wrong.

    :example:

    >>> from crops.elements import sequences as ces
    >>> my_oligoseq = ces.oligoseq(oligomer_id='exampleID')
    >>> my_oligoseq.add_monomer
    >>> my_sequence.add_monomer('header_example','GATTACA',nid='mychain')
    >>> my_sequence.add_monomer('another_header','TACATACA')
    >>> my_sequence.nchains()
    2
    >>> my_sequence.length('mychain')
    7
    >>> my_sequence.write('/path/to/output/dir/')
    >>> print(my_sequence)
    docs Protein/polynucleotide sequence object: (id='example_id', # chains = 2)
    >>> my_sequence.purge()
    >>> my_sequence.nchains()
    0

    """
    _kind = 'Multiple sequence'
    __slots__ = ['id', 'imer']

    def __init__(self, oligomer_id=None, imer=None):
        if not isinstance(oligomer_id, str) and oligomer_id is not None:
            logging.critical("'oligomer_id' should be a string.")
            raise TypeError
        if not isinstance(imer, dict) and imer is not None:
            logging.critical("Sequence container 'imer' should be a dictionary.")
            raise TypeError
        elif isinstance(imer, dict):
            for val in imer.values():
                if not isinstance(val, sequence):
                    logging.critical("Sequence container 'imer' should only "
                                    "contain :class:`~crops.elements.sequences.sequence` objects.")
                    raise TypeError
        self.id = oligomer_id
        self.imer = imer if imer is not None else {}

    def __repr__(self):
        string = self._kind+" object: (id="+ str(self.id) + ", sequences = "+str(self.imer)+")"
        return string

    def __getitem__(self, key):
        return self.imer[key]

    def __len__(self):
        return len(self.imer)

    def __iter__(self):
        return iter(self.imer.values())

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    def purge(self):
        """Clear the object's content without deleting the object itself."""
        self.id = None
        self.imer.clear()

    def add_sequence(self, newseq):
        """Add a new :class:`crops.elements.sequences.sequence` to the object.

        :param newseq: Sequence object.
        :type newseq: :class:`crops.elements.sequences.sequence`

        :raises `TypeError`: If `newseq` is not a :class:`crops.elements.sequences.sequence` object.
        :raises `Exception`: If sequence content is incompatible with that in oligoseq (oligomer id, other sequences, etc).

        """
        addall = None
        errormsg = ('Sequence content is incompatible with oligoseq ' +
                    self.id + '.')
        if self.id is not None and newseq.oligomer_id is not None:
            if self.id.upper() != newseq.oligomer_id.upper():
                logging.critical(errormsg)
                raise ValueError

        if newseq.name is not None:
            if newseq.name in self.imer:
                if self.imer.seqs['mainseq'] == newseq.seqs['mainseq']:
                    addall = False
                else:
                    logging.critical(errormsg)
                    raise ValueError
            else:
                for seq in self.imer.values():
                    if seq.seqs['mainseq'] == newseq.seqs['mainseq']:
                        logging.critical(errormsg)
                        raise ValueError
                addall = True
        else:
            for seq in self.imer.values():
                if seq.seqs['mainseq'] == newseq.seqs['mainseq']:
                    addall = False
                    newseq.name = seq.name
                    break
            if addall is not False:
                addall = True

        if addall is True:
            for ch in newseq.chains:
                for seq in self.imer.values():
                    if ch in seq.chains:
                        logging.critical(errormsg)
                        raise ValueError
            if newseq.name is None:
                n = 1
                while True:
                    if str(n) in self.imer:
                        n += 1
                    else:
                        newseq.name = str(n)
                        break
            self.imer[newseq.name] = newseq
            if self.id is None and newseq.oligomer_id is not None:
                self.id = newseq.oligomer_id.upper()
                for seq in self.imer.values():
                    seq.oligomer_id = newseq.oligomer_id.upper()
            elif self.id is not None and newseq.oligomer_id is None:
                self.imer[newseq.name].oligomer_id = self.id.upper()
        else:
            for ch in newseq.chains:
                for seq in self.imer.values():
                    if ch in seq.chains and seq.name != newseq.name:
                        logging.critical(errormsg)
                        raise ValueError
                self.imer[newseq.name].chains.add(ch)

            for header in newseq.source_headers:
                if header not in self.imer[newseq.name].source_headers:
                    self.imer[newseq.name].source_headers.append(header)
            if newseq.source != self.imer[newseq.name].source:
                self.imer[newseq.name].source = 'Diverse'

        self.imer[newseq.name].update_cropsheader()

        return

    def del_sequence(self, seqid):
        """Remove the selected :class:`crops.elements.sequences.sequence` from the object.

        :param seqid: Doomed sequence's identifier.
        :type seqid: str

        :raises `TypeError`: If `seqid` is not a string.

        """
        if isinstance(seqid, int):
            seqid = str(seqid)
        if not isinstance(seqid, str):
            raise TypeError("'seqid' should be a string.")

        if seqid in self.imer:
            self.imer.pop(seqid)
        else:
            logging.warning('Sequence named ' + seqid + ' not found in oligoseq.')

        return

    def set_cropmaps(self, mapdict, cropmain=False):
        """Sets the parsed cropmaps from :class:`crops.iomod.parsers.parsemapfile`.

        :param mapdict: Parsed maps for this specific object.
        :type mapdict: dict [str, dict [str, dict [int, int]]]
        :param cropmain: If True, it will crop 'mainseq' and generate 'fullseq' and 'cropseq'. If 'mainseq' has been edited before this operation will yield wrong results, defaults to False.
        :type cropmain: bool, optional

        :raises `TypeError`: When `mapdict` has not the appropriate format.

        """
        if not isinstance(mapdict, dict):
            logging.critical("'mapdict' should be a dictionary.")
            raise TypeError

        for seqid in mapdict:
            if not isinstance(seqid, str):
                logging.critical("Values in 'mapdict' should be strings.")
                raise TypeError
            if seqid in self.imer:
                if ('cropmap' not in mapdict[seqid] or
                        'cropbackmap' not in mapdict[seqid]):
                    logging.critical("'mapdict' is not a crop map.")
                    raise TypeError
                self.imer[seqid].cropmap = copy.deepcopy(mapdict[seqid]['cropmap'])
                self.imer[seqid].cropbackmap = copy.deepcopy(mapdict[seqid]['cropbackmap'])
                self.imer[seqid].intervals = intinterval(description=self.id+'_'+str(seqid))
                for resc, res0 in mapdict[seqid]['cropbackmap'].items():
                    self.imer[seqid].intervals = self.imer[seqid].intervals.union(other=res0)
                if cropmain is True:
                    self.imer[seqid].seqs['fullseq'] = self.imer[seqid].seqs['mainseq']
                    self.imer[seqid].seqs['mainseq'] = ''
                    self.imer[seqid].seqs['cropseq'] = ''
                    for n in range(len(self.imer[seqid].seqs['fullseq'])):
                        if self.imer[seqid].cropmap[n+1] is None:
                            self.imer[seqid].seqs['cropseq'] += '+'
                        else:
                            self.imer[seqid].seqs['mainseq'] += self.imer[seqid].seqs['fullseq'][n]
                            self.imer[seqid].seqs['cropseq'] += self.imer[seqid].seqs['fullseq'][n]
                self.imer[seqid].infostring += '|' + self.imer[seqid].cropinfo()
                self.imer[seqid].update_cropsheader()

        return

    def write(self, outdir, infix="", split=False, oneline=False):
        """Write all :class:`crops.elements.sequences.sequence` objects to .fasta file or string.

        :param outdir: Output directory.
        :type outdir: str
        :param infix: Filename tag to distinguish from original input file, defaults to "".
        :type infix: str, optional
        :param split: If True, identical sequences are dumped for each chain, defaults to False.
        :type split: bool, optional
        :param oneline: If True, sequences are not split in 80 residue-lines, defaults to False.
        :type oneline: bool, optional

        :raises `FileNotFoundError`: Output directory not found.

        """
        if not os.path.isdir(outdir) and outdir != 'string':
            logging.critical(outdir + ' directory not found.')
            raise FileNotFoundError

        if outdir == 'string':
            outpath = 'string'
            outstring = ""
            for seq in self.imer.values():
                outstring += seq.dump(outpath, split=split, oneline=oneline)

            return outstring
        else:
            outpath = os.path.join(outdir, self.seq_id + infix + ".fasta")
            for seq in self.imer.values():
                seq.dump(outpath, split=split, oneline=oneline)

            return

    def length(self, seqid):
        """Return the length of a certain sequence.

        :param seqid: ID of :class:`crops.elements.sequences.sequence`.
        :type seqid: str

        :raises `TypeError`: When 'seqid' is not a string.
        :raises `KeyError`: Specific sequence not found in :class:`crops.elements.sequences.oligoseq`.

        :return: Length of :class:`crops.elements.sequences.sequence`.
        :rtype: int

        """
        if isinstance(seqid, int):
            seqid = str(seqid)
        if not isinstance(seqid, str):
            logging.critical('chain input must be a string.')
            raise TypeError
        if seqid in self.imer:
            return self.imer[seqid].length()
        else:
            logging.critical(seqid+' monomer not found in sequence.')
            raise KeyError

    def nchains(self):
        """Return number of chains in object, counting all sequence objects contained.

        :return: Number of chains in object, counting al :class:`crops.elements.sequences.sequence` contained.
        :rtype: int
        """
        n = 0
        for seqid in self.imer:
            n += len(self.imer[seqid].chains)

        return n

    def nseqs(self):
        """Return number of sequence objects in object.

        :return: Number of :class:`crops.elements.sequences.sequence` objects in object.
        :rtype: int
        """
        return len(self.imer)

    def chainlist(self):
        """Return a set with all the chain names in the object.

        :return: Chain names in :class:`crops.elements.sequences.oligoseq`.
        :rtype: set [str]

        """
        newset = set()
        for seqid in self.imer:
            newset = newset.union(self.imer[seqid].chains)

        return newset

    def whatseq(self, chain):
        """Return the sequence number corresponding to a given chain.

        :param chain: The chain ID.
        :type chain: str

        :return: The :class:`crops.elements.sequences.sequence` of that chain.
        :rtype: str

        """
        myseq = None
        for seqid in self.imer:
            if chain in self.imer[seqid].chains:
                myseq = seqid
                break

        return myseq
