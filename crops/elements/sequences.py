from crops.about import __prog__, __description__, __author__, __date__, __version__

import os
import io
import copy
import logging

from crops.io.taggers import retrieve_id
from crops.io.taggers import makeheader
from crops.libs.rescodes import reslist
from crops.libs.rescodes import nuclist

def guess_type(inseq):
    """Returns the biological type of the sequence as guessed from residue types.

    :param inseq: Sequence to be evaluated.
    :type inseq: str
    :return: Sequence type ('Protein' or 'DNA' or 'RNA' or 'Unknown').
    :rtype: str

    """
    if not isinstance(inseq, str):
        raise TypeError("Sequence 'inseq' should be a string.")

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
    """A :class:`~crops.elements.sequences.sequence` object representing a single chain sequence.
    The :class:`~crops.elements.sequences.sequence` class represents a data structure to hold all
    sequence versions and other useful information characterising it.
    It contains functions to store, manipulate and organise sequence versions.

    :param seqid: Sequence identifier. Can be used alone or together with oligomer ID.
    :type seqid: str
    :param oligomer: Oligomer identifier. Sometimes as important as seqid.
    :type oligomer: str, optional
    :param seq: Sequence string.
    :type seq: str, optional

    :param header: Standard .fasta header, starting with ">".
    :type header: str, optional
    :ivar info: Useful information of the :class:`~crops.elements.sequence.monomer_sequence`.
    :vartype info: dict [str, any]
    :ivar seqs: The set of sequences, including default "mainseq", in :class:`~crops.elements.sequence.monomer_sequence`.
    :vartype seqs: dict [str, str]

    :example:

    >>> from crops.elements import sequences as csq
    >>> myseq = csq.sequence('exampleID')
    >>> mysq.mainseq('GATTACA')
    >>> myseq.mainseq()
    'GATTACA'
    >>> myseq.addseq('gapseq','GAT--C-')
    >>> myseq.addseq('cobra','TACATACA')
    >>> myseq.length()
    7
    >>> myseq.ngaps('gapseq')
    3
    >>> myseq.guess_biotype()
    'DNA'
    >>> print(myseq)
    Sequence object: ('>exampleID|Chain exampleID', seq='GATTACA', type='DNA', length=7)

    :example:

    >>> from crops.elements import sequences as csq
    >>> from crops.io import parsers as csp
    >>> myseq = csp.parseseqfile('7M6C.fasta')
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
                 'infostring']
    def __init__(self, seqid=None, oligomer=None, seq=None, chains=None, source=None,
                 header=None, biotype=None, extrainfo=None):
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

        if seqid is not None:
            if isinstance(seqid, str):
                self.name = seqid
            else:
                raise TypeError("Sequence ID 'seqid' should be a string.")

        if seq is not None:
            if isinstance(seq, str):
                self.seqs['mainseq'] = seq
            else:
                raise TypeError("Chain sequence 'seq' should be a string.")
        else:
            self.seqs['mainseq'] = ''

        if oligomer is not None:
            if isinstance(oligomer, str):
                self.oligomer_id = oligomer
            else:
                raise TypeError("Oligomer ID 'oligomer' should be a string.")

        if chains is not None:
            if isinstance(chains, set):
                for ch in chains:
                    if isinstance(ch, str):
                        self.chains.add(ch)
                    else:
                        raise TypeError("Chain IDs in 'chains' set should be strings.")
            else:
                raise TypeError("Argument 'chains' should be a set of strings.")

        if source is not None:
            if isinstance(source, str):
                self.source = source
            else:
                raise TypeError("Argument 'source' should be a string.")

        if header is not None:
            if isinstance(header, str):
                self.source_headers.append(header)
            else:
                raise TypeError("Argument 'source' should be a string.")

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
                raise TypeError("Argument 'extrainfo' should be a string.")

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
            raise ValueError("'mainseq' sequence not found.")
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

    def addseq(self,newid,newseq):
        """Add sequence to :class:`~crops.elements.sequences.sequence`.

        :param newid: New sequence's identifier.
        :type newid: str
        :param newseq: New sequence.
        :type newseq: str
        :raises TypeError: If newid is not a string.
        :raises KeyError: If sequence is not a string.

        """
        if not isinstance(newid, str):
            raise TypeError("New sequence ID 'newid' should be a string.")
        if not isinstance(newseq, str):
            raise TypeError("New sequence string 'newseq' should be a string.")
        if newid in self.seqs:
            raise KeyError("Key name 'newid' already exists.")

        self.seqs[newid] = newseq

    def delseq(self, delid=None, wipeall=False):
        """Deletes sequence(s) from :class:`~crops.elements.sequences.sequence`.

        :param delid: ID of sequence to be deleted, defaults to None.
        :type delid: str, optional
        :param wipeall: If True, all the sequences are deleted, defaults to False.
        :type wipeall: bool, optional
        :raises TypeError: If delid is not a string.

        """
        if not isinstance(delid, str):
            raise TypeError("Sequence ID 'delid' should be a string.")
        if not isinstance(wipeall, bool):
            raise TypeError("Boolean switch 'wipeall' is neither True nor False.")

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
        """Returns or modifies the main sequence.

        :param add: If given main sequence is changed to 'add' sequence, defaults to None.
        :type add: str, optional
        :raises TypeError: If 'add' is given and is not a string.
        :return: If 'add' is None, the main sequence is returned.
        :rtype: str

        """
        if not isinstance(add, str) and add is not None:
            raise TypeError("If included, sequence 'add' should be a string.")

        if add is not None:
            self.seqs['mainseq'] = add

        return self.seqs['mainseq']

    def guess_biotype(self):
        """Saves the guessed biotype and returns it.

        :return: Guessed biotype.
        :rtype: str

        """
        if self.seqs['mainseq'] is None:
            self.biotype = None
        else:
            self.biotype = guess_type(self.seqs['mainseq'])

        return self.biotype

    def dump(self, out, split=False, oneline=False):
        """Writes header and main sequence to a file. If file exists, output is appended.

        :param out: An output filepath (str) or an open file.
        :type out: str, file
        :param split: If True, identical sequences are dumped for every chain, defaults to False.
        :type split: bool, optional
        :param oneline: If True, sequences are not split in 80 residue-lines, defaults to False.
        :type oneline: bool, optional
        :raises TypeError: If out is neither a string nor an open file.

        """
        if not isinstance(out, str) and not isinstance(out, io.IOBase):
            raise TypeError("Argument 'out' should be a string or a file.")

        outheader = []

        if split:
            if (self.chains is None or
                    (isinstance(self.chains, set) and len(self.chains) == 0)):
                raise KeyError('No chains defined in sequence.')

            for ch in self.chains:
                outheader.append(makeheader(mainid=self.oligomer_id,
                                            seqid=self.name,
                                            chains={ch},
                                            source=self.source,
                                            extrainfo=self.infostring))
        else:
            outheader.append(self.crops_header)

        if not oneline:
            lenseq = len(self.seqs['mainseq'])
            nlines = int((lenseq-1)/80)+1
        for header in outheader:
            if isinstance(out, io.IOBase):
                out.write(header+'\n')
                if oneline:
                    out.write(self.seqs['mainseq']+'\n')
                else:
                    for n in range(nlines):
                        out.write(self.seqs['mainseq'][n*80:(n+1)*80]+'\n')
            else:
                outpath = out
                op = 'a' if os.path.isfile(outpath) else 'w'
                with open(outpath, op) as out:
                    out.write(header+'\n')
                    if oneline:
                        out.write(self.seqs['mainseq']+'\n')
                    else:
                        for n in range(nlines):
                            out.write(self.seqs['mainseq'][n*80:(n+1)*80]+'\n')

        return

    def dumpmap(self, out, split=False):
        """Writes header and cropmap to a file. If file exists, output is appended.

        :param out: An output filepath (str) or an open file.
        :type out: str, file
        :param backmap: If True, the output will be self.cropbackmap, defaults to False.
        :type backmap: bool, optional
        :param split: If True, identical sequences are dumped for every chain, defaults to False.
        :type split: bool, optional
        :raises TypeError: If out is neither a string nor an open file.
        :raises ValueError: If self.cropmap / self.cropbackmap does not exist.

        """
        if not isinstance(out, str) and not isinstance(out, io.IOBase):
            raise TypeError("Argument 'out' should be a string or a file.")

        if self.cropmap is None:
            stringerr = "Cropmap not found in sequence."
            raise ValueError(stringerr)

        outheader = []

        if split:
            if (self.chains is None or
                    (isinstance(self.chains, set) and len(self.chains) == 0)):
                raise KeyError('No chains defined in sequence.')

            for ch in self.chains:
                outheader.append(makeheader(mainid=self.oligomer_id,
                                            seqid=self.name,
                                            chains={ch},
                                            source=self.source,
                                            extrainfo=self.infostring))
        else:
            outheader.append(self.crops_header)

        for header in outheader:
            if isinstance(out, io.IOBase):
                out.write(header+'\n')
                for key, value in self.cropmap.items():
                    if value is not None:
                        out.write(str(key)+'  '+str(value)+'\n')
                    else:
                        out.write(str(key)+'  0\n')
            else:
                outpath = out
                op = 'a' if os.path.isfile(outpath) else 'w'
                with open(outpath, op) as out:
                    out.write(header+'\n')
                    for key, value in self.cropmap.items():
                        if value is not None:
                            out.write(str(key)+'  '+str(value)+'\n')
                        else:
                            out.write(str(key)+'  0\n')

        return

    def length(self):
        """Returns the length of the main sequence.

        :return: Length of the main sequence.
        :rtype: int

        """

        return len(self.seqs['mainseq'])

    def full_length(self):
        """Returns the length of the full sequence. If not found, the length of the main sequence is returned.

        :return: Length of the full sequence.
        :rtype: int

        """
        if 'fullseq' not in self.seqs:
            self.seqs['fullseq'] = self.seqs['mainseq']

        return len(self.seqs['fullseq'])

    def ngaps(self,seqid='gapseq'):
        """Returns the number of gaps ('-') in a sequence.

        :param seqid: The ID of the sequence containing the gaps, defaults to 'gapseq'.
        :type seqid: str, optional
        :raises TypeError: If seqid is not a string.
        :return: Number of gaps in seqid. If seqid not found, ngaps=0.
        :rtype: int

        """
        if not isinstance(seqid, str):
            raise TypeError("Sequence ID 'seqid' should be a string.")
        n = 0
        if seqid in self.seqs:
            for char in self.seqs[seqid]:
                if char == '-':
                    n += 1
        return n

    def ncrops(self, seqid='cropseq', offterminals=False, offmidseq=False):

        """Returns the number of cropped elements ('+','*') in a sequence.

        :param seqid: The ID of the sequence containing the cropped elements, defaults to 'cropseq'.
        :type seqid: str, optional
        :param offterminals: Count those removed from terminals only, defaults to False.
        :type offterminals: bool, optional
        :param offmidseq: Count those removed NOT from terminals only, defaults to False.
        :type offmidseq: bool, optional
        :raises TypeError: If seqid is not a string.
        :return: Number of cropped elements in seqid according to interval chosen. If seqid not found, ncrops=0.
        :rtype: int

        """
        if not isinstance(seqid, str):
            raise TypeError("Sequence ID 'seqid' should be a string.")

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
        if self.oligomer_id is None:
            self.crops_header = makeheader(mainid='NOID',
                                           seqid=self.name,
                                           chains=self.chains,
                                           source=self.source,
                                           extrainfo=self.infostring)
        else:
            self.crops_header = makeheader(mainid=self.oligomer_id,
                                           seqid=self.name,
                                           chains=self.chains,
                                           source=self.source,
                                           extrainfo=self.infostring)

class oligoseq:
    """A :class:`~crops.elements.sequences.oligoseq` object grouping several
    :class:`~crops.elements.sequences.sequence` objects pertaining to a
    common oligomer.

    :param oligomer_id: Oligomer identifier (e.g. PDB id), defaults to None.
    :type oligomer_id: str
    :param imer: Container of several :class:`~crops.elements.sequences.sequence` objects making up the oligomer, defaults to empty dict.
    :type imer: dict [str, :class:`~crops.elements.sequences.sequence`], optional
    :ivar seq_id: Sequence identifier (e.g. PDB id).
    :vartype seq_id: str
    :ivar imer: Container of several :class:`~crops.elements.sequence.monomer_sequence` making up the oligomer.
    :vartype imer: dict [str, :class:`~crops.elements.sequence.monomer_sequence`]
    :ivar source: Information concerning the source of the :class:`~crops.elements.sequence.Sequence` (e.g. Uniprot).
    :vartype source: str

    :example:

    >>> from crops.elements import Sequence
    >>> my_sequence = Sequence(seq_id='example_id', source='docs')
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
            raise TypeError("'oligomer_id' should be a string.")
        if not isinstance(imer, dict) and imer is not None:
            raise TypeError("Sequence container 'imer' should be a dictionary.")
        elif isinstance(imer, dict):
            for val in imer.values():
                if not isinstance(val, sequence):
                    raise TypeError("Sequence container 'imer' should only " +
                                    "contain :class:`~crops.elements.sequences.sequence` objects.")
        self.id = oligomer_id
        self.imer = imer if imer is not None else {}

    def __repr__(self):
        string = self._kind+" object: (id="+ str(self.id) + ", sequences = "+str(self.imer)+")"
        return string

    def __iter__(self):
        return iter(self.imer.values())

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    def purge(self):
        """Clears the :class:`~crops.elements.sequences.oligoseq` without deleting the object itself.

        """
        self.id = None
        self.imer.clear()

    def add_sequence(self, newseq):
        """Adds a new :class:`~crops.elements.sequences.sequence` to the :class:`~crops.elements.sequences.oligoseq`.

        :param seq: Sequence object.
        :type seq: :class:`~crops.elements.sequences.sequence`
        :raises TypeError: If 'seq' is not a :class:`~crops.elements.sequences.sequence` object.
        :raises Exception: If sequence content is incompatible with that in oligoseq.
        """
        addall = None
        errormsg = ('Sequence content is incompatible with oligoseq ' +
                    self.id + '.')
        if (newseq.oligomer_id is not None and self.id is not None and
                newseq.oligomer_id.lower() != self.id):
            raise Exception(errormsg)

        if newseq.name is not None:
            if newseq.name in self.imer:
                if self.imer.seqs['mainseq'] == newseq.seqs['mainseq']:
                    addall = False
                else:
                    raise Exception(errormsg)
            else:
                for seq in self.imer:
                    if seq.seqs['mainseq'] == newseq.seqs['mainseq']:
                        raise Exception(errormsg)
                addall = True
        else:
            for seq in self.imer:
                if seq.seqs['mainseq'] == newseq.seqs['mainseq']:
                    addall = False
                    newseq.name = seq.name
                    break
            if addall is not False:
                addall = True

        if self.id is not None and newseq.oligomer_id is not None:
            if self.id != newseq.oligomer_id.lower():
                raise Exception(errormsg)

        if addall is True:
            for ch in newseq.chains:
                for seq in self.imer:
                    if ch in seq.chains:
                        raise Exception(errormsg)
            if newseq.name is None:
                while True:
                    n = 1
                    if str(n) in self.imer:
                        n += 1
                    else:
                        newseq.name = str(n)
                        break
            self.imer[newseq.name] = newseq
            if self.id is None and newseq.oligomer_id is not None:
                self.id = newseq.oligomer_id.lower()
                for seq in self.imer:
                    seq.oligomer_id = newseq.oligomer_id.lower()
            elif self.id is not None and newseq.oligomer_id is None:
                self.imer[newseq.name].oligomer_id = self.id
        else:
            for ch in newseq.chains:
                for seq in self.imer:
                    if ch in seq.chains and seq.name != newseq.name:
                        raise Exception(errormsg)
                self.imer[newseq.name].chains.add(ch)

            for header in newseq.source_headers:
                if header not in self.imer[newseq.name].source_headers:
                    self.imer[newseq.name].source_headers.append(header)
            if newseq.source != self.imer[newseq.name].source:
                self.imer[newseq.name].source = 'Diverse'

        self.imer[newseq.name].update_cropsheader()

        return

    def del_sequence(self, seqid):
        """Removes the selected :class:`~crops.elements.sequences.sequence` from the :class:`~crops.elements.sequences.oligoseq`.

        :param seqid: Doomed sequence's identifier.
        :type seqid: str
        :raises TypeError: When 'seqid' is not a string.

        """
        if not isinstance(seqid, str):
            raise TypeError("'seqid' should be a string.")

        if seqid in self.imer:
            self.imer.pop(seqid)
        else:
            logging.warning('Sequence named ' + seqid + ' not found in oligoseq.')

        return

    def set_cropmaps(self, mapdict):
        """Sets the parsed cropmaps from :class:`~crops.io.parsers.parsemapfile`.

        :param mapdict: Parsed maps for this specific :class:`~crops.elements.sequences.oligoseq`.
        :type mapdict: dict [str, dict [str, dict [int, int]]]
        :raises TypeError: When 'mapdict' has not the appropriate format.

        """
        if not isinstance(mapdict, dict):
            raise TypeError("'mapdict' should be a dictionary.")

        for seqid in mapdict:
            if not isinstance(seqid, str):
                raise TypeError("Values in 'mapdict' should be strings.")
            if seqid in self.imer:
                if ('cropmap' not in mapdict[seqid] or
                        'cropbackmap' not in mapdict[seqid]):
                    raise TypeError("'mapdict' is not a crop map.")
                self.imer[seqid].cropmap = copy.deepcopy(mapdict[seqid]['cropmap'])
                self.imer[seqid].cropackmap = copy.deepcopy(mapdict[seqid]['cropbackmap'])
        return

    def write(self, outdir, infix="", split=False, oneline=False):
        """Writes all :class:`~crops.elements.sequences.sequence` objects to .fasta file.

        :param outdir: Output directory.
        :type outdir: str
        :param infix: Filename tag to distinguish from original input file, defaults to "".
        :type infix: str, optional
        :param split: If True, identical sequences are dumped for each chain, defaults to False.
        :type split: bool, optional
        :param oneline: If True, sequences are not split in 80 residue-lines, defaults to False.
        :type oneline: bool, optional
        :raises FileNotFoundError: Output directory not found.
        """

        if not os.path.isdir(outdir):
            raise FileNotFoundError(outdir + ' directory not found.')

        outpath = os.path.join(outdir, self.seq_id + infix + ".fasta")
        for seq in self.imer:
            seq.dump(outpath, split=split, oneline=oneline)

        return

    def length(self, seqid):
        """Returns the length of a certain sequence.

        :param seqid: ID of :class:`~crops.elements.sequences.sequence`.
        :type seqid: str
        :raises TypeError: When 'seqid' is not a string.
        :raises KeyError: Specific sequence not found in :class:`~crops.elements.sequences.oligoseq`.
        :return: Length of :class:`~crops.elements.sequences.sequence`.
        :rtype: int

        """
        if not isinstance(seqid, str):
            raise TypeError('chain input must be a string.')
        if seqid in self.imer:
            return self.imer[seqid].length()
        else:
            raise KeyError(seqid+' monomer not found in sequence.')

    def nchains(self):
        """Returns number of chains in :class:`~crops.elements.sequences.oligoseq`.

        :return: Number of chains in all :class:`~crops.elements.sequences.sequence` of :class:`~crops.elements.sequences.oligoseq`.
        :rtype: int
        """
        n = 0
        for seqid in self.imer:
            n += len(self.imer[seqid].chains)

        return n

    def nseqs(self):
        """Returns number of :class:`~crops.elements.sequences.sequence` objects in :class:`~crops.elements.sequences.oligseq`.

        :return: Number of :class:`~crops.elements.sequences.sequence` objects in :class:`~crops.elements.sequences.oligoseq`.
        :rtype: int
        """
        return len(self.imer)

    def chainlist(self):
        """
        Returns a set with all the chain names in the :class:`~crops.elements.sequences.oligseq`.

        :return: Chain names in :class:`~crops.elements.sequences.oligseq`.
        :rtype: set [str]

        """
        newset = set()
        for seqid in self.imer:
            newset = newset.union(self.imer[seqid].chains)

        return newset

    def whatseq(self, chain):
        """Returns the sequence number corresponding to a chain.

        :param chain: The chain ID.
        :type chain: str
        :return: The :class:`~crops.elements.sequences.sequence` of that chain.
        :rtype: str

        """
        for seqid in self.imer:
            if chain in self.imer[seqid].chains:
                myseq = seqid
                break

        return myseq
