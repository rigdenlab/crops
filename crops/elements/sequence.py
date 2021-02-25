from crops.about import __prog__, __description__, __author__, __date__, __version__

import os
import io
import copy
import logging

from crops.io.taggers import retrieve_id

class monomer_sequence:
    """A sequence object representing a single chain sequence.
    The monomer sequence class represents a data structure to hold all
    sequences and other useful information characterising the monomer.
    It contains functions to store, manipulate and organise sequences.

    :param chid: Chain ID.
    :type chid: str
    :param seq: Sequence string.
    :type seq: str, optional
    :param header: Standard .fasta header, starting with ">".
    :type header: str, optional
    :ivar info: Useful information of the :class:`~crops.elements.sequence.monomer_sequence`.
    :vartype info: dict [str, any]
    :ivar seqs: The set of sequences, including default "mainseq", in :class:`~crops.elements.sequence.monomer_sequence`.
    :vartype seqs: dict [str, str]

    :example:

    >>> from crops.elements import monomer_sequence
    >>> my_monomer = monomer_sequence('example_id')
    >>> my_monomer.mainseq('GATTACA')
    >>> my_monomer.mainseq()
    'GATTACA'
    >>> my_monomer.addseq('gapseq','GAT--C-')
    >>> my_monomer.addseq('cobra','TACATACA')
    >>> my_monomer.length()
    7
    >>> my_monomer.ngaps('gapseq')
    3
    >>> print(my_monomer)
    Single chain sequence object: (id='example_id', seq='GATTACA', length=7)

    """

    _kind='Single chain sequence'
    __slots__=['info','seqs']
    def __init__(self,chid,seq=None,header=None):
        self.info={}
        self.seqs={}
        if chid is not None:
            if isinstance(chid,str):
                self.info['chain_id']=chid
            else:
                raise TypeError("Chain ID 'chid' should be a string.")
        if seq is not None:
            if isinstance(seq,str):
                self.seqs['mainseq']=seq
            else:
                raise TypeError("Chain sequence 'seq' should be a string.")
        else:
            self.seqs['mainseq']=''
        if header is not None:
            if isinstance(header,str):
                self.info['header']=header
                self.info['oligomer_id']=retrieve_id(header)[0].lower()
            else:
                raise TypeError("Sequence chain header 'header' should be a string.")
        else:
            self.info['header']=None
            self.info['oligomer_id']=None

    def __repr__(self):
        if 'mainseq' not in self.seqs:
            raise ValueError('mainseq sequence not found.')
        showseq=self.seqs['mainseq'] if len(self.seqs['mainseq'])<=20 else self.seqs['mainseq'][:10]+'[...]'+self.seqs['mainseq'][len(self.seqs['mainseq'])-10:]
        if self.info['oligomer_id'] is not None:
            string=self._kind+" object: ( id="+ self.info['oligomer_id']+'_'+self.info['chain_id'] +", seq="+str(showseq)+", length="+str(len(self.seqs['mainseq']))+" )"
        else:
            string=self._kind+" object: ( id="+ self.info['chain_id'] +", seq="+str(showseq)+", length="+str(len(self.seqs['mainseq']))+" )"
        return string

    def __iter__(self):
        return iter(self.seqs.values())

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    def addseq(self,newid,newseq):
        """Add sequence to :class:`~crops.elements.sequence.monomer_sequence`.

        :param newid: New sequence's identifier.
        :type newid: str
        :param newseq: New sequence.
        :type newseq: str
        :raises TypeError: If newid is not a string.
        :raises KeyError: If sequence is not a string.

        """
        if not isinstance(newid,str):
            raise TypeError("New sequence ID 'newid' should be a string.")
        if not isinstance(newseq,str):
            raise TypeError("New sequence string 'newseq' should be a string.")
        if newid in self.seqs:
            raise KeyError("Key name 'newid' already exists.")

        self.seqs[newid]=newseq

    def delseq(self,delid=None,wipeall=False):
        """Deletes sequence(s) from :class:`~crops.elements.sequence.monomer_sequence`.

        :param delid: ID of sequence to be deleted, defaults to None.
        :type delid: str, optional
        :param wipeall: If True, all the sequences are deleted, defaults to False.
        :type wipeall: bool, optional
        :raises TypeError: If delid is not a string.

        """
        if wipeall:
            self.seqs={}
            self.seqs['mainseq']=''
            return
        if delid is None:
            return

        if not isinstance(delid,str):
            raise TypeError("Sequence ID 'delid' should be a string.")
        if not isinstance(wipeall,bool):
            raise TypeError("Boolean switch 'wipeall' is neither True nor False.")
        if delid=='mainseq':
            self.seqs['mainseq']=''
        else:
            self.seqs.pop(delid)

    def mainseq(self,add=None):
        """Returns or modifies the main sequence.

        :param add: If given main sequence is changed to 'add' sequence, defaults to None.
        :type add: str, optional
        :raises TypeError: If 'add' is given and is not a string.
        :return: If 'add' is None, the main sequence is returned.
        :rtype: str

        """
        if not isinstance(add,str) and add is not None:
            raise TypeError("If included, sequence 'add' should be a string.")
        if add is None:
            return self.seqs['mainseq']
        else:
            self.seqs['mainseq']=add

    def dump(self, out):
        """Writes header and main sequence to a file. If file exists, output is appended.

        :param out: An output filepath (str) or an open file.
        :type out: str, file
        :raises TypeError: If out is neither a string nor an open file.

        """
        if not isinstance(out,str) and not isinstance(out,io.IOBase):
            raise TypeError("Argument 'out' should be a string or a file.")

        if self.oligomer_id() is not None:
            header='>'+self.oligomer_id().upper()+'_'+self.chain_id()
            if self.header() is not None:
                header+=retrieve_id(self.header(),extrainfo=True)
        else:
            if self.header() is not None:
                ids=retrieve_id(self.header())
                header='>'+ids[0].upper()+'_'+self.chain_id()+retrieve_id(self.header(),extrainfo=True)
            else:
                header='>'+self.chain_id()+':CROPS'

        if isinstance(out,io.IOBase):
            out.write(header+'\n')
            lenseq=len(self.seqs['mainseq'])
            nlines=int((lenseq-1)/80)+1
            for n in range(nlines):
                out.write(self.seqs['mainseq'][n*80:(n+1)*80]+'\n')
        else:
            outpath=out
            op='a' if os.path.isfile(outpath) else 'w'
            with open(outpath, op) as out:
                out.write(header+'\n')
                lenseq=len(self.seqs['mainseq'])
                nlines=int((lenseq-1)/80)+1
                for n in range(nlines):
                    out.write(self.seqs['mainseq'][n*80:(n+1)*80]+'\n')

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
            self.seqs['fullseq']=self.seqs['mainseq']

        return len(self.seqs['fullseq'])

    def ngaps(self,seqid='gapseq'):
        """Returns the number of gaps ('-') in a sequence.

        :param seqid: The ID of the sequence containing the gaps, defaults to 'gapseq'.
        :type seqid: str, optional
        :raises TypeError: If seqid is not a string.
        :return: Number of gaps in seqid. If seqid not found, ngaps=0.
        :rtype: int

        """
        if not isinstance(seqid,str):
            raise TypeError("Sequence ID 'seqid' should be a string.")
        n=0
        if seqid in self.seqs:
            for char in self.seqs[seqid]:
                if char=='-':
                    n+=1
        return n

    def ncrops(self,seqid='cropseq', offterminals=False, offmidseq=False):

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
        if not isinstance(seqid,str):
            raise TypeError("Sequence ID 'seqid' should be a string.")

        n=0
        if seqid not in self.seqs:
            return n

        for char in self.seqs[seqid]:
            if char=='+' or char=='*':
                n+=1

        if (offterminals==False and offmidseq==False) or (offterminals==True and offmidseq==True):
            return n
        else:
            nterms=0
            for char in self.seqs[seqid]:
                if char=='+' or char=='*':
                    nterms+=1
                else:
                    break
            for char in reversed(self.seqs[seqid]):
                if char=='+' or char=='*':
                    nterms+=1
                else:
                    break

        if offterminals==False and offmidseq==True:
            return n-nterms
        elif offterminals==True and offmidseq==False:
            return nterms

    def header(self):
        """Returns the header identifying the sequence in a fasta file.

        :return: Fasta format header.
        :rtype: str

        """
        if 'header' not in self.info:
            self.info['header']=None
        return self.info['header']

    def oligomer_id(self):
        """Returns the oligomer ID.

        :return: Oligomer ID.
        :rtype: str

        """
        if 'oligomer_id' not in self.info:
            self.info['oligomer_id']=None
        return self.info['oligomer_id']

    def chain_id(self):
        """Returns the chain ID.

        :return: Chain ID.
        :rtype: str

        """
        if 'chain_id' not in self.info:
            self.info['chain_id']=None
        return self.info['chain_id']

class Sequence:
    """A :class:`~crops.elements.sequence.Sequence` object grouping several chain sequence objects.
    The :class:`~crops.elements.sequence.Sequence` class represents a data structure to hold
    all :class:`~crops.elements.sequence.monomer_sequence` and other
    useful information characterising an oligomer.

    :param seq_id: Sequence identifier (e.g. PDB id), defaults to None.
    :type seq_id: str
    :param imer: Container of several :class:`~crops.elements.sequence.monomer_sequence` objects making up the oligomer, defaults to empty dict.
    :type imer: dict [str, :class:`~crops.elements.sequence.monomer_sequence`], optional
    :param source: Information concerning the source of the :class:`~crops.elements.sequence.Sequence` (e.g. Uniprot), defaults to None.
    :type source: str, optional
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
    _kind='Protein/polynucleotide sequence'
    __slots__ = ['seq_id', 'imer', 'source']
    def __init__(self, seq_id=None, imer=None, source=None):

        if not isinstance(seq_id,str) and seq_id is not None:
            raise TypeError("Sequence ID 'seq_id' should be a string.")
        if not isinstance(imer,dict) and imer is not None:
            raise TypeError("Monomer container 'imer' should be a dictionary.")
        elif isinstance(imer,dict):
            for val in imer.values():
                if not isinstance(val,monomer_sequence):
                    raise TypeError("Monomer container 'imer' should only contain elements of monomer_sequence class.")
        self.seq_id = seq_id
        self.imer = imer if imer is not None else {}
        self.source = source

    def __repr__(self):
        if self.source is not None:
            string=self.source+" "+self._kind+" object: (id="+ self.seq_id+", # chains = "+str(len(self.imer))+")"
        else:
            string=self._kind+" object: (id="+ str(self.seq_id)+", # chains = "+str(len(self.imer))+")"
        return string

    def __iter__(self):
        return iter(self.imer.values())

    def copy(self):
        return copy.copy(self)

    def deepcopy(self):
        return copy.deepcopy(self)

    def purge(self):
        """Clears the :class:`~crops.elements.sequence.Sequence` without deleting the object itself.

        """
        self.seq_id=None
        self.source=None
        self.imer.clear()

    def add_monomer(self, nheader, nseq,  nid=None,forceentry=False):
        """Adds a new :class:`~crops.elements.sequence.monomer_sequence` to the :class:`~crops.elements.sequence.Sequence`.

        :param nheader: Standard .fasta header, starting with ">".
        :type nheader: str
        :param nseq: New sequence.
        :type nseq: str
        :param nid: New chain's identifier, defaults to None.
        :type nid: str, optional
        :param forceentry: Switch to force entry under a new ID if nid already in, defaults to False.
        :type forceentry: bool, optional
        :raises KeyError: When :class:`~crops.elements.sequence.monomer_sequence` ID already in :class:`~crops.elements.sequence.Sequence` and forceentry=False.

        """
        if nid is None:
            nid=retrieve_id(nheader)
            nid=nid[1]
        else:
            nid=[nid]

        for iid in nid:
            if iid in self.imer:
                if forceentry:
                    logging.warning('Chain named '+iid+' already exists in Sequence '+self.seq_id+'. Entry forced in with new name.')
                    while True:
                        iid += 'R'
                        if iid not in self.imer:
                            break
                else:
                    raise KeyError('add_monomer ERROR: Chain named '+iid+' already exists in Sequence '+self.seq_id+".")

            self.imer[iid]=monomer_sequence(chid=iid,seq=nseq,header=nheader)


    def del_monomer(self, nid):
        """Removes the selected :class:`~crops.elements.sequence.monomer_sequence` from the :class:`~crops.elements.sequence.Sequence`.

        :param nid: Doomed chain's identifier.
        :type nid: str
        :raises TypeError: When nid is not a string.

        """
        if not isinstance(nid,str):
            raise TypeError('nid should be a string.')

        if nid in self.imer:
            self.pop(nid)
        else:
            logging.warning('Chain named '+ nid+' not found in Sequence.')

    def write(self, outdir, infix="",single=None):
        """Writes all :class:`~crops.elements.sequence.monomer_sequence` objects to .fasta file.

        :param outdir: Output directory.
        :type outdir: str
        :param infix: Filename tag to distinguish from original input file, defaults to "".
        :type infix: str, optional
        :param single: If not None, only :class:`~crops.elements.sequence.monomer_sequence` with ID given by single is written, defaults to None.
        :type single: str, optional
        :raises FileNotFoundError: Output directory not found.
        :raises TypeError: Argument 'single' should be a string.
        :raises KeyError: Specific single ID not found in :class:`~crops.elements.sequence.Sequence`.
        """

        if not os.path.isdir(outdir):
            raise FileNotFoundError(outdir + ' directory not found.')

        if single is None:
            outpath=os.path.join(outdir,self.seq_id+infix+".fasta")
            for monomer in self.imer:
                monomer.dump(outpath)
        else:
            if not isinstance(single,str):
                try:
                    single=str(single)
                except:
                    raise TypeError('single='+single+' should be a string.')
            if single not in self.imer:
                raise KeyError(single + " not found in Sequence "+ self.seq_id)
            outpath=os.path.join(outdir,self.seq_id+"_"+single+infix+".fasta")
            self.imer[single].dump(outpath)

    def length(self,chain):
        """Returns the length of a certain sequence.

        :param chain: ID of :class:`~crops.elements.sequence.monomer_sequence`.
        :type chain: str
        :raises TypeError: When 'chain' is not a string.
        :raises KeyError: Specific chain not found in :class:`~crops.elements.sequence.Sequence`.
        :return: Length of :class:`~crops.elements.sequence.monomer_sequence`.
        :rtype: int

        """
        if not isinstance(chain,str):
            raise TypeError('chain input must be a string.')
        if chain in self.imer:
            return self.imer[chain].length()
        else:
            raise KeyError(chain+' monomer not found in sequence.')

    def nchains(self):
        """Returns number of :class:`~crops.elements.sequence.monomer_sequence` objects in :class:`~crops.elements.sequence.Sequence`.

        :return: Number of :class:`~crops.elements.sequence.monomer_sequence` objects in :class:`~crops.elements.sequence.Sequence`.
        :rtype: int
        """
        return len(self.imer)
