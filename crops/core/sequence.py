# -*- coding: utf-8 -*-

__prog__="CROPS"
__description__="Cropping and Renumbering Operations for PDB structure and Sequence files"
__author__ = "J. Javier Burgos-Mármol"
__date__ = "May 2020"
__version__ = "0.3.0"

import os
import io
import copy
from warnings import warn

def retrieve_id(seqheader,extrainfo=False):
    """
    Extract sequence IDs from a standard .fasta header

    Parameters
    ----------
    seqheader : str
        Standard .fasta header, starting with ">"
    extrainfo : bool, optional
        If True extra information string is returned, instead of sequence IDs. The default is False.

    Raises
    ------
    ValueError
        If seqheader is not a string.

    Returns
    -------
    TYPE
        A list with the two sequence identifiers (e.g. [pdb ID, chain ID]) or a single string if extrainfo==True.

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
                k=0 if seqheader[j+7]==' ' else 1
                newchid=''
                for jj in range(j+7+k+1,len(seqheader)):
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

class monomer_sequence:
    """A sequence object representing a single chain sequence.
    The monomer sequence class represents a data structure to hold all
    sequences and other useful information characterising the monomer.
    It contains functions to store, manipulate and organise sequences.
    Examples
    --------
    >>> from crops.core import monomer_sequence
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
    
    Attributes
    ----------
    info : dict
        Useful information of the :obj:`~crops.core.sequence.monomer_sequence`.
    seqs : dict
        The set of sequences, including default "mainseq", in :obj:`~crops.core.sequence.monomer_sequence`.
    """
    
    kind='Single chain sequence'
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
        showseq=len(self.seqs['mainseq']) if self.seqs['mainseq']<20 else self.seqs['mainseq'][:10]+'[...]'+self.seqs['mainseq'][len(self.seqs['mainseq'])-10:]
        if 'oligomer_id' in self.info:
            string=self.kind+" object: ( id="+ self.info['oligomer_id']+'_'+self.info['chain_id'] +", seq="+str(showseq)+", length="+str(len(self.seqs['mainseq']))+" )"
        else:
            string=self.kind+" object: ( id="+ self.info['chain_id'] +", seq="+str(showseq)+", length="+str(len(self.seqs['mainseq']))+" )"
        return string

    def __iter__(self):
        return iter(self.seqs.values())

    def copy(self):
        return copy.copy(self)
    
    def deepcopy(self):
        return copy.deepcopy(self)

    def addseq(self,newid,newseq):
        """
        Add sequence to :obj:`~crops.core.sequence.monomer_sequence`.

        Parameters
        ----------
        newid : str
            New sequence's identifier.
        newseq : str
            New sequence.

        Raises
        ------
        TypeError
            If newid is not a string.
        KeyError
            If sequence is not a string.

        Returns
        -------
        None.

        """
        if not isinstance(newid,str):
            raise TypeError("New sequence ID 'newid' should be a string.")
        if not isinstance(newseq,str):
            raise TypeError("New sequence string 'newseq' should be a string.")
        if newid in self.seqs:
            raise KeyError("Key name 'newid' already exists.")
        
        self.seqs[newid]=newseq
        
    def delseq(self,delid=None,wipeall=False):
        """
        Deletes sequence(s) from :obj:`~crops.core.sequence.monomer_sequence`.
        Note:
        
        Parameters
        ----------
        delid : str, optional
            ID of sequence to be deleted. The default is None.
        wipeall : bool, optional
            If True, all the sequences are deleted. The default is False.

        Raises
        ------
        TypeError
            If delid is not a string.
        TypeError
            If wipeall is not a boolean.

        Returns
        -------
        None.

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
        """
        Returns or modifies the main sequence.

        Parameters
        ----------
        add : str, optional
            If included, main sequence is changed to 'add' sequence instead of returned. The default is None.

        Raises
        ------
        TypeError
            If 'add' is not a string.

        Returns
        -------
        str
            If 'add' not included, the main sequence is returned.

        """
        if not isinstance(add,str) and add is not None:
            raise TypeError("If included, sequence 'add' should be a string.")
        if add is None:
            return self.seqs['mainseq']
        else:
            self.seqs['mainseq']=add

    def dump(self, out):
        """
        Writes header and main sequence to a file. If file exists, output is appended.

        Parameters
        ----------
        out : str, file
            An output filepath (str) or an open file.

        Raises
        ------
        TypeError
            If out is neither a string nor an open file.

        Returns
        -------
        None.

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

    @property
    def length(self):
        """
        Returns the length of the main sequence.

        Returns
        -------
        int
            Length of the main sequence.

        """
        return len(self.seqs['mainseq'])
    @property
    def ngaps(self,seqid):
        """
        Returns the number of gaps ('-') in a sequence.

        Parameters
        ----------
        seqid : str
            The ID of the sequence containing the gaps.

        Raises
        ------
        TypeError
            If seqid is not a string.

        Returns
        -------
        n : int
            Number of gaps in seqid.

        """
        if not isinstance(seqid,str):
            raise TypeError("Sequence ID 'seqid' should be a string.")
        n=0
        for char in self.seqs[seqid]:
            if char=='-':
                n+=1
        return n

    @property
    def ncrops(self,seqid):
        """
        Returns the number of cropped elements ('+','*') in a sequence.

        Parameters
        ----------
        seqid : str
            The ID of the sequence containing the cropped elements.

        Raises
        ------
        TypeError
            If seqid is not a string.

        Returns
        -------
        n : int
            Number of cropped elements in seqid.

        """
        if not isinstance(seqid,str):
            raise TypeError("Sequence ID 'seqid' should be a string.")
        n=0
        for char in self.seqs[seqid]:
            if char=='+' or char=='*':
                n+=1
        return n
    
    @property
    def header(self):
        """
        Returns the header identifying the sequence in a fasta file.

        Returns
        -------
        str
            Fasta format header.

        """
        if 'header' not in self.info:
            self.info['header']=None
        return self.info['header']

    @property
    def oligomer_id(self):
        """
        Returns the oligomer ID.

        Returns
        -------
        str
            Oligomer ID.

        """
        if 'oligomer_id' not in self.info:
            self.info['oligomer_id']=None
        return self.info['oligomer_id']

    @property
    def chain_id(self):
        """
        Returns the chain ID.

        Returns
        -------
        str
            Chain ID.

        """
        if 'chain_id' not in self.info:
            self.info['chain_id']=None
        return self.info['chain_id']

class Sequence:
    """A sequence object grouping several chain sequence objects.
    The Sequence class represents a data structure to hold all
     :obj:`~crops.core.sequence.monomer_sequence` and other useful information
    characterising an oligomer.
    Examples
    --------
    >>> from crops.core import Sequence
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
    
    Attributes
    ----------
    seq_id : str
        Sequence identifier
    imer : dict
        Container of several :obj:`~crops.core.sequence.monomer_sequence` making up the oligomer.
    source : str
        Information concerning the source of the :obj:`~crops.core.sequence.Sequence`.
    """
    kind='Protein/polynucleotide sequence'
    __slots__ = ['seq_id', 'imer', 'source']
    def __init__(self, seq_id=None, imer={}, source=None):

        if not isinstance(seq_id,str) and seq_id is not None:
            raise TypeError("Sequence ID 'seq_id' should be a string.")
        if not isinstance(imer,dict):
            raise TypeError("Monomer container 'imer' should be a dictionary.")
        else:
            for val in imer.values():
                if not isinstance(val,monomer_sequence):
                    raise TypeError("Monomer container 'imer' should only contain elements of monomer_sequence class.")
        self.seq_id = seq_id
        self.imer = imer#.copy()
        self.source = source

    def __repr__(self):
        if self.source is not None:
            string=self.source+" "+self.kind+" object: (id="+ self.seq_id+", # chains = "+str(len(self.imer))+")"
        else:
            string=self.kind+" object: (id="+ self.seq_id+", # chains = "+str(len(self.imer))+")"
        return string
    
    def __iter__(self):
        return iter(self.imer.values())

    def copy(self):
        return copy.copy(self)
    
    def deepcopy(self):
        return copy.deepcopy(self)

    def purge(self):
        """
        Empties the :obj:`~crops.core.sequence.Sequence` without deleting the object itself.

        Returns
        -------
        None.

        """
        self.seq_id=None
        self.source=None
        self.imer.clear()
        
    def add_monomer(self, nheader, nseq,  nid=None,forceentry=False): # INCLUDE DEL_MONOMER
        """
        Adds a new :obj:`~crops.core.sequence.monomer_sequence` to the :obj:`~crops.core.sequence.Sequence`

        Parameters
        ----------
        nheader : str
            Standard .fasta header, starting with ">".
        nseq : str
            New sequence.
        nid : str, optional
            New sequence's identifier. The default is None.
        forceentry : bool, optional
            Switch to force entry under a new ID if already in. The default is False.

        Raises
        ------
        KeyError
            When ID already in dictionary of :obj:`~crops.core.sequence.monomer_sequence` and forceentry=False.

        Returns
        -------
        None.

        """
        if nid is None:
            nid=retrieve_id(nheader)
            nid=nid[1]
        else:
            nid=[nid]
        
        for iid in nid:
            if iid in self.imer and not forceentry:
                raise KeyError('add_monomer ERROR: Chain named '+iid+' already exists in Sequence '+self.seq_id+".")
            elif iid in self.imer and forceentry:
                warn('add_monomer WARNING: Chain named '+iid+' already exists in Sequence '+self.seq_id+'. Entry forced in with new name.')
                while True:
                    iid += 'R'
                    if iid not in self.imernamed:
                        break

            newmonomer=monomer_sequence(iid,nseq,nheader)
            self.imer[iid]=newmonomer

    def write(self, outdir, infix="",single=None):
        """
        Write all :obj:`~crops.core.sequence.monomer_sequence` to .fasta file.

        Parameters
        ----------
        outdir : str
            Output directory.
        infix : str, optional
            Mark to distinguish from original input file. The default is "".
        single : str, optional
            If not None, only :obj:`~crops.core.sequence.monomer_sequence` with ID given by single is written. The default is None.

        Raises
        ------
        FileNotFoundError
            Output directory not found.
        TypeError
            Argument 'single' should be a string.
        KeyError
            Specific single ID not found in :obj:`~crops.core.sequence.Sequence`.

        Returns
        -------
        None.

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
            
    @property
    def length(self,chain):
        """
        Returns the length of a certain sequence.

        Parameters
        ----------
        chain : str
            ID of :obj:`~crops.core.sequence.monomer_sequence`.

        Raises
        ------
        TypeError
            When 'chain' is not a string.
        KeyError
            Specific single ID not found in :obj:`~crops.core.sequence.Sequence`.

        Returns
        -------
        int
            Length of :obj:`~crops.core.sequence.monomer_sequence`.

        """
        if not isinstance(chain,str):
            raise TypeError('chain input must be a string.')
        if chain in self.imer:
            return self.imer[chain].length()
        else:
            raise KeyError(chain+' monomer not found in sequence.')
    @property
    def nchains(self):
        """
        Returns number of :obj:`~crops.core.sequence.monomer_sequence` in :obj:`~crops.core.sequence.Sequence`.

        Returns
        -------
        int
            Number of :obj:`~crops.core.sequence.monomer_sequence` in :obj:`~crops.core.sequence.Sequence`.

        """
        return len(self.imer)