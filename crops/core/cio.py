# -*- coding: utf-8 -*-

__prog__="CROPS"
__description__="Cropping and Renumbering Operations for PDB structure and Sequence files"
__author__ = "J. Javier Burgos-Mármol"
__date__ = "Jul 2020"
__version__ = "0.3.1"

import gemmi
import os
import argparse
import csv
import copy
from crops.sequence import Sequence
from crops.sequence import retrieve_id
from crops.intervals import intinterval

#import sys
#from io import StringIO  # Python3


def check_path(path,typeofpath=None):
    """Returns full path if correct.
    
    :param path: Input (local) path.
    :type path: str
    :param typeofpath: The type of path, 'dir' or 'file', defaults to None.
    :type typeofpath: str, optional
    :raises ValueError: When given typeofpath is neither 'dir' nor 'file'.
    :raises argparse: If wrong path given.
    :return: Complete checked path.
    :rtype: str

    """

    pathok=False
    if typeofpath is not None and typeofpath!='file' and typeofpath!='dir':
        raise ValueError("Input string 'typeofpath' should be either 'dir' or 'file'.")
    if typeofpath=='dir' or typeofpath is None:
       if os.path.isdir(path):
           pathok=True
    if typeofpath=='file' or typeofpath is None:
        if os.path.isfile(path):
           pathok=True
    if pathok:
        return os.path.abspath(path)
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

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
    :param terms: Are only terminal ends discarded?, defaults to False.
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
        "crop" : ".crops"+cut,
        "renumber" : ".crops.seq"}
    
    return infix_out

def import_db(inpath,pdb_in=None):
    """Imports intervals database.
    
    :param inpath: Path to interval database used.
    :type inpath: str
    :param pdb_in: Chain ID, defaults to None.
    :type pdb_in: str, list, optional
    :raises TypeError: When pdb_in is given and is neither a string nor a list of strings.
    :return: dict
    :rtype: A dictionary of :obj:`~crops.core.intervals.intinterval`.

    """
    database_out={}
    if isinstance(pdb_in,str):
        pdb_in=[pdb_in]
    
    if isinstance(pdb_in,list):
        for element in pdb_in:
            if not isinstance(element,str):
                raise TypeError('Argument should be a list of strings or a string')
            element=element.lower()

    if os.path.basename(inpath)=='pdb_chain_uniprot.csv':
        mol=0
        chain=1
        up=2
        leftend=3
        rightend=4
    else:
        mol=0
        chain=1
        leftend=2
        rightend=3
        up=None

    csv_chain_file = open(inpath)
    csv_chain = csv.reader(csv_chain_file)   

    for entry in csv_chain:
        if entry[0][0] != "#" and entry[0] !="PDB":
            if pdb_in is None or entry[mol] in pdb_in:
                if entry[mol] not in database_out:
                    database_out[entry[mol]]={}
                if entry[chain] not in database_out[entry[mol]]:
                    database_out[entry[mol]][entry[chain]]=intinterval(description=entry[mol]+'_'+entry[chain])
                    if up is not None:
                        database_out[entry[mol]][entry[chain]].tags['uniprot']={}
                database_out[entry[mol]][entry[chain]].union([int(entry[leftend]),int(entry[rightend])])
                if up is not None:
                    if entry[up] in database_out[entry[mol]][entry[chain]].tags['uniprot']:
                        database_out[entry[mol]][entry[chain]].tags['uniprot'][entry[up]][1]=int(entry[rightend])
                    else:
                        database_out[entry[mol]][entry[chain]].tags['uniprot'][entry[up]]=[int(entry[leftend]),int(entry[rightend])]

    return database_out


def parsestrfile(STR_INPATH):
    """Returns dictionary containing gemmi structures and another one with the file names.
    
    :param STR_INPATH: Either a directory or file path.
    :type STR_INPATH: str
    :raises KeyError: More than one structure file containing same identifier.
    :return strdict: A dictionary containing imported gemmi structures.
    :rtype strdict: dict
    :return filedict: A dictionary containing file names.
    :rtype filedict: dict

    """
    strdict={}
    filedict={}
    if os.path.isfile(STR_INPATH):
        structure=gemmi.read_structure(STR_INPATH)
        pdbid=structure.name.lower()
        strdict[pdbid]=copy.deepcopy(structure)
        filedict[pdbid]=os.path.basename(STR_INPATH)
    elif os.path.isdir(STR_INPATH):
        filelist=os.listdir(STR_INPATH)
        for file in filelist:
            if os.isfile(file):
                try:
                    structure=gemmi.read_structure(file)
                    pdbid=structure.name.lower()
                    if pdbid in strdict:
                        raise KeyError('Structure '+pdbid+' loaded more than once. Check files in directory and remove duplicates.')
                    strdict[pdbid]=copy.deepcopy(structure)
                    filedict[pdbid]=os.path.basename(STR_INPATH)
                except:
                    pass

    return strdict, filedict

def parseseqfile(inpath):
    """Returns dictionary containing :class:`~crops.core.sequence.Sequence`.
    
    :param inpath: File path.
    :type inpath: str
    :return: A dictionary containing parsed :class:`~crops.core.sequence.Sequence`.
    :rtype: dict

    """
    newseqs={}
    newid=[]
    head=''
    chain=''
    with open(inpath,'r') as f:    
        indx=-1
        while True:
            line=f.readline().rstrip()
            if not line:
                try:
                    line=f.readline().rstrip()
                    if not line:
                        break
                    else:
                        if line.startswith(">"):
                            if indx >= 0:
                                if newid[0] not in newseqs:
                                    newseqs[newid[0].lower()]=Sequence(seq_id=newid[0].lower(),source=os.path.basename(inpath))
                                for iid in newid[1]:
                                    newseqs[newid[0].lower()].add_monomer(head,chain,iid)
                            newid=retrieve_id(line)
                            head=line
                            indx += 1
                            chain = ''
                        elif line.startswith("#") or line.startswith(' #'):
                            pass
                        else:
                            chain += str(line)
                except:
                    break
            else:
                if line.startswith(">"):
                    if indx >= 0:
                        if newid[0] not in newseqs:
                            newseqs[newid[0].lower()]=Sequence(seq_id=newid[0].lower(),source=os.path.basename(inpath))
                        for iid in newid[1]:
                            newseqs[newid[0].lower()].add_monomer(head,chain,iid)
                    newid=retrieve_id(line)
                    head=line
                    indx += 1
                    chain = ''
                elif line.startswith("#") or line.startswith(' #'):
                    pass
                else:
                    chain += str(line)

    return newseqs

def outpath(globaldir,subdir=None,filename=None,mksubdir=False):
    """Returns the desired output filepath.
    
    :param globaldir: General output dir.
    :type globaldir: str
    :param subdir: Additional subdirectory, defaults to None.
    :type subdir: str, optional
    :param filename: File name, defaults to None.
    :type filename: str, optional.
    :param mksubdir: Create directory if not existing, defaults to False.
    :type mksubdir: bool, optional
    :raises FileNotFoundError: Directory does not exist.
    :return: Output filepath.
    :rtype: str

    """
    newpath=check_path(globaldir,'dir')
    
    if subdir is not None:
        newpath=os.path.join(newpath,subdir)
        if not os.path.isdir(newpath):
            if mksubdir:
                os.path.mkdir(newpath)
            else:
                raise FileNotFoundError('Directory does not exist')
    if filename is not None:
        newpath=os.path.join(newpath,filename)
    
    return newpath
        