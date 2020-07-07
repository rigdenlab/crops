# -*- coding: utf-8 -*-

__prog__="CROPS"
__description__="Cropping and Renumbering Operations for PDB structure and Sequence files"
__author__ = "J. Javier Burgos-Mármol"
__date__ = "May 2020"
__version__ = "0.3.0"

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
    """
    Returns full path if correct.

    Parameters
    ----------
    path : str
        Input (local) path.
    typeofpath : str, optional
        The type of path, 'dir' or 'file'. The default is None.

    Raises
    ------
    ValueError
        When typeofpath is neither 'dir' nor 'file'.
    argparse
        If wrong path given.

    Returns
    -------
    str
        Complete path.

    """
    pathok=False
    if typeofpath=='dir' or typeofpath is None:
       if os.path.isdir(path):
           pathok=True
    elif typeofpath=='file' or typeofpath is None:
        if os.path.isfile(path):
           pathok=True
    else:
        raise ValueError("Input string 'typeofpath' should be either 'dir' or 'file'.")
    if pathok:
        return os.path.abspath(path)
    else:
        raise argparse.ArgumentTypeError(f"readable_dir:{path} is not a valid path")

def target_format(inpath,terms=False,th=0):
    """
    Returns extra information for .fasta headers

    Parameters
    ----------
    inpath : str
        Path to interval database used.
    terms : bool, optional
        Are only terminal ends discarded?. The default is False.
    th : int, float, optional
        Uniprot threshold. The default is 0.

    Returns
    -------
    outcome : TYPE
        DESCRIPTION.

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
    """
    Returns filename tag for outputs.

    Parameters
    ----------
    inpath : str
        Path to interval database used.
    terms : bool, optional
        Are only terminal ends discarded?. The default is False.

    Returns
    -------
    infix_out : str
        Filename tag.

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
    """
    Imports intervals database.

    Parameters
    ----------
    inpath : str
        Path to interval database used.
    pdb_in : str, list, optional
        Chain identifier. The default is None.

    Raises
    ------
    TypeError
        When pdb_in is given and is neither a string nor a list of strings.

    Returns
    -------
    database_out : dict
        A list of :obj:`~crops.core.intervals.intinterval`.

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
    """
    Returns dictionary containing gemmi structures and another one with the file names.

    Parameters
    ----------
    STR_INPATH : str
        Either a directory or file path.

    Raises
    ------
    KeyError
        More than one structure file containing same identifier.

    Returns
    -------
    strdict : dict
        A dictionary containing imported gemmi structures.
    filedict : dict
        A dictionary containing file names.

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
    """
    Returns dictionary containing :obj:`~crops.core.sequence.Sequence`.

    Parameters
    ----------
    inpath : str
        File path.

    Returns
    -------
    newseqs : dict
        A dictionary containing parsed :obj:`~crops.core.sequence.Sequence`.

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
    """
    Returns the desired output filepath.

    Parameters
    ----------
    globaldir : str
        General output dir.
    subdir : str, optional
        Additional subdirectory. The default is None.
    filename : str, optional
        File name. The default is None.
    mksubdir : bool, optional
        Create directory if not existing. The default is False.

    Raises
    ------
    FileNotFoundError
        Directory does not exist.

    Returns
    -------
    newpath : str
        Output filepath.

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
        