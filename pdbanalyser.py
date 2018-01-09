# a collection of functions to analyse pdb files
# author: Boris Shilov
import re
from math import sqrt, pow
#below imports used for plotting only
import matplotlib.pyplot as plt
import pandas as pd

def parse_pdb(filename):
    """
    This function parses the PDB file and returns a tuple of three lists - one containing the amino acid names in order,
    one containing the coordinates, and one containing raw line strings for HETATM.
    :param filename: absolute filepath to the the PDB file
    :return: ([str], [[float, float, float]], [str])
    """

    #open file and dump every line as a string into list
    with open(filename) as file:
        PDBfile = [line.rstrip('\n') for line in file]

    #retrieve only the ATOM CA lines
    CAregex = re.compile(r'^ATOM\s*[0-9]*\s*CA')
    CAlist = []
    for i in range(len(PDBfile)):
        if re.findall(CAregex, PDBfile[i]):
            CAlist.append(PDBfile[i])

    #retrieve just the coordinate string
    provisionalcoordinatelist = []
    provisionalcoordinateregex = re.compile(r'(?<=[0-9]\s{6})[0-9].*(?=\s{2}[0-1])')
    for i in range(len(CAlist)):
        provisionalcoordinatelist.append(provisionalcoordinateregex.findall(CAlist[i]))

    #convert coordinate string list into a list of lists of x, y, z float coordinates
    provisionalcoordinatelist = [item for sublist in provisionalcoordinatelist for item in sublist]
    coordinatelist = [item.split() for item in provisionalcoordinatelist]
    for coordinate in coordinatelist:
        coordinate[0:3] = map(float, coordinate[0:3])

    #retrieve the amino acid names in order
    residuelist = []
    residueregex = re.compile(r'(?<=CA\s\s)[A-Z]{3,4}(?=\s*A)')
    residueregex2 = re.compile(r'(?<=CA\s)[A-Z]{3,4}(?=\s*A)')
    for i in range(len(CAlist)):
        residuelist.append(''.join(residueregex.findall(CAlist[i]) or residueregex2.findall(CAlist[i])))

    #filter out any traces of conformation B from both the amino acid name list and the corresponding coordinates
    #replace specially marked A variant amino acids with standard designators
    wrongconformationlist = ['BVAL', 'BARG', 'BASP', 'BLYS', 'BGLU', 'BHIS']
    replacementdict = {'AVAL':'VAL', 'AARG':'ARG', 'AASP':'ASP', 'ALYS':'LYS', 'AGLU':'GLU', 'AHIS':'HIS'}
    for i in reversed(range(len(residuelist))):
        if residuelist[i] in wrongconformationlist:
            del residuelist[i]
            del coordinatelist[i]
        elif residuelist[i] in replacementdict.keys():
            residuelist[i] = replacementdict[residuelist[i]]

    #retrieve the HET lines
    HETregex = re.compile(r'^HETATM')
    HETlist = []
    for i in range(len(PDBfile)):
        if re.findall(HETregex, PDBfile[i]):
            HETlist.append(PDBfile[i])
    return coordinatelist, residuelist, HETlist

def centercoordinate(coordinatelist):
    """
    Computes the center of the protein. Returns a protein center vector.
    :param coordinatelist: list of ordered amino acid coordinates
    :param aminoacidlist: list of ordered amino acids
    :return: [x, y, z]
    """
    protein_center = [0, 0, 0]
    for i in range(len(coordinatelist)):
        protein_center[0] = (protein_center[0] + coordinatelist[i][0])
        protein_center[1] = (protein_center[1] + coordinatelist[i][1])
        protein_center[2] = (protein_center[2] + coordinatelist[i][2])
    protein_center[0] = protein_center[0] / len(coordinatelist)
    protein_center[1] = protein_center[1] / len(coordinatelist)
    protein_center[2] = protein_center[2] / len(coordinatelist)
    return protein_center

def datawriter(aminoacidlist, distance_to_center, outputname):
    """
    Writes down a tab separated file with three columns: amino acid, hydrophobicity, distance to the center of the protein.
    :param aminoacidlist: ordered list of amino acids
    :param distance_to_center: ordered list of distances to the center of the protein
    :param outputname: name of the file to be written
    :return: void
    """
    with open(outputname, 'w') as file:
        file.write("AA\tHP\tDIST\n")
        for i in range(len(distance_to_center)):
            if hydrophobechecker(aminoacidlist[i]):
                file.write(aminoacidlist[i] + "\t" + '1' + '\t'
                           + str(distance_to_center[i])
                           + '\n')
            else:
                file.write(aminoacidlist[i] + "\t" + '0' + '\t'
                           + str(distance_to_center[i])
                           + '\n')

def hydrophobechecker(residue):
    """
    Checks if a residue is hydrophobic.
    :param residue: str
    :return: boolean
    """
    hydrophobelist = ['ALA', 'CYS', 'PHE', 'ILE', 'LEU', 'MET', 'PRO', 'VAL', 'TRP']
    if residue in hydrophobelist:
        return True
    else:
        return False

def distance2center(coordinatelist):
    """
    Computes distance to the center of the protein based on the locations of the alpha carbon atoms.
    :param coordinatelist: ordered list of coordinates of the carbon atoms
    :param centercoordinate: center of the protein, a list [x, y, z]
    :return: ordered list of distances to the center
    """
    protein_center = centercoordinate(coordinatelist)
    distance_to_center = []
    for i in range(len(coordinatelist)):
        distance_to_center.append(sqrt(
                                    pow((protein_center[0] - coordinatelist[i][0]), 2)
                                  + pow((protein_center[1] - coordinatelist[i][1]), 2)
                                  + pow((protein_center[2] - coordinatelist[i][2]), 2)))
    return distance_to_center

def distance2atom(coordinatelist, atomcoordinate):
    """
    Computes distances to a given atomic coordinate vector.
    :param coordinatelist: ordered list of coordinate vectors
    :param atomcoordinate: ordered list of amino acids corresponding to the coordinate
    :return: distances to the given atom as a list of lists
    """
    distance_to_atom = []
    for i in range(len(coordinatelist)):
        distance_to_atom.append(sqrt(
            pow((atomcoordinate[0] - coordinatelist[i][0]), 2)
            + pow((atomcoordinate[1] - coordinatelist[i][1]), 2)
            + pow((atomcoordinate[2] - coordinatelist[i][2]), 2)))
    return distance_to_atom

def closestfiveneighbours(aminoacidlist, distance_to_atom):
    """
    Returns the five atoms with the least computed distance
    :param aminoacidlist: ordered list of strings of amino acid names
    :param distance_to_atom: ordered list of computed distances to a given atom (distance2atom computes this)
    :return: tuple of the list of strings of top amino acids and ordered list of distances
    """
    topfivedistance =[]
    topfiveamino = []
    while len(topfivedistance) < 5:
        localminindex = distance_to_atom.index(min(distance_to_atom))
        topfivedistance.append(distance_to_atom[localminindex])
        topfiveamino.append(aminoacidlist[localminindex])
        del distance_to_atom[localminindex], aminoacidlist[localminindex]
    return topfiveamino, topfivedistance


def search_heterogen(HET, atomname):
    """
    Searches a list of strings for a HETATM with a given identity (Such as 'Fe'), returns the list of all coordinates
    where it is found
    :param HET: list of strings of HETATM lines from a PDB file
    :param atomname: string, name of atom to be searched for
    :return: list of coordinates where this atom is found
    """
    coordinateregex = re.compile(r'(?<=[0-9]\s{6})[0-9].*(?=\s{2}[0-1])')
    coordinatelist =[]
    for i in range(len(HET)):
        if re.match(r'^HETATM.[0-9]{4}\s*' + atomname, HET[i]):
            coordinatelist.append(coordinateregex.findall(HET[i]))
    coordinatelist = [item for sublist in coordinatelist for item in sublist]
    coordinatelist = [item.split() for item in coordinatelist]
    for coordinate in coordinatelist:
        coordinate[0:3] = map(float, coordinate[0:3])
    return coordinatelist

def plotter(filename, outputfilename):
    """
    Uses matplotlib to plot a histogram of two distributions of coordinates read from file.
    :param filename:
    :return: void
    """
    dataset = pd.read_csv(filename, sep='\t')
    pd.options.mode.chained_assignment = None
    fig = plt.figure(1)
    hydrophobic = dataset[dataset['HP'] == 1]
    hydrophobic.drop('HP', 1, inplace=True)
    ax1 = plt.subplot(211)
    ax1.set_title("Hydrophobic")
    hydrophobic.plot(kind='hist', ax=plt.gca())
    nothydrophobic = dataset[dataset['HP'] == 0]
    nothydrophobic.drop('HP', 1, inplace=True)
    plt.xlabel('Distance to center')
    ax2 = plt.subplot(212)
    ax2.set_title("Non-hydrophobic")
    nothydrophobic.plot(kind='hist', ax=plt.gca())
    plt.xlabel('Distance to center')
    fig.tight_layout()
    plt.savefig(outputfilename)
