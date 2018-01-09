#!/usr/local/bin/python3
# analysis of an example pdb file
# author: Boris Shilov
import sys
from pdbanalyser import *

def main():
    """
    Computes some tasks.
    """
    pdbfile = sys.argv[1]

    coordinatelist, residuelist, HET = parse_pdb(pdbfile)
    print('%s file successfully gutted.' % pdbfile)

    centdistances = distance2center(coordinatelist)

    centerdistancesfile = 'distance2centeroutput.tab'
    datawriter(residuelist, centdistances, centerdistancesfile)
    print('Distances from center of protein for all residues dumped in file %s.' % centerdistancesfile)

    plotoutput='hydrodistributions.pdf'
    plotter(centerdistancesfile, plotoutput)
    print('Histogram of the distribution of hydrophobic and non-hydrophobic residues by distance dumped in %s.' % plotoutput)

    iron = 'FE'
    atomcoords = search_heterogen(HET, iron)
    atomdistances = distance2atom(coordinatelist, atomcoords[0])
    topfivenames, topfivedistances = closestfiveneighbours(residuelist, atomdistances)

    print('For the heterogen atom %s the nearest five neighbours have been calculated to be: ' % iron)
    for i in range(len(topfivenames)):
        if hydrophobechecker(topfivenames[i]):
            print('%s, hydrophobic, at a distance of %s' % (topfivenames[i], topfivedistances[i]))
        else:
            print('%s, at a distance of %s' % (topfivenames[i], topfivedistances[i]))

main()
