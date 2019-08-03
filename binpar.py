import numpy as np
from scipy.spatial import distance
#import qcportal as ptl

import os
from openforcefield.topology import Molecule
from openforcefield.topology import Topology
from openforcefield.typing.engines.smirnoff import ForceField
from openeye import oechem
import copy

from openeye.oechem import *
from openeye.oeomega import * # conformer generation
from openeye.oequacpac import * #for partial charge assignment
from simtk import openmm, unit
from simtk.openmm import app
#from utils import openmmTop_to_oemol
#import utils


#from 1dplotting import *

#import smirnoff99frosst as ff


def get_energy(system, positions):
    """
    Return the potential energy.

    Parameters
    ----------
    system : simtk.openmm.System
        The system to check
    positions : simtk.unit.Quantity of dimension (natoms,3) with units of length
        The positions to use
    Returns
    ---------
    energy
    """

    integrator = openmm.VerletIntegrator(2.0 * unit.femtoseconds)
    context = openmm.Context(system, integrator)
    context.setPositions(positions)
    state = context.getState(getEnergy=True)
    energy = state.getPotentialEnergy().in_units_of(unit.kilocalories_per_mole)
    return energy

def minimize(system, positions, topology):
    # Even though we're just going to minimize, we still have to set up an integrator, since a Simulation needs one
    integrator = openmm.VerletIntegrator(2.0*unit.femtoseconds)
    # Prep the Simulation using the parameterized system, the integrator, and the topology
    simulation = app.Simulation(topology, system, integrator)
    # Copy in the positions
    simulation.context.setPositions( positions)

    # Get initial state and energy; print
    state = simulation.context.getState(getEnergy = True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    print("Energy before minimization (kcal/mol): %.2g" % energy)

    # Minimize, get final state and energy and print
    simulation.minimizeEnergy()
    state = simulation.context.getState(getEnergy=True, getPositions=True)
    energy = state.getPotentialEnergy() / unit.kilocalories_per_mole
    print("Energy after minimization (kcal/mol): %.2g" % energy)
    newpositions = state.getPositions()
    newpositions = newpositions.value_in_unit(unit.angstrom)

    return energy, newpositions

def minimizeStructure(oemol, ff, tag):
        #prep both force fields, create system and topology, get MM energies and store in data tag
    #read in original oemol
    working_copy = copy.deepcopy(oemol)

    molecule = Molecule.from_openeye(oemol, allow_undefined_stereo=True)

    topology = Topology.from_molecules([molecule])

    #print(topology)
    #ff = ForceField(FF)
    system = ff.create_openmm_system(topology)

    positions = molecule.conformers[0]

    #positions = molecule.getPositions()
    #simulation.context.getPositions()

    #positions = extractPositionsFromOEMol(oemol)
    #print("we are minimizing with FF:" + str(FF))
    result = minimize(system, positions, topology)

    # Write out a PDB
    #print('result1', result[1], len(result[1]))
    #print('max atom', oemol.GetMaxAtomIdx())
    out_pos = []
    for i in result[1]:
        for coord in i:
            out_pos.append(coord)
    #print('out_pos', out_pos, len(out_pos))
    working_copy.SetCoords(out_pos)

    #outmol = openmmTop_to_oemol( topology, result[1])


    #print("this is the final MM energy: " + str(result[0]))
    #working_copy.SetData(tag, result[0])
    #oemol.SetData(tag, result[0]._value)


    return working_copy


def GetMM(oemol, FF, tag):
    """
    Description:
    Takes in an oemol and calculates the MM energies with two forcefields, one which has removed
    nitrogen improper parameters. The data from these calculations is stored in the oemol object
    as MM and MMRmNit. The units of energy are KCal/mol.

    Input:
    oemol: A single oemol object
    FF: .offxml file of smirnoff99Frosst.offxml
    tag: The name to tag the energy as, ex: "MM"

    Return:
    oemol: An oemol with the MM energies stored in tag in units kcal/mol
    """

    #prep both force fields, create system and topology, get MM energies and store in data tag
    #read in original oemol
    working_copy = copy.deepcopy(oemol)


    molecule = Molecule.from_openeye(oemol, allow)



    topology = Topology.from_molecules([molecule])


    ff = ForceField(FF)
    system = ff.create_openmm_system(topology)


    positions = molecule.conformers[0]


    #positions = molecule.getPositions()
    #simulation.context.getPositions()

    #positions = extractPositionsFromOEMol(oemol)
    energy = get_energy(system, positions)
    print(type(energy))
    print("this is the MM energy: " + str(energy))
    oemol.SetData(tag, energy._value)

    return oemol




def checkimpropers(mol, FF):
    impropers=[]
    working_copy = copy.deepcopy(mol)
    molecule = Molecule.from_openeye(mol, allow_undefined_stereo=True)
    topology = Topology.from_molecules([molecule])

    molecule_force_list = FF.label_molecules(topology)
    for mol_idx, mol_forces in enumerate(molecule_force_list):
        print(f'Forces for molecule {mol_idx}')
        for force_tag, force_dict in mol_forces.items():
            for (atom_indices, parameter) in force_dict.items():
                if "i" in str(parameter.id) and str(parameter.id) != 'i1' and str(parameter.id) != 'i2' and str(parameter.id) != 'i4':
                    impropers.append(str(parameter.id))
    print("This molecule contains these improper parameters:" +str(impropers))
    return impropers



def makeFile(mol, filename):
    """
    mol: oemol
    filename: name of the file, String
    """
    ofile1 = oemolostream(filename+'.mol2')
    OEWriteConstMolecule(ofile1, mol)
    ofile1.close()

def sortparam(paramDict, bound, cutoff=True):
    """
    input: Parameter dictionary of [RMSD]=smiles. RMSD are the key for the dictionary
    bound: The RMSD lower bound for sorting
    cutoff=True or False. True = upper bound, False = lower bound

    return: A string of the molecules above the upper bound, seperated my periods for loading into picto.
    """
    smilespic=str()

    smilesList=[]
    #removing repeats
    for rmsd, smiles in paramDict.items():
        smilesList.append(smiles)

    originals=set(smilesList)

    for rmsd, smiles in paramDict.items():
        if smiles in originals:
            if cutoff==True:
                if rmsd > bound:
                    smilespic+=str(smiles)+"."

            if cutoff==False:
                if rmsd < bound:
                    smilespic+=str(smiles)+"."


    return smilespic

#from openmoltools import *

#from openmoltools.forcefield_generators import generateTopologyFromOEMol

#open input xyz file
ifs = oechem.oemolistream()
ifs.open('moltest.oeb')

#generate empty list of oemols
oemolList=list()



#iterate through oemols and set coordinates to original oemol

#ofile = oemolostream('1.xyz')



#for mol in ifs.GetOEGraphMols():
    #print(mol.GetData('QM'))


    #oemolList.append(oechem.OEGraphMol(mol))
    #print(mol.GetCoords())
    #OEWriteConstMolecule(ofile, mol)
    #ofile.close()


    #print(mol.GetConfs().next())
#print("here")
improved=0
totalmols=0
errormol=0
#print(dir(ifs))


newissue=0
origissue=0


#comparison variables
smalldiff=0
diffcheck1=0
diffcheck2=0
diffcheck3=0
diffcheck4=0

notimproved=0
ismalldiff=0
idiffcheck1=0
idiffcheck2=0
idiffcheck3=0
idiffcheck4=0


oldff = ForceField('smirnoff99Frosst.offxml')
newff = ForceField('finalff.offxml')


compareDict={}
errorList=[]




i5match=0
i6match=0
i7match=0
i8match=0
i9match=0
i10match=0
i11match=0
i12match=0
i5dict={}
i6dict={}
i7dict={}
i8dict={}
i9dict={}
i10dict={}
i11dict={}
i12dict={}

for mol in ifs.GetOEMols():
    print("***************!Start of analysis of new molecule!****************")



    overlay= True
    qmmol1 = copy.deepcopy(mol)
    qmmol2 = copy.deepcopy(mol)
    mol1 = copy.deepcopy(mol)
    mol2 = copy.deepcopy(mol)

    mol3 = copy.deepcopy(mol)


    try:
        imps = checkimpropers(mol, newff)
        smiles = oechem.OEMolToSmiles(mol)
    #switched order of newmol and oldmol...let's see if this changes our RMSD error.
        newmol = minimizeStructure(mol1, newff, 'newFF')
        newmol3=copy.deepcopy(newmol)
        newFFRMSD=oechem.OERMSD(qmmol1, newmol, overlay)
        print("The new FF RMSD:" + str(newFFRMSD))


        oldmol = minimizeStructure(mol2, oldff, 'originalFF')
        oldmol3=copy.deepcopy(oldmol)
        origFFRMSD=oechem.OERMSD(qmmol2, oldmol, overlay)
        print("The orig FF RMSD:" + str(origFFRMSD))





        rmsdcompare= origFFRMSD-newFFRMSD
        imprmsd=str(rmsdcompare)+ "       " +str(imps)
        compareDict[imprmsd]=smiles


        #analyze parameters
        for param in imps:
            print("this is param for test1:" + str(param))
            if str(param)=="i5":
                i5match+=1
                i5dict[rmsdcompare]=smiles
                print(str(param) + " match! total" + str(i5match))
            elif str(param)=="i6":
                i6match+=1
                i6dict[rmsdcompare]=smiles
                print(str(param) + " match! total"+ str(i6match))
            elif str(param)=="i7":
                i7match+=1
                i7dict[rmsdcompare]=smiles
                print(str(param) + " match! total"+ str(i7match))
            elif str(param)=="i8":
                i8match+=1
                i8dict[rmsdcompare]=smiles
                print(str(param) + " match! total"+ str(i8match))
            elif str(param)=="i9":
                i9match+=1
                i9dict[rmsdcompare]=smiles
                print(str(param) + " match! total"+ str(i9match))
            elif str(param)=="i10":
                i10match+=1
                i10dict[rmsdcompare]=smiles
                print(str(param) + " match! total"+ str(i10match))
            elif str(param)=="i11":
                i11match+=1
                i11dict[rmsdcompare]=smiles
                print(str(param) + " match! total"+ str(i11match))
            elif str(param)=="i12":
                i12match+=1
                i12dict[rmsdcompare]=smiles
                print(str(param) + " match! total"+ str(i12match))
            else:
                print("no matches, something is off. :/")
            print("this is param for test2:" + str(param))

    except:
        pass

print("i5:" + str(i5dict))

print("i6:" + str(i6dict))

print("i7:" + str(i7dict))

print("i8:" + str(i8dict))

print("i9:" + str(i9dict))

print("i10:" + str(i10dict))

print("i11:" + str(i11dict))

print("i12:" + str(i12dict))

print("upper 1.0, good cases i8:")
print(sortparam(i8dict, 1.0, True))

print("lower 1.0, bad cases i8:")
print(sortparam(i8dict, -1.0, False))



print("upper 0.5, good cases i8:")
print(sortparam(i8dict, 0.5, True))


print("lower 0.5, bad cases i8:")
print(sortparam(i8dict, -.5, False))


print("upper 1.0, good cases i6:")
print(sortparam(i6dict, 1.0, True))

print("lower 1.0, bad cases i6:")
print(sortparam(i6dict, -1.0, False))


print("upper 0.5, good cases i6:")
print(sortparam(i6dict, 0.5, True))

print("lower 0.5, bad cases i6:")
print(sortparam(i6dict, 0.5, False))
