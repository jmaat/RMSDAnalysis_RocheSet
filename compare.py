import numpy as np
from scipy.spatial import distance
import qcportal as ptl
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

i5dictm={}
i6dictm={}
i7dictm={}
i8dictm={}
i9dictm={}
i10dictm={}
i11dictm={}
i12dictm={}

fullDict={}

for mol in ifs.GetOEMols():
    print("***************!Start of analysis of new molecule!****************")



    overlay= True
    automorph= True
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
        newFFRMSD=oechem.OERMSD(qmmol1, newmol, automorph, overlay)
        print("The new FF RMSD:" + str(newFFRMSD))


        oldmol = minimizeStructure(mol2, oldff, 'originalFF')
        oldmol3=copy.deepcopy(oldmol)
        origFFRMSD=oechem.OERMSD(qmmol2, oldmol, automorph, overlay)

        print("The orig FF RMSD:" + str(origFFRMSD))





        rmsdcompare= origFFRMSD-newFFRMSD
        imprmsd=str(rmsdcompare)+ "       " +str(imps)
        compareDict[imprmsd]=smiles
        smilesImp=[]
        smilesImp.append(smiles)
        smilesImp.append(imps)

        fullDict[rmsdcompare]=smilesImp


        #analyze parameters

        if len(imps)==1:
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
        else:
            for param in imps:
                print("this is param for test1:" + str(param))
                if str(param)=="i5":
                    i5dictm[rmsdcompare]=smiles
                elif str(param)=="i6":
                    i6dictm[rmsdcompare]=smiles
                elif str(param)=="i7":
                    i7dictm[rmsdcompare]=smiles
                elif str(param)=="i8":
                    i8dictm[rmsdcompare]=smiles
                elif str(param)=="i9":
                    i9dictm[rmsdcompare]=smiles
                elif str(param)=="i10":
                    i10dictm[rmsdcompare]=smiles
                elif str(param)=="i11":
                    i11dictm[rmsdcompare]=smiles
                elif str(param)=="i12":
                    i12dictm[rmsdcompare]=smiles
                else:
                    print("no matches, something is off. :/")
                print("this is param for test2:" + str(param))

        totalmols+=1
        smiles = oechem.OEMolToSmiles(mol)
        print("***Comparing the RMSD's***")
        if origFFRMSD>newFFRMSD:
            improved+=1
            print("This many cases have an improved energy:" + str(improved) + " out of total:" + str(totalmols))
            if (origFFRMSD-newFFRMSD)<0.1:
                smalldiff+=1
                print("Small RMSD diff. The difference is less than 0.1, total:" + str(smalldiff))
            if (origFFRMSD-newFFRMSD)>0.5:
                diffcheck3+=1
                print("The difference is greater than .5!!!, total:"+ str(diffcheck3))
                makeFile(oldmol3, "good_HighposDiff_RMSD_"+ str(rmsdcompare)+ "_old_" + str(imps) + "_" + str(smiles))
                makeFile(mol3, "good_HighposDiff_RMSD_"+ str(rmsdcompare)+ "_QM_" + str(imps) + "_" + str(smiles))
                makeFile(newmol3, "good_HighposDiff_RMSD_"+ str(rmsdcompare)+ "_new_" + str(imps) + "_" + str(smiles))
            if (origFFRMSD-newFFRMSD)>0.3:
                diffcheck1+=1
                makeFile(oldmol3, "kindagood_HighposDiff_RMSD_"+ str(rmsdcompare)+ "_old_" + str(imps) + "_" + str(smiles))
                makeFile(mol3, "kindagood_HighposDiff_RMSD_"+ str(rmsdcompare)+ "_QM_" + str(imps) + "_" + str(smiles))
                makeFile(newmol3, "kindagood_HighposDiff_RMSD_"+ str(rmsdcompare)+ "_new_" + str(imps) + "_" + str(smiles))
                print("The difference is greater than 0.3, total:" + str(diffcheck1))


        if origFFRMSD<newFFRMSD:
            notimproved+=1
            print("This many cases do NOT have improved energy:" + str(notimproved) + " out of total:" + str(totalmols))
            if (newFFRMSD-origFFRMSD)<0.1:
                ismalldiff+=1
                print("i Small RMSD diff. The difference is less than 0.1, total:" + str(ismalldiff))
            if (origFFRMSD-newFFRMSD)>0.5:
                diffcheck3+=1
                print("The difference is greater than 1!!!, total:"+ str(diffcheck3))
                makeFile(oldmol3, "bad_HighnegDiff_RMSD_"+ str(rmsdcompare)+ "_old_" + str(imps) + "_" + str(smiles))
                makeFile(mol3, "bad_HighnegDiff_RMSD_"+ str(rmsdcompare)+ "_QM_" + str(imps) + "_" + str(smiles))
                makeFile(newmol3, "bad_HighnegDiff_RMSD_"+ str(rmsdcompare)+ "_new_" + str(imps) + "_" + str(smiles))
            if (origFFRMSD-newFFRMSD)>0.3:
                diffcheck1+=1
                makeFile(oldmol3, "kindabad_HighposDiff_RMSD_"+ str(rmsdcompare)+ "_old_" + str(imps) + "_" + str(smiles))
                makeFile(mol3, "kindabad_HighposDiff_RMSD_"+ str(rmsdcompare)+ "_QM_" + str(imps) + "_" + str(smiles))
                makeFile(newmol3, "kindabad_HighposDiff_RMSD_"+ str(rmsdcompare)+ "_new_" + str(imps) + "_" + str(smiles))
                print("The difference is greater than 0.1, total:" + str(diffcheck1))


    except Exception as e:
        print("There's an error....")
        print(e)
        errorList.append(smiles)
        pass


print("LOOP FINISHED, woo.")
print("The final stats are:")
print("***Analysis for the improved molecules***")
print("This many cases have an improved energy:" + str(improved) + " out of total:" + str(totalmols))
print("Small RMSD diff. The difference is less than 0.1, total:" + str(smalldiff))
print("The difference is greater than 0.1, total:" + str(diffcheck1))
print("The difference is even greater than 0.5, total:" + str(diffcheck2))
print("The difference is greater than 1!!!, total:"+ str(diffcheck3))
print("The difference is greater than 2!!!, total:"+ str(diffcheck4))

print("***Analysis for the non-improved molecules***")
print("This many cases have an improved energy:" + str(notimproved) + " out of total:" + str(totalmols))
print("Small RMSD diff. The difference is less than 0.1, total:" + str(ismalldiff))
print("The difference is greater than 0.1, total:" + str(idiffcheck1))
print("The difference is even greater than 0.5, total:" + str(idiffcheck2))
print("The difference is greater than 1!!!, total:"+ str(idiffcheck3))
print("The difference is greater than 2!!!, total:"+ str(idiffcheck4))

print("-----The final total---------")

print("total molecules:" +  str(totalmols))

with open('rmsds.txt','w') as fid:
    for mol, rmsd in compareDict.items():
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")


with open('errormols.txt','w') as fid2:
    for error in errorList:
        fid2.write(str(error) + "\n")


print("This is the final list of RMSDS:")
print(str(compareDict))


with open('i8.txt','w') as fid:
    fid.write("#The total number of mols using i8:" + str(i8match))
    fid.write("#single param dict: \n" + str(i8dict) + "\n")
    fid.write("#multi param dict: \n" + str(i8dictm) + "\n")
    for mol, rmsd in i8dict.items():
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")

with open('i6.txt','w') as fid:
    fid.write("#The total number of mols using i6:" + str(i6match))
    fid.write("#single param dict: \n" + str(i6dict) + "\n")
    fid.write("#multi param dict: \n " + str(i6dictm) + "\n" )
    for mol, rmsd in i6dict.items():
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")

import pickle

single={}
single['i5']=i5dict
single['i6']=i6dict
single['i7']=i7dict
single['i8']=i8dict
single['i9']=i9dict
single['i10']=i10dict
single['i11']=i11dict
single['i12']=i12dict

multi={}
multi['i5']=i5dictm
multi['i6']=i6dictm
multi['i7']=i7dictm
multi['i8']=i8dictm
multi['i9']=i9dictm
multi['i10']=i10dictm
multi['i11']=i11dictm
multi['i12']=i12dictm

with open('rmsds.pickle','wb') as fid:
    all_dicts = {'all': fullDict, 'single':single, 'multi':multi}
    pickle.dump(all_dicts, fid)


#with open('i6.pickle','rb') as fid:
    #mynewdict=pickle.load(fid)


"""
with open('i5.txt','w') as fid:
    fid.write("The total number of mols using i5:" + str(i5match))
    for mol, rmsd in i5dict.items():
        fid.write("single param dict:" + str(i5dict))
        fid.write("multi param dict:" + str(i5dictm))
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")

with open('i6.txt','w') as fid:
    fid.write("The total number of mols using i6:" + str(i6match))
    for mol, rmsd in i6dict.items():
        fid.write("single param dict:" + str(i6dict))
        fid.write("multi param dict:" + str(i6dictm))
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")

with open('i7.txt','w') as fid:
    fid.write("The total number of mols using i7:" + str(i7match))
    for mol, rmsd in i7dict.items():
        fid.write("single param dict:" + str(i7dict))
        fid.write("multi param dict:" + str(i7dictm))
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")

with open('i8.txt','w') as fid:
    fid.write("The total number of mols using i8:" + str(i8match))
    for mol, rmsd in i8dict.items():
        fid.write("single param dict:" + str(i8dict))
        fid.write("multi param dict:" + str(i8dictm))
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")

with open('i9.txt','w') as fid:
    fid.write("The total number of mols using i9:" + str(i9match))
    for mol, rmsd in i9dict.items():
        fid.write("single param dict:" + str(i9dict))
        fid.write("multi param dict:" + str(i9dictm))
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")


with open('i10.txt','w') as fid:
    fid.write("The total number of mols using i10:" + str(i10match))
    for mol, rmsd in i10dict.items():
        fid.write("single param dict:" + str(i10dict))
        fid.write("multi param dict:" + str(i10dictm))
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")

with open('i11.txt','w') as fid:
    fid.write("The total number of mols using i11:" + str(i11match))
    for mol, rmsd in i11dict.items():
        fid.write("single param dict:" + str(i11dict))
        fid.write("multi param dict:" + str(i11dictm))
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")

with open('i12.txt','w') as fid:
    fid.write("The total number of mols using i12:" + str(i12match))
    for mol, rmsd in i12dict.items():
        fid.write("single param dict:" + str(i12dict))
        fid.write("multi param dict:" + str(i12dictm))
        fid.write(str(mol)+ "       " + str(rmsd)+"\n")

with open('param_total.txt', 'w') as fid:
    fid.write("The total number of mols using i5:" + str(i5match)+"\n")
    fid.write("The total number of mols using i6:" + str(i6match)+"\n")
    fid.write("The total number of mols using i7:" + str(i7match)+"\n")
    fid.write("The total number of mols using i8:" + str(i8match)+"\n")
    fid.write("The total number of mols using i9:" + str(i9match)+"\n")
    fid.write("The total number of mols using i10:" + str(i10match)+"\n")
    fid.write("The total number of mols using i11:" + str(i11match)+"\n")
    fid.write("The total number of mols using i12:" + str(i12match)+"\n")

"""

"""
        if newFFRMSD == -1:
            newissue+=1
            print("Issue with newFFRMSD")
            print("This many cases also have issues with the newFF:"+str(newissue))
            smiles = oechem.OEMolToSmiles(mol)
            print("there's a case of -1 RMSD!!!!! .mol2 files are being generated for molecule:" + str(smiles))
            smiles = oechem.OEMolToSmiles(mol)
            ofile1 = oemolostream(str(smiles)+'_qm.mol2')
            ofile2 = oemolostream(str(smiles)+'_originalFF.mol2')
            ofile3 = oemolostream(str(smiles)+'_newFF_issuewith.mol2')
            OEWriteConstMolecule(ofile1, qmmol)
            OEWriteConstMolecule(ofile2, oldmol)
            OEWriteConstMolecule(ofile3, newmol)
            ofile1.close()
            ofile2.close()
            ofile3.close()


        if origFFRMSD == -1:
            origissue+=1
            print("Issue with origFFRMSD")
            print("This many cases also have issues with the origFF:"+str(origissue))

            print("there's a case of -1 RMSD!!!!! .mol2 files are being generated for molecule:" + str(smiles))
            smiles = oechem.OEMolToSmiles(mol)
            ofile1 = oemolostream(str(smiles)+'_qm.mol2')
            ofile2 = oemolostream(str(smiles)+'_originalFF_issuewith.mol2')
            ofile3 = oemolostream(str(smiles)+'_newFF.mol2')
            OEWriteConstMolecule(ofile1, qmmol)
            OEWriteConstMolecule(ofile2, oldmol)
            OEWriteConstMolecule(ofile3, newmol)
            ofile1.close()
            ofile2.close()
            ofile3.close()
"""
