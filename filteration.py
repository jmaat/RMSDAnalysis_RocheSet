import numpy as np
from scipy.spatial import distance
#import qcportal as ptl

import os
from openforcefield.topology import Molecule
from openforcefield.topology import Topology
from openforcefield.typing.engines.smirnoff import ForceField
from openeye import oechem
import copy
import matplotlib.pyplot as plt
from openeye.oechem import *
from openeye.oeomega import * # conformer generation
from openeye.oequacpac import * #for partial charge assignment
from simtk import openmm, unit
from simtk.openmm import app
#from utils import openmmTop_to_oemol
#import utils
import i8

#from 1dplotting import *

#import smirnoff99frosst as ff

i8={-0.15158957347357088: 'CC(=O)N(C)c1ccccc1', -1.083962135605259: 'CC(=O)N(C)c1ccccc1', -0.12184230953250469: 'CC(=O)N(C)c1ccccn1', -0.02427922149511097: 'CC(=O)N(C)c1ccccn1', -0.005664215290238139: 'CC(=O)N(C)c1ccccn1', -0.002175025757409299: 'CC(=O)N(C)c1ccccn1', 0.24023525790319067: 'CC(=O)N(C)c1ccccn1', 0.08776204144789335: 'CC(=O)N(C)c1ncccn1', -0.01813144317786136: 'CC(=O)N(C)c1ncccn1', 0.2563693837725532: 'CC(=O)Nc1c(cccn1)Cl', -0.33124174699087394: 'CC(=O)Nc1ccccc1', -0.3725960200400711: 'CC(=O)Nc1ccccc1Cl', -0.10717703380372434: 'CC(=O)Nc1ccccn1', -0.4692860623495799: 'CC(=O)Nc1ccccn1', -0.5922167617622718: 'CC(=O)Nc1ncccn1', -0.32885369012781607: 'CC(C)(C)Nc1ccccn1', 0.004852063718096933: 'CC(C)Nc1ccccn1', -0.018229273844392907: 'CC1=NCC(=O)N(c2c1cccc2)C', -0.016691681019886218: 'CC1=NCC(=O)N(c2c1cccc2)C', -0.07950034615859587: 'CC1=NCC(=O)Nc2c1cccc2', -0.061172046559091774: 'CCNc1ccccn1', 0.04552495889220198: 'CC[N@@](C)c1ccccn1', 0.035352371135619266: 'CC[N@](C)c1ccccn1', 0.04560332585100407: 'CC[N@](C)c1ccccn1', -0.0017135346510105329: 'CN(C)c1ccccc1', 0.03780796139164366: 'CN(C)c1ccccc1C(=O)c2ccccc2', 0.056853682810940276: 'CN(C)c1ccccc1C(=O)c2ccccc2', 0.031214537910813134: 'CN(C)c1ccccc1c2ccccn2', 0.10509150421955471: 'CN(C)c1ccccc1c2ccccn2', -0.0010151142464005676: 'CN(C)c1ccccn1', -1.514269147878038: 'CN(c1ccccc1)C(=O)c2ccccc2', -1.7321453449725053: 'CN(c1ccccc1)C(=O)c2ccccc2', -0.04980826907593647: 'CN(c1ccccc1)C(=O)c2ccccc2', -0.045005996322098035: 'CN(c1ccccc1)C(=O)c2ccccc2', -0.058689819989384856: 'CN(c1ccccc1)c2ccccc2', -0.13874454021339566: 'CN(c1ccccc1)c2ccccn2', -1.1749998789398237: 'CN(c1ccccc1)c2ccccn2', -0.1284274891501707: 'CN(c1ccccc1)c2ccccn2', -0.21923825492072035: 'CN(c1ccccc1)c2nccs2', -0.0639381298091195: 'CN(c1ccccc1)c2nccs2', -0.06653555381795029: 'CN(c1ccccn1)C(=O)c2ccccc2', 0.07157641048024588: 'CN(c1ccccn1)C(=O)c2ccccc2', 0.3970027323468929: 'CN(c1ccccn1)C(=O)c2ccccc2', -0.030570322128534222: 'CN(c1ccccn1)C(=O)c2ccccc2', 0.041531979473580005: 'CN(c1ccccn1)C(=O)c2ccccc2', -0.07330381974720435: 'CN(c1ccccn1)C(=O)c2ccccc2', -0.0807629221189472: 'CN(c1ccccn1)C(=O)c2ccccc2', -0.003255027418607337: 'CN(c1ccccn1)c2nccs2', -0.13622154182492974: 'CN(c1ccccn1)c2nccs2', -0.23226662710952528: 'CN(c1ccccn1)c2nccs2', -0.17004231673319803: 'CN(c1ccccn1)c2nccs2', -0.1563204691314931: 'CN(c1ccccn1)c2nccs2', -0.003259853944416205: 'CN(c1ccccn1)c2nccs2', -0.0006561966501717198: 'CN1c2ccccc2CCCC1=O', -0.09517763272283523: 'CNc1ccccc1', -0.0027416473507555425: 'CNc1ccccc1C(=O)OC', -0.24414527780938736: 'CNc1ccccc1C(=O)OC', 0.23892009040581008: 'CNc1ccccc1C(=O)c2ccccc2', -0.5232779773771629: 'CNc1ccccc1C(=O)c2ccccc2', -0.11073206072239687: 'CNc1ccccc1F', -0.1286504475102871: 'CNc1ccccc1O', 0.05933826681123425: 'CNc1ccccc1c2ccccn2', 0.11531935664114965: 'CNc1ccccc1c2ccccn2', -0.12573414763508778: 'CNc1ccccc1c2ccccn2', -0.13369298678880584: 'CNc1ccccn1', -0.2419240729496695: 'C[C@@H](c1ccccc1)Nc2ccccc2', -0.11103268604936556: 'C[C@@H](c1ccccc1)Nc2ccccc2', -0.18717287636969512: 'C[C@@H](c1ccccc1)Nc2ccccc2', -0.2115866596651615: 'C[C@@H](c1ccccc1)Nc2ccccc2', -0.11941407343025694: 'C[C@@H]1CCC[N@]1c2ccccc2', 0.048890958460085565: 'C[C@@H]1C[NH2+]CC[N@@]1c2ccccc2', 0.4839925493202816: 'C[C@@H]1C[NH2+]CC[N@@]1c2ccccc2', -0.12454374373546728: 'CN(Cc1ccccc1)c2ccccc2', 0.7825138430824782: 'C[N@@](Cc1ccccc1)c2ccccc2', -0.11771425371082234: 'C[N@@](Cc1ccccc1)c2ccccc2', 0.05479321081724786: 'Cc1cccc(c1N2CC[NH2+]CC2)C', -0.009960657293851982: 'Cc1cccc(c1NC)C', -0.0366622557399864: 'Cc1cccc(c1NC)F', 0.03891546367548432: 'Cc1cccc(c1NC)O', 0.002507671525530511: 'Cc1cccc(c1NC)OC', -0.043185774215065253: 'Cc1ccccc1C2=NCC(=O)Nc3c2cccc3', -0.011319215233398544: 'Cc1ccccc1C2=NCC(=O)Nc3c2cccc3', -0.04607561786427802: 'Cc1ccccc1C2=NCC(=O)Nc3c2cccc3', -0.007723267298359371: 'Cc1ccccc1C2=NCC(=O)Nc3c2cccc3', -0.16939117123486097: 'Cc1ccccc1CNc2ccccc2', -0.0541021820717984: 'Cc1ccccc1CNc2ccccc2', -0.02267903557289075: 'Cc1ccccc1CNc2ccccc2', -0.35092009017135034: 'Cc1ccccc1CNc2ccccc2', 1.4368176763250806: 'Cc1ccccc1CNc2ccccc2', 1.4412573192813427: 'Cc1ccccc1CNc2ccccc2', -0.18182015971364268: 'Cc1ccccc1N2CCCC2', 0.020613176300380565: 'Cc1ccccc1N2CC[NH2+]CC2', -0.12933775933664024: 'Cc1ccccc1NC', 0.006199491146322511: 'Cc1ccccc1NC(=O)c2ccccc2', -0.25075838182535554: 'Cc1ccccc1NC(=O)c2ccccc2', -0.26487665433574836: 'Cc1ccccc1NC(=O)c2ccccc2', 0.07619240472267685: 'c1ccc(c(c1)N2CCCC2)F', 0.003306183801560003: 'c1ccc(c(c1)c2ccccn2)N', 0.019373990172065314: 'c1ccc(c(c1)c2ccccn2)N', 0.05143122813040174: 'c1ccc(cc1)C(=O)Nc2ccccc2', -0.1621322952533456: 'c1ccc(cc1)C(=O)Nc2ccccc2O', -0.15122113916216628: 'c1ccc(cc1)C(=O)Nc2ccccc2O', -0.26971584633327916: 'c1ccc(cc1)C(=O)Nc2ccccc2O', -0.44901272317135715: 'c1ccc(cc1)C(=O)Nc2ccccn2', 0.03468044376121188: 'c1ccc(cc1)C(=O)c2ccccc2N', 0.033990744076082546: 'c1ccc(cc1)C(=O)c2ccccc2N', 0.22944590823622274: 'c1ccc(cc1)C2(CC2)Nc3ccccc3', 0.6842023334777197: 'c1ccc(cc1)C2(CC2)Nc3ccccc3', 0.02333794889751756: 'c1ccc(cc1)C2(CC2)Nc3ccccc3', 0.03392505510298627: 'c1ccc(cc1)C2(CC2)Nc3ccccc3', -0.03156361008415581: 'c1ccc(cc1)C2=NCC(=O)Nc3c2c(ccc3)F', -0.026747935780210258: 'c1ccc(cc1)C2=NCC(=O)Nc3c2c(ccc3)F', -0.01301567856574637: 'c1ccc(cc1)C2=NCC(=O)Nc3c2cccc3', -0.02022079051462411: 'c1ccc(cc1)C2=NCC(=O)Nc3c2cccc3', -0.15274486282325772: 'c1ccc(cc1)CNc2ccccc2', -0.34574847479856574: 'c1ccc(cc1)CNc2ccccc2', -0.36353808513514896: 'c1ccc(cc1)CNc2ccccc2', -0.07299604038617605: 'c1ccc(cc1)N2CCC2', 0.010179111798401808: 'c1ccc(cc1)N2CCCC2', -0.008737131981171098: 'c1ccc(cc1)N2CCCCC2=O', 0.03758080207002765: 'c1ccc(cc1)N2CC[NH2+]CC2', -0.47078851962857404: 'c1ccc(cc1)NC(=O)c2ccccn2', -0.5231866899465725: 'c1ccc(cc1)NC(=O)c2ccccn2', -0.4673318941656156: 'c1ccc(cc1)NC(=O)c2ncccn2', -0.3363467298255647: 'c1ccc(cc1)NCc2ccccc2F', -0.09637260084113075: 'c1ccc(cc1)NCc2ccccc2F', 1.2673709568125893: 'c1ccc(cc1)NCc2ccccc2F', 0.08280823825969386: 'c1ccc(cc1)NCc2ccccc2F', 1.2891022843530486: 'c1ccc(cc1)NCc2ccccc2F', -0.05763863573653871: 'c1ccc(cc1)Nc2ccccc2', -0.009800544179793601: 'c1ccc(cc1)Nc2ccccn2', -0.5380648721809155: 'c1ccc(cc1)Nc2ccccn2', -0.001181100055885155: 'c1ccc(cc1)Nc2ncco2', -0.6004621288939342: 'c1ccc(cc1)Nc2ncco2', -0.0033258768508203987: 'c1ccc(cc1)Nc2nccs2', -0.2679173838033128: 'c1ccc(cc1)Nc2nccs2', 0.030328178798145228: 'c1ccc(cc1)c2ccccc2N', 0.09538297630877304: 'c1ccc(cc1)c2ccccc2N3CCCCC3', -0.02140115621745592: 'c1ccc2c(c1)C(=NCC(=O)N2)c3ccccc3Cl', -0.030524723943135906: 'c1ccc2c(c1)C(=NCC(=O)N2)c3ccccc3Cl', -0.022001821660362963: 'c1ccc2c(c1)C(=NCC(=O)N2)c3ccccc3Cl', -0.0375517839079656: 'c1ccc2c(c1)C(=NCC(=O)N2)c3ccccc3Cl', -0.32480847870576496: 'c1ccc2c(c1)CCCC(=O)N2', 0.009637417192828603: 'c1ccc2c(c1)CCc3ccccc3N2', -0.6473497162650346: 'c1ccnc(c1)C(=O)Nc2ccccn2', -0.7728290795908569: 'c1ccnc(c1)C(=O)Nc2ccccn2', -0.12058011424437876: 'c1ccnc(c1)N2CCCC2', -0.6628915760437857: 'c1ccnc(c1)Nc2ncco2', -1.2030302348495876: 'c1ccnc(c1)Nc2ncco2', -0.8531444194763805: 'c1ccnc(c1)Nc2ncco2', -0.6475241201717872: 'c1ccnc(c1)Nc2ncco2', 0.011799886576359046: 'c1ccnc(c1)Nc2nccs2', -0.7158309039723869: 'c1ccnc(c1)Nc2nccs2', -0.7971907047573209: 'c1ccnc(c1)Nc2nccs2', -0.6857453607239624: 'c1ccnc(c1)Nc2nccs2', -0.03710543075835693: 'c1cnc(nc1)N2CCCCC2=O', -0.07808543678860129: 'c1cnc(nc1)N2CC[NH2+]CC2'}




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
                    smilespic+=str(smiles)+" | "

            if cutoff==False:
                if rmsd < bound:
                    smilespic+=str(smiles)+" | "


    return smilespic







def filter(parDict, filtSmiles):
    """
    parDict: Dictionary of RMSD as keys and smiles
    filtSmiles: The smiles strings that need to be removed from the dictionary
    """
    names = ""
    nomatch=""
    finDict={}
    nomatchDict={}
    mol = oechem.OEGraphMol()
    fragSearch=[]
    before=[]
    after=[]
    singlemolbefore={}
    singlemolafter={}

    for frag in filtSmiles:
        fragSearch.append(oechem.OESubSearch(frag))

    for rmsd, smiles in parDict.items():
        #before.append(rmsd)
        singlemolbefore[smiles]=rmsd
        mol_copy= copy.deepcopy(mol)
        oechem.OESmilesToMol(mol_copy, smiles)
        #create substructure search
        count=0
        for search in fragSearch:
            oechem.OEPrepareSearch(mol_copy, search)
            if search.SingleMatch(mol_copy):
                count+=1
                #print("match")
        if count ==0:
            singlemolafter[smiles]=rmsd
            #after.append(rmsd)
            finDict[rmsd]=smiles
            #finDict = {key:val for key, val in parDict.items() if val != smiles }
            names+=str(smiles) + " | "
        if count>0:
            nomatch+=str(smiles) + " | "
            nomatchDict[rmsd]=smiles

    for smiles, rmsd in singlemolafter.items():
        after.append(rmsd)

    for smiles, rmsd in singlemolbefore.items():
        before.append(rmsd)

    f,a = plt.subplots(2,1)
    a[0].hist(after, 34)
    a[0].set_ylabel("Total mols after filteration")

    a[0].set_title("total:" +str(len(after)))
    #a[0].set_title("After filter histogram: " +str(parDict)+ " " + str(len(after)))
    a[1].hist(before, 34)
    a[1].set_title("total:" +str(len(before)))
    a[1].set_xlabel('RMSD Diff: RMSD(QM, origFF)-RMSD(QM, newFF)')
    a[1].set_ylabel("Total mols before filteration")

    a[0].set_ylim(0, 100)
    a[1].set_ylim(0, 100)

    a[1].grid(True)
    a[0].grid(True)
    a[0].set_xlim(-1.7, 1.7)
    a[1].set_xlim(-1.7, 1.7)
    plt.show()

    #print(nomatch)
    #print(nomatchDict)

    return finDict

def writeDict(dictionary, fileName):
    with open(fileName +".txt", 'w') as fid:
        for RMSD, smiles, in dictionary.items():
                fid.write(str(smiles)+ "       " + str(RMSD)+"\n")
#print(i6)
"""
i6dict = filter(i6, ["o", "O"])
print(len(i6dict))
print(len(i6))
writeDict(i6dict, 'i6Filteredrmsds')
writeDict(i6, 'i6_nonfiltered_rmsds')

print("no filteR:")
print(sortparam(i6, 0, cutoff=False))
print(sortparam(i6dict, 0, cutoff=False))
print("no filteR:")
print(sortparam(i8, -.5, cutoff=False))

"""
filterTests={

i8dict = i8
print("2 filteR:")
print(sortparam(i8dict, -.5, cutoff=False))

print("good cases:")
print(sortparam(i8dict, .5, cutoff=True))


i8dict = filter(i8, ['c1ncccn1', 'c2nccs2', 'C(=O)c2ccccc2'])
print("1 filteR:")
print(sortparam(i8dict, -.5, cutoff=False))

i8dict = filter(i8, ['c1ncccn1', 'c2nccs2', 'C(=O)c2ccccc2', 'c2ncco2'])
print("2 filteR:")
print(sortparam(i8dict, -.5, cutoff=False))



i8dict = filter(i8, ['c1ncccn1', 'c2nccs2', 'C(=O)c2ccccc2', 'c2ncco2', 'c2ccccn2'])
print("2 filteR:")
print(sortparam(i8dict, -.5, cutoff=False))


i8dict = filter(i8, ['c1ncccn1', 'c2nccs2', 'C(=O)c2ccccc2', 'c2ncco2', 'c2ccccn2', 'C(=O)N'])
print("2 filteR:")
print(sortparam(i8dict, -.5, cutoff=False))

"""
i8dict = filter(i8, ['o', 'O', 's', 'N2CCCC2', 'N2CCC2'])
print("2 filteR:")
print(sortparam(i8dict, 0, cutoff=False))
print(sortparam(i8dict, 0, cutoff=True))
#i8dict = filter(i8, ['c1ncccn1', 'c2nccs2', 'c2ncco2', 'Cc1ccccc1'])
#i8dict = filter(i8, ['c1ncccn1', 'c2nccs2', 'c2ncco2', 'Cc1ccccc1', 'c1ccccn1'])
#print(len(i8))
#print(len(i8dict))
writeDict(i8dict, 'i8Filteredrmsds')
writeDict(i8, 'i8_nonfiltered_rmsds')
"""
