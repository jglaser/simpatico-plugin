from hoomd_plugins.simpatico import _simpatico

from hoomd_script import globals
from hoomd_script import pair
from hoomd_script import bond
from hoomd_script import analyze
from hoomd_script import util
from hoomd_script import data

import hoomd

import sys
import weakref

class diagnostic(analyze._analyzer):
    ## \brief Launches and initializes Simpatico
    #
    # \param diagnostic_params  parameter file section for diagnostics
    def __init__(self, period, diagnostic_params=None, out_dir=None, queue_limit=5*1024^3):
        util.print_status_line();

        # initalize base class
        analyze._analyzer.__init__(self);

        if out_dir == None:
           out_dir = "./";
        self.out_dir = out_dir;
        self.species = self.find_distinct_species();
        self.diagnostic_params = diagnostic_params;

        # create c++ mirror class
        self.cpp_analyzer = _simpatico.Simpatico(globals.system_definition, self.generate_parameters, queue_limit);
        self.setupAnalyzer(period);
        #hack to prevent cylic reference
        self.cpp_analyzer = None

    def find_distinct_species(self):
        # nodes correspond to particles
        # bonds are connections between nodes
        nodes = []; 

        system = data.system_data(globals.system_definition)
        nparticles = len(system.particles)
        for p in range(0,nparticles):
            nodes.append([[], False])
        
        bdata_snapshot = hoomd.SnapshotBondData(len(system.bonds))

        globals.system_definition.getBondData().takeSnapshot(bdata_snapshot)

        for i in range(0,len(bdata_snapshot.bonds)):
            bond = bdata_snapshot.bonds[i]
            type = bdata_snapshot.type_id[i]

            nodes[bond.x][0].append((bond.y, type));
            nodes[bond.y][0].append((bond.x, type));

        # find distinct connected components
        species = {};
       
        for (nodeID, (connections, visited)) in enumerate(nodes):
            if visited == False:
                (this_species, natoms) = self.recursivelyMark(nodeID, nodes)
                this_species = tuple([tuple(sublist) for sublist in this_species])
                if this_species not in species.keys():
                    species[this_species]=[1, natoms+1]
                else:
                    species[this_species][0]+=1;

        return species; 

    ## \internal
    def recursivelyMark(self,nodeID, nodes, local_bonds=None, local_id=None, bond_type=None, local_types=None):
        system = data.system_data(globals.system_definition)
        (connections, visited) = nodes[nodeID]
        if visited:
            return

        connectedNodeID = local_id

        if local_id is None:
            local_id = -1

        local_id += 1

        if not local_types:
            local_types = []

        if not local_bonds:
            local_bonds = []

        if connectedNodeID is not None:
            local_bonds.append((connectedNodeID,local_id,bond_type))
        
        # mark as visited 
        nodes[nodeID][1] = True
         
        local_types.append(system.particles[nodeID].type)
        
        for (next_connectedNodeID,next_bond_type) in connections:
            ret = self.recursivelyMark(next_connectedNodeID, nodes, local_bonds, local_id, next_bond_type, local_types)
            if ret is not None:
                ((local_bonds,local_types),local_id) = ret
        
        return ((local_bonds,local_types),local_id)
 
    ## \internal
    # \brief Generates the parameter file contents
    def generate_parameters(self):
        system = data.system_data(globals.system_definition)
        parameters = "MdSimulation{\n"
        parameters += " FileMaster{\n"
        parameters += "  commandFileName paramfile\n"
        parameters += "  inputPrefix ./\n"
        parameters += "  outputPrefix " + str(self.out_dir) + "/\n"
        parameters += " }\n"
        parameters += " nAtomType "+ str(len(system.particles.types)) + "\n"
        parameters += " nBondType "+ str(system.bonds.bdata.getNBondTypes()) + "\n"
        parameters += " atomTypes "
        for cur_type in system.particles.types: 
            parameters += " " +  cur_type + " 1.0"; # mass 1.0
        parameters += "\n"
        parameters += " maskedPairPolicy"
#FIXME
        parameters += "  MaskNone\n"
#        if (globals.neighbor_list.cpp_nlist.countExclusions() > 0):
#             parameters += " MaskBonded\n"
#        else:
#            parameters += " MaskNone\n"

        parameters += " SpeciesManager{\n"
        for cur_species in self.species.keys():
            parameters += "  Species{\n"
            parameters += "  moleculeCapacity " + str(self.species[cur_species][0]) + "\n"
            parameters += "  nAtom " + str(self.species[cur_species][1]) + "\n"
            parameters += "  nBond " + str(len(cur_species[0])) + "\n"
            parameters += "  atomTypeIds\n"
            for cur_type in cur_species[1]:
                parameters += str(system.particles.types.index(cur_type)) + "\n"
            parameters += "  speciesBonds\n"
            
            for cur_bond in cur_species[0]:
                parameters += str(cur_bond[0]) + " " + str(cur_bond[1]) + " " + str(cur_bond[2]) + "\n"
        parameters += "  }\n"
        parameters += " }\n"
        parameters += " Random{\n"
        parameters += "   seed 0\n"
        parameters += " }\n"
        parameters += " MdSystem{\n"

        # count pair forces
        num_pair_force = 0;   
        for force in globals.forces:
            if force.enabled and isinstance(force, pair.pair):
                num_pair_force+=1
                pair_force = force;       

        if num_pair_force == 0 or num_pair_force > 1:
            print >> sys.stderr, "\n***Error! Only one pair force is supported in Simpatico";
            raise RuntimeError('Error creating Simpatico parameter file');

        parameters += "  pairStyle"
        if (isinstance(pair_force, pair.lj)):
            parameters += " LJPair\n"
            pair_parameters = "   epsilon"
            for typei in system.particles.types:
                for typej in system.particles.types:
                    pair_parameters += " " + str(pair_force.pair_coeff.get(typei,typej,'epsilon'))
            pair_parameters += "\n"
            pair_parameters += "   sigma"
            for typei in system.particles.types:
                for typej in system.particles.types:
                    pair_parameters += " " + str(pair_force.pair_coeff.get(typei,typej,'sigma'))
            pair_parameters += "\n"
            pair_parameters += "   cutoff"
            for typei in system.particles.types:
                for typej in system.particles.types:
                    pair_parameters += " " + str(pair_force.pair_coeff.get(typei,typej,'r_cut'))
            pair_parameters += "\n"
        elif (isinstance(pair_force, pair.dpd_conservative)):
            parameters += " DpdPair\n"
            pair_parameters = "   epsilon"
            for typei in system.particles.types:
                for typej in system.particles.types:
                    pair_parameters += " " + str(pair_force.pair_coeff.get(typei,typej,'A'))
            pair_parameters += "\n"
            pair_parameters += "   sigma"
            for typei in system.particles.types:
                for typej in system.particles.types:
                    pair_parameters += " " + str(pair_force.pair_coeff.get(typei,typej,'r_cut'))
            pair_parameters += "\n"
        else:
            parameters += " LJPair\n"
            pair_parameters = "   epsilon"
            for typei in system.particles.types:
                for typej in system.particles.types:
                    pair_parameters += " 1.0"
            pair_parameters += "\n"
            pair_parameters += "   sigma"
            for typei in system.particles.types:
                for typej in system.particles.types:
                    pair_parameters += " 1.0"
            pair_parameters += "\n"
            pair_parameters += "   cutoff"
            for typei in system.particles.types:
                for typej in system.particles.types:
                    pair_parameters += " 1.0"
            pair_parameters += "\n"
 
            print >> sys.stderr, "\n***Warning! Unsupported pair potential. Substituting default values.";

        # map bond potential
        parameters += "  bondStyle"
        bond_force = None;
        for force in globals.forces:
            if isinstance(force,bond.harmonic) or isinstance(force,bond.fene):
                if bond_force is not None:
                    print >> sys.stderr, "\n***Error! Only one bond force is supported";
                    raise RuntimeError('Error creating Simpatico parameter file');
                bond_force = force

        nbondtypes = globals.system_definition.getBondData().getNBondTypes();
        bond_type_list = [];
        for i in range(0,nbondtypes):
            bond_type_list.append(globals.system_definition.getBondData().getNameByType(i));

        if isinstance(bond_force,bond.harmonic):
            parameters += " HarmonicBond\n"
            bond_parameters = "   kappa"
            for type in bond_type_list:
                bond_parameters += " " + str(bond_force.bond_coeff.get(type,'k'))
            bond_parameters += "\n"
            bond_parameters += "   length"
            for type in bond_type_list:
                bond_parameters += " " + str(bond_force.bond_coeff.get(type,'r0'))
            bond_parameters += "\n"
        else:
            print >> sys.stderr, "\n***Error! Unsupported bond potential";
            raise RuntimeError('Error creating Simpatico parameter file');

        parameters += "  MdPairPotential{\n"
        parameters += pair_parameters
        parameters += "   maxBoundary orthorhombic " + str(system.box[0]) + " " + str(system.box[1]) + " " + str(system.box[2]) + "\n"
        parameters += "  PairList{\n"
        parameters += "    atomCapacity " + str(len(system.particles)) + "\n"
        # adjust neighbors_per_atom if necessary
        neighbors_per_atom = 50
        parameters += "    pairCapacity " + str(len(system.particles)*neighbors_per_atom) + "\n"
        parameters += "    skin  .001\n " # any value > 0 
        parameters += "   }\n"
        parameters += "  }\n"
        parameters += "  BondPotential{\n"
        parameters += bond_parameters
        parameters += "  }\n"

        parameters += "  EnergyEnsemble{\n"
        parameters += "   type isothermal\n" # for post-processing we can choose any energy ensemble
        parameters += "   temperature 1.0\n"
        parameters += "  }\n"
        parameters += "  BoundaryEnsemble{\n"
        parameters += "   type rigid\n" # choose default boundary ensemble
        parameters += "  }\n"
        parameters += "  NvtNhIntegrator{\n"
        parameters += "   dt 0.001\n" # standard values
        parameters += "   tauT 1.0\n" 
        parameters += "  }\n"
        parameters += " }\n"

        # append definition of diagnostic manager
        parameters += " DiagnosticManager{\n"
        parameters += "  baseInterval 1\n"
        if self.diagnostic_params is not None:
             parameters += self.diagnostic_params
        parameters += " }\n"
        parameters += "}\n"

        return parameters
