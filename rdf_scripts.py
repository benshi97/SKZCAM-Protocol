from ase import io
from ase.units import Bohr
from ase.visualize import view
from ase.neighborlist import natural_cutoffs
import numpy as np
import matplotlib.pyplot as plt



def get_test_cluster(filename='',output_file_append='',M_element='Mg',O_element='O',centre_atom='O',centre_atom_index=0,fmt='pun',test_cluster_size=1000,two_dimensional=False, \
    pbc=False,supercell=[7,7,7]):
    """Create a test cluster of .xyz format centred around the chosen central atom"""
    # The purpose of this function is to generate a .xyz file cluster which contains the first test_cluster_size atoms around the central atom, given by centre_atom_index.

    # If we have a .pun file created by py-ChemShell, then it is expected that the first atom is the central atom so centre_atom_index doesn't matter
    # The below script parses out the first test_cluster_size atoms into a .xyz file
    if fmt == 'pun':
        max = 1000
        f = open('{0}cluster_test.xyz'.format(output_file_append),'w')
        f.write('{0}\n test\n'.format(max))
        with open(filename,'r') as infile:
            total_counter = 0
            i = 0
            while total_counter < max:

                a = infile.readline().split()
                if a[0] == M_element:
                    f.write('{0} {1:20.8f} {2:20.8f} {3:20.8f}\n'.format('{0}'.format(M_element),float(a[1])*Bohr,float(a[2])*Bohr, float(a[3])*Bohr))
                    total_counter += 1
                elif a[0] == O_element:
                    f.write('{0} {1:20.8f} {2:20.8f} {3:20.8f}\n'.format('O ',float(a[1])*Bohr,float(a[2])*Bohr, float(a[3])*Bohr))
                    total_counter += 1
                i += 1
        f.close()
    else:
    # If we instead have a file that is not of .pun and a traditional format accepted by ASE, then the below code applies.

        if pbc == True and two_dimensional == True:
            # This code pertains to the case where a repeating surface unit cell is provided.
            # It is expected that the vacuum region is along the z-direction.
            supercell = [supercell[0],supercell[1],1] # The supercell option should have 1 along the z-direction since it's the vacuum region.
            a = io.read(filename,format=fmt)
            sc = a*supercell # Increase the slab into a large supercell for cleaving out the test cluster
            centre_atom_sc_indices = [x for x in range(len(sc)) if sc.get_chemical_symbols()[x]==centre_atom] # Find the indices of all atoms with the same element as the centre_atom
            maxzpos = np.max(sc.get_positions()[centre_atom_sc_indices,2]) # We assume that the central atom (e.g. vacancy site) is located at the very top surface.
            # The central atom can be modified with the centre_atom_index
            top_centre_atom_indices = [x for x in centre_atom_sc_indices if abs(sc[x].position[2]- maxzpos) < 0.05]
            top_centre_atom_distance = np.linalg.norm(sc.get_positions()[top_centre_atom_indices]-[(sc.get_cell()[0][:2]+\
                sc.get_cell()[1][:2])[0]/2,(sc.get_cell()[0][:2]+sc.get_cell()[1][:2])[1]/2,maxzpos],axis=1)
            val, idx = min((val, idx) for (idx, val) in enumerate(top_centre_atom_distance))
            centre_atom_index = top_centre_atom_indices[idx] # Find index of the central atom
            sc.pbc = [False,False,False] # Remove pbc to create a cluster
            sc.translate(-sc[centre_atom_index].position) # Translate the atoms so that the central atom is at [0.0,0.0,0.0] position
            sc_sorted = sc[np.argsort(np.linalg.norm(sc.get_positions(),axis=1))[:test_cluster_size]]
            io.write('{0}cluster_test.xyz'.format(output_file_append),sc_sorted)
        elif pbc == True and two_dimensional == False: # For bulk materials
            a = io.read(filename,format=fmt)
            sc = a*supercell # Unit cell must be increased to a sufficiently large unit cell to encapsulate all RDF clusters
            centre_atom_sc_indices = [x for x in range(len(sc)) if sc.get_chemical_symbols()[x]==centre_atom]
            atom_distances_from_centre = np.linalg.norm(sc.get_positions()[centre_atom_sc_indices]-[(sc.get_cell()[0]+sc.get_cell()[1]+sc.get_cell()[2])/2],axis=1)
            val, idx = min((val, idx) for (idx, val) in enumerate(atom_distances_from_centre))
            centre_atom_index = centre_atom_sc_indices[idx]
            sc.pbc = [False,False,False] # Turn supercell into a cluster
            sc.translate(-sc[centre_atom_index].position) # Translate central atom to [0.0,0.0,0.0] position
            sc_sorted = sc[np.argsort(np.linalg.norm(sc.get_positions(),axis=1))[:test_cluster_size]] # Sort by distance from centre
            io.write('{0}cluster_test.xyz'.format(output_file_append),sc_sorted) # Write to file
        else: # For the case whereby a cluster is already created with a chosen centre atom position
            a = io.read(filename,format=fmt) # Read file
            a.translate(-a[centre_atom_index].position) # Translate central atom to [0.0,0.0,0.0] position

            # Sorts out the atoms according to distance in the next few steps
            atom_distances_from_centre = np.linalg.norm(sc.get_positions(),axis=1)
            if len(a) > test_cluster_size:
                sc_sorted = a[np.argsort(np.linalg.norm(sc.get_positions(),axis=1))[:test_cluster_size]]
            else:
                sc_sorted = a[np.argsort(np.linalg.norm(sc.get_positions(),axis=1))]
            io.write('{0}cluster_test.xyz'.format(output_file_append),a)  # Write to file

def get_atoms_list(M_list,a,O_element,MO_dist,filename):
    """Finds the coordination shell of O anions around a metal cation"""

    # MO_dist can be set manually or automatically found using in-built neighbourlist function of ASE
    O_list = []

    for i in range(len(M_list)):
        for j in range(len(a)):
            if a[j].symbol == O_element and np.linalg.norm(a[j].position - a[M_list[i]].position) < (MO_dist + 0.2):
                O_list += [j]

    total_list = list(set(M_list + O_list))
    io.write(filename,a[total_list],format='xyz') 
    return total_list

def get_rdf_cluster(filename='cluster_test.xyz',output_filename_append='',centre_index=0,MO_dist=0.5,M_element='Mg',O_element='O',num_rdf_clusters=10,prec=3):
    """Creates the series of clusters according the SKZCAM approach"""
    
    a = io.read(filename)

    # Set an unphysical default M-O bond distance so that if that happens, we will automatically find M-O distance using the ASE:
    if MO_dist < 1:
        b = natural_cutoffs(a)
        MO_dist = np.max([b[i] for i in [x for x in range(len(a)) if a.get_chemical_symbols()[x]==M_element]]) + \
            np.max([b[i] for i in [x for x in range(len(a)) if a.get_chemical_symbols()[x]==O_element]])



    M_dist_list = []
    M_indices_list = []

    # The prec command controls to how many decimal places the distance from the central atom should be accurate to.
    # If prec is too high, then it may not obtain the correct RDF shells due to inbuilt noise in the distances.
    # Alternatively, if prec is too low, it may include more than one RDF shell. prec=3 is a good default that we have used.

    precision = 10**prec
    for i in range(len(a)):
        if a[i].symbol == M_element:
            atomic_dist = np.linalg.norm(a[i].position - a[centre_index].position)
            M_dist_list += ['{0:.3f}'.format(np.floor(atomic_dist * precision)/precision)]
            M_indices_list += [i]

    # Find unique distances from the central atom corresponding to RDF peaks
    unique_M_dist_list = [float(i) for i in list(set(M_dist_list))]
    unique_M_dist_list.sort()


    # Find the number of metal ions in each RDF peak
    total_M_list = np.zeros((len(unique_M_dist_list)),dtype=int)
    num_M_list = np.zeros((len(unique_M_dist_list)),dtype=int)
    indices_M_list = []

    for i in range(len(unique_M_dist_list)):
        dummy_indices_list = []
        for j in range(len(M_dist_list)):
            if (abs(float(M_dist_list[j]) - unique_M_dist_list[i]) < 0.0015):
                num_M_list[i] += 1
                dummy_indices_list += [M_indices_list[j]]
        indices_M_list += [dummy_indices_list]
        total_M_list[i] = total_M_list[i-1] + num_M_list[i]

    # Generate the clusters corresponding to each RDF peak
    cumulative_M_list = []
    cumulative_M_indices = []
    for i in range(num_rdf_clusters):
        cumulative_M_list = cumulative_M_list + indices_M_list[i]
        cumulative_M_indices += [get_atoms_list(cumulative_M_list,a,O_element,MO_dist,\
            'Structures/{1}cluster_rdf_{0}.xyz'.format(i+1,output_filename_append))]

    # Plot out RDF plot
    fig, axs = plt.subplots(figsize=(6.69,3),dpi=300, constrained_layout=True)
    
    peak_dist = [-0.05, 0.0, 0.05]
    max_ylim = np.max(num_M_list[:num_rdf_clusters]) + 2*np.ceil(np.max(num_M_list[:num_rdf_clusters]) + 2)/10

    for i in range(num_rdf_clusters):
        dummy2 = axs.plot([x + total_M_list[i] for x in peak_dist],[0,num_M_list[i],0],linewidth=1)
        axs.text(total_M_list[i],num_M_list[i]+0.8*max_ylim/10,'{0}'.format(i+1),horizontalalignment='center',verticalalignment='center',color=dummy2[0].get_color())

        

    # axs.set_xlabel(r'Distance from central atom, $r$ (\AA{})')
    axs.set_xlabel('Quantum cluster size (# of {0} ions)'.format(M_element))

    axs.set_ylabel('RDF (# of {0} ions)'.format(M_element))
    if len(output_filename_append) > 0:
        axs.set_title('{0} system'.format(output_filename_append[:-1]))

    axs.set_xticks(total_M_list[:num_rdf_clusters])
    axs.set_xticklabels(total_M_list[:num_rdf_clusters],fontsize=7,rotation=90)


    secax = axs.secondary_xaxis('top', functions=(lambda x: 1*x, lambda x: 1*x))
    secax.set_xticks(total_M_list[:num_rdf_clusters])
    secax.set_xlabel(r'Distance from central atom (Å)')
    secax.set_xticklabels(['{0:.2f}'.format(x) for x in np.array(unique_M_dist_list)[:num_rdf_clusters]],fontsize=7,rotation=90)

    axs.yaxis.get_major_locator().set_params(integer=True)
    axs.set_xlim([0,total_M_list[num_rdf_clusters]])
    axs.set_ylim([0,max_ylim])


    print('YOU HAVE BEEN SKZCAMMED')
    plt.savefig('{0}cluster_rdf_plot.png'.format(output_filename_append),format='png')
    plt.show()
    np.save('{0}cluster_rdf_indices.npy'.format(output_filename_append),cumulative_M_indices)

def generate_SKZCAM_clusters(filename='',output_file_append='',M_element='Mg',O_element='O',centre_atom='O',centre_atom_index=0,fmt='pun',test_cluster_size=1000,two_dimensional=True, \
    pbc=False,supercell=[7,7,7],MO_dist=0.5,num_rdf_clusters=10,prec=3):
    """Single function to create RDF clusters using get_test_cluster followed by get_rdf_cluster"""

    if len(output_file_append) > 0:
        output_file_append = output_file_append + '_'
    get_test_cluster(filename,output_file_append,M_element,O_element,centre_atom,centre_atom_index,fmt,test_cluster_size,two_dimensional,pbc,supercell)

    get_rdf_cluster(filename='{0}cluster_test.xyz'.format(output_file_append),output_filename_append=output_file_append,centre_index=0,MO_dist=MO_dist,M_element=M_element\
        ,O_element=O_element,num_rdf_clusters=num_rdf_clusters,prec=prec)

