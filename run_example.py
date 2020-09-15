from Bio import SeqIO
from frag_factory import FragFactory
from model_builder import ModelSuite
from shutil import copyfile as cp
import os, argparse
from inspect import cleandoc

#################################
# Command-line argument parsing #
#################################
parser = argparse.ArgumentParser(description='PolyPrep Module I - Polyprotein sequence fragmentat factory.')
parser.add_argument('-p','--polyprotfile',
                    help='File containing fasta sequences for target polyproteins.',
                    default=False,
                    required=True)
parser.add_argument('-i','--interf',
                    help='File containing fasta sequences for a number of positve interfaces.',
                    default=False,
                    required=True)
parser.add_argument('-l','--fraglen',
                    help='Fragment lengths to be produced.',
                    nargs='+',
                    default=False,
                    required=True, type=int)
parser.add_argument('-r','--receptor',
                    help='Receptor NAME (with no ".pdb"), for modeling',
                    default=False,
                    required=False)
parser.add_argument('-d','--dockingprep',
                    help='Docking preparation method.',
                    default=False,
                    required=False)

#################################################
# Argument parser object, produces a namespace. #
#################################################
args = parser.parse_args()                      #
#################################################

########################################
# Modelling procedures variables setup #
########################################
receptor_name = args.receptor if args.receptor else False
docking_prep_method = args.dockingprep if args.dockingprep else False

######################
# Seq IO preparation #
######################
polyprot_IO = [(fasta.id, str(fasta.seq)) for fasta in SeqIO.parse(open(args.polyprotfile),'fasta')]
polyprot_ids, polyprot_seqs = zip(*polyprot_IO)
interf_IO = [(fasta.id, str(fasta.seq)) for fasta in SeqIO.parse(open(args.interf),'fasta')]
interf_ids, interf_seqs = zip(*interf_IO)

###################################
# Fragment Factory inner workings #
###################################
factory = FragFactory()
factory.setup_factory(seq_list=polyprot_seqs, interfaces=interf_seqs, len_frags=args.fraglen)
factory.make_all_frags()

##################
# Modeling Suite #
##################

##########################
# Model Building Utility #
#################################################
# In future versions, it would be good to move  #
# all this part to modelSuite as a stati method.#
#################################################

# Check for receptor file inside current folder
if receptor_name:
    if not input(cleandoc('''WARNING! The number of fragments modelled will be huge,
                 and this could take several minutes to run.
                 To continue insert any key, to stop just press Enter.''')):
        # This is a (badly implemented) checkpoint. If the user decides that it
        # does not want to move on and create a ton of folders and wait several minutes,
        # it just interrupts the program before creating any folder.
        raise KeyboardInterrupt('Program interrupted by the user')

    
    assert receptor_name+'.pdb' in os.listdir(), 'No receptor file found on current directory'

    # if everythong goes well, create one directory for each Sequence on polyprot file
    pth = os.getcwd() # grabs cwd
    frag_files = [f for f in os.listdir() if f.startswith('frags_')] # grabs all frag_files
    run_dirs = [os.path.join(pth,'job_{}'.format(frag.replace('.fas',''))) for frag in frag_files] # makes dir names for each

    #########################################
    # **Frag File Batch execution**         #
    # Will loop through all frags and:      #
    # 1- Make a directory for it;           #
    # 2- copy receptor and frag files to    #
    # corresponding job folders             #
    # 3- Performs modelling of all frags    #
    # 4- Delete receptor file               #
    #########################################
    for d, frag in zip(run_dirs,frag_files):
        if not os.path.exists(d):
            # creates a folder that is at most 10 char long
            os.makedirs(d)
            frag_final_pth = '{}/{}'.format(d,frag)
            cp(receptor_name+'.pdb','{}/{}.pdb'.format(d,receptor_name))
            cp(frag,frag_final_pth)            
            os.chdir(d)
            m = ModelSuite(frag_file=frag_final_pth,pdb_in=receptor_name)
            _ = m.bulk()
            os.remove('{}/{}.pdb'.format(d,receptor))
            os.chdir(pth)
else:
    print("only performed frag'ing")
