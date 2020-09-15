from Bio import SeqIO
from modeller import *
import subprocess
from os import path, makedirs,chdir,getcwd,listdir,remove
from modeller.automodel import *
import traceback
#########################################################################
#  This script consists of a reimplementation of                        #
#  Modeller's tutorial script. Since it uses default                    #
#  parameters, and those worked for most protocols I ran.               #
#  Have in mind that this script was written around 3 years             #
#  ago and consists of one of my first "useful" scripts.                #
#  In June 2020 I edited this code but did not change much              #
#  of it's structures. All the rights for the modelling parts           #
#  are reserved for Modeller's developers, and so the license for the   #
#  current implementation follows the same rules. Feel free to          #
#  change all the parameters adopted below.                             #
#  Just for the record, the changes made here regard running the        #
#  original scripts on a batch fashion.                                 #
#########################################################################

log.none()
class ModelSuite:
    def __init__(self, **kwargs):        
        self.__dict__ = dict(kwargs)
        self.setup_msuite(**self.__dict__)
        self.frag_list=[]
        self.absolute_pth= getcwd()
        self.env = environ()
        self.env.libs.topology.read('$(LIB)/top_heav.lib')
        self.env.io.atom_files_directory = self.absolute_pth+'/'
        self.cwd = getcwd()
        self.receptor=self.absolute_pth+'/'+self.pdb_in+'.pdb' #

    def setup_msuite(self, **kwargs):
        self.frag_file = self.parse_entry('frag_file', kwargs)
        self.pdb_in = self.parse_entry('pdb_in',kwargs)
        self.prep_4_docking = self.parse_entry('prep_4_docking',kwargs)

    def parse_prep_parm(self, arguments):
        if 'prep_4_docking' in arguments.keys():
            docking_procedure = arguments['prep_4_docking'].upper()
            assert docking_procedure == 'AD4' or docking_procedure == "VINA", "Wrong Docking type procedure, please determine wether Autodock4 ou Vina should be adopted"
            if arguments[prep_4_docking].upper() == 'AD4':
                self.prep_4_docking = ['--AddPolarH', '--partialcharge gasteiger']                                            
        else:
            self.prep_4_docking = 'VINA'
            #defaults to Vina prep
        
            
    def parse_entry(self, entry, arguments):        
        return arguments[entry] if entry in arguments.keys() else self.seq_list if entry in vars(self).keys() else ''

    def getSeqs(self):
        
        seq_buff = SeqIO.parse(open(self.frag_file),'fasta')
        for fasta in seq_buff:
            name, sequence = fasta.id, str(fasta.seq)
            self.frag_list.append((name,sequence))
        
    
    def create_ali(self,frag=('',''),pth='./'):
        ali_file = path.join(pth,'alinhamento.ali')

        with open(ali_file, 'w') as outfile:
            outfile.write('>P1;'+self.pdb_in+'.pdb\n')
            outfile.write('structureX:{}.pdb: FIRST:A : LAST :A ::::\n'.format(self.pdb_in))
            outfile.write('*\n')
            outfile.write('>P1;'+frag[0]+'\n')
            outfile.write('sequence:::::::::\n')
            outfile.write(frag[1]+'*')

    def align_core(self,seq):
    #if self.cwd==self.path:
        log.none()
        alvo = seq[0]
        template = self.pdb_in+'.pdb'
        aln = alignment(self.env)
        self.env.io.atom_files_directory = ['./',self.cwd,self.absolute_pth]
        aln_pth = path.join(self.cwd,"alinhamento.ali")
        aln.append(file=aln_pth, align_codes=template)
        aln_block = len(aln)
        aln.append(file=aln_pth, align_codes=alvo)
        aln.align2d(overhang=0, gap_penalties_1d=(-100, 0),
        gap_penalties_2d=(3.5, 3.5, 3.5, 0.2, 4.0, 6.5, 2.0, 0., 0.),
        align_block=aln_block)
        aln_file = path.join(self.cwd,'align2d_{}_'.format(alvo))
        aln.write(file=aln_file+'.ali', alignment_format='PIR')
        aln.write(file=aln_file+'.pap', alignment_format='PAP',
        alignment_features='INDICES HELIX BETA STRAIGHTNESS ' + \
        'ACCESSIBILITY CONSERVATION')
        aln.check()
        aln = alignment(self.env)
        aln.append(file=aln_file+'.ali', align_codes=(template, alvo),
        alignment_format='PIR', remove_gaps=True)
        return aln_file+'.ali'

    def modeling(self,aln_file,frag):
        chdir(path.join(self.cwd,'MODEL/'))
        self.env.io.atom_files_directory = ['./',self.cwd,self.absolute_pth, path.join(self.cwd,'MODEL/')]
        engine = automodel(self.env, alnfile = path.join(self.cwd,aln_file),
                      knowns = self.pdb_in+'.pdb', # codes of the templates
                        sequence = frag[0],
                      assess_methods=(assess.DOPE, assess.GA341))
        engine.starting_model= 1
        engine.ending_model = 1
        # start the engine
        engine.make()
        # Get a list of all successfully built models from the model engine.outputs
        models = list(filter(lambda x: x['failure'] is None, engine.outputs))
        models = models[0]
        
        filelist = [ f for f in listdir(".") if f.endswith(".pdb") ]
        for f in filelist:
            if f != models['name']:
                remove(f)
        return (models['name'], models['DOPE score'])

    def bulk(self):
        try:
            self.getSeqs()
        except Exception as er:
            print(er)
            traceback.print_exc()
            print("Problems encountered at sequence parsing step, please review the provided file paths")
        models=[]
        receptor=self.receptor
        for i in range(len(self.frag_list)):
            pth = path.join(self.absolute_pth,'Seq_{}/'.format(self.frag_list[i][0]))
            if not path.exists(pth):
                makedirs(path.join(pth,'MODEL'))
            self.create_ali(self.frag_list[i], pth)
            self.cwd = pth
            try:
                aln_file = self.align_core(self.frag_list[i])
                models.append(self.modeling(aln_file,self.frag_list[i]))
            except Exception as er:
                print("Alignment and modeling steps interrupted")
                print(er)
                traceback.print_exc()
        
            '''
            Model preparation for Docking follows a simple rule for Autodock Vina:
            PolarH and file conversion to pdbqt (the first is not mandatory, but we
            used to run this is cript for Autodock4)
            '''
            
            cmd_line = 'babel -ipdb {0}.pdb -opdbqt {0}.pdbqt'.format(models[-1][0].replace('.pdb',''))
            if type(self.prep_4_docking) == list:
                parm = self.prep_4_docking
                cmd_line +='{} {}'.format(*parm)
            
            _ = subprocess.getstatusoutput(cmd_line)
            print(_)
            chdir(self.absolute_pth)
            

        return models



