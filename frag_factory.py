class FragFactory(object):
    def __init__(self,**kwargs):
        
        self.__dict__ = dict(kwargs)
        self.setup_factory(**self.__dict__)
        if self.seq_list == []:
            print("The initialized Sequence list is empty, use setup factory to update it")
        if self.len_frags == []:
            print("Initialization of Fragments Lengths not set, please refer to setup_factory.")
        
    
    def setup_factory(self, **kwargs):
        self.seq_list = self.parse_entry('seq_list', kwargs)
        self.interfaces = self.parse_entry('interfaces',kwargs)
        self.len_frags = self.parse_entry('len_frags',kwargs)
        # Forcce write)files to be always defaulted to true
        self.write_files = kwargs['write_files'] if 'write_files' in kwargs.keys() else True 

    def parse_entry(self, entry, arguments):
        
        return arguments[entry] if entry in arguments.keys() else self.seq_list if entry in vars(self).keys() else []
    
    @staticmethod        
    def make_frags(seq,len_frag, interfaces, filename='',to_file=True):
        assert len_frag >=3, "Cannot operate on fragments less than 3 residues long!"
        if len_frag > 8: print("WARNING! This script wasn't tested for fragment lengths greater than 8")
        # Asserts if the interfaces are of even length
        # otherwisethe script cannot center on a residue pair and terminates
        assert len(interfaces[0]) % 2 == 0, "You provided interfaces of even length, cannot center on a pair of residues"
        cutoff = len_frag - (1+len(interfaces[0])//2)
            
        pos_list=[]
        f = open(filename,'w')
        
        index=0
        pre_list = []
        for j in range(len(interfaces)):
            m = seq.find(interfaces[j])
            n = m+len(interfaces[j])
            frag = seq[m-cutoff:n+cutoff] if m != -1 else ''
            pre_list.append(frag)

        j=0
        final_list = []
        with open(filename, 'w') as f:
            for sequence in pre_list:
                for i in range(len(sequence)-len_frag+1):
                    if to_file:
                        header = ">pep_POS_"+str(len_frag)+'aa_'+str(j)+'-'+str(i)+'\n'
                        f.write(header)
                        f.write(sequence[i:i+len_frag]+'\n')
                    final_list.append(sequence[i:i+len_frag])
                j+=1
        return final_list

    #@staticmethod
    def make_neg(self,seq,pos,len_frag,filename='',to_file=True):
        # change it here, to get rid of the self.interfaces and provide an interface size
        cutoff = len_frag - (1+len(self.interfaces[0])//2)
        neg_list = []
        l_seq = len(seq)
        for i in range(l_seq-cutoff):
            frag = seq[i:i+len_frag]
            if frag not in pos and len(frag)==len_frag:
                neg_list.append(frag)
        if to_file:
            with open(filename,'w') as f_out:
                for i in range(len(neg_list)):
                    header = ">pep_NEG_"+str(len_frag)+'aa_'+str(i)+'\n'
                    f_out.write(header)
                    f_out.write(neg_list[i]+'\n')
        return neg_list
    
    def make_all_frags(self):
        # Execution loop!
        # Assertions regarding input data
        assert len(self.seq_list) > 0, "Cannot operate on an empty sequence list."
        assert self.interfaces != [], "No interfaces set. The program will operate but on non-discriminatory fragmentation."
        self.frags = {}
        i=1
        for seq_i,seq in enumerate(self.seq_list):
            for l in self.len_frags:                    
                pos_frags = self.make_frags(seq,l, self.interfaces, "Length_{}_POS_seq_{}.fas".format(l,seq_i))
                neg_frags = self.make_neg(seq,pos_frags, l,"Length_{}_NEG_seq_{}.fas".format(l,seq_i))
                self.frags['seq_{}_{}aa'.format(i,l)] = [pos_frags, neg_frags]
