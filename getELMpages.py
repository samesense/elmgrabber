#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
"""
This is a wrapper for getELMpages_runner.
The program takes a list of proteins to 
annotate, a FASTA file, and an output directory.
It calls ckechELMpages to find which proteins
need annotations, and splits the protein list
into 30 pieces and calls the runner on the pieces.
"""
import utils_motif, utils_fasta, utils_scripting, utils_graph
import sys, os
from threading import Thread

req_args = ['list of genes',
            'fasta file',
            'output dir']
examples = ['../../Data/Network/Human/HPRD/hprd.intr.ls',
            '../../Data/FASTA/Human/hprd.intr.fasta',
            '../../Data/ELM/Human/HTML/']

utils_scripting.checkStart(sys.argv, req_args, examples, 3, True)

fasta_file = sys.argv[2]
html_dump_dir = sys.argv[3]

class elm_thread(Thread):
    def __init__(self, redo_file):
        Thread.__init__(self)
        self.file = redo_file

    def run(self):
        os.system('python getELMpages_runner.py '
                  + self.file + ' '
                  + fasta_file + ' '
                  + html_dump_dir)
        
def split_seq(seq_dict, size):
    """ Split a {} into (size) equal pieces and return [] of pieces. """
    
    seq = []
    for item in seq_dict.keys():
        seq.append(item)
    newseq = []
    splitsize = 1.0/size*len(seq)
    for i in range(size):
        newseq.append(seq[int(round(i*splitsize)):int(round((i+1)*splitsize))]) 
    return newseq

def ls2dict(ls):
    """ Convert ls to {}. """
    
    new_dict = {}
    for item in ls:
        new_dict[item] = True
    return new_dict

gene_ls_file = sys.argv[1]

os.system('python checkELMpages.py ' + gene_ls_file
          + ' ' + html_dump_dir + ' > redo.protein.ls')

protein_ls = utils_graph.getNodes('redo.protein.ls')
protein_pieces_ls = split_seq(protein_ls, 30)
process_ls = []
for index in xrange(len(protein_pieces_ls)):
    protein_dict = ls2dict(protein_pieces_ls[index])
    redo_file = 'redo.protein.ls_' + str(index)
    utils_graph.dumpNodes(redo_file, protein_dict)
    current = elm_thread(redo_file)
    process_ls.append(current)
    current.start()
    #process_ls.append(subprocess.Popen('python getELMpages_runner.py redo.protein.ls_'
    #                                   + str(index) + ' ' + fasta_file + ' '
    #                                   + html_dump_dir, shell=True))
for elm_process in process_ls:
    elm_process.join()
