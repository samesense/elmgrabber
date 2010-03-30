#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
"""
This interacts with the ELM website to get
ELM HTML pages.  The ELM HTML pages are not
'proper' HTML, so tidy is called to fix them.
The program takes a list of proteins to 
annotate, a FASTA file, and an output directory.
"""
import utils_motif, utils_fasta, utils_scripting, utils_graph, elm_tools
import sys, os, time

req_args = ['list of genes',
            'fasta file',
            'output dir']
examples = ['../../Data/Network/Human/HPRD/hprd.intr.ls',
            '../../Data/FASTA/Human/hprd.intr.fasta',
            '../../Data/ELM/Human/HTML/']

utils_scripting.checkStart(sys.argv, req_args, examples, 3, True)
html_dump_dir = sys.argv[3]

genes = utils_graph.getNodes(sys.argv[1])
fasta = utils_fasta.loadFASTA(sys.argv[2])

for gene in genes.keys():
    seq = fasta[gene]
    then_time = time.clock()
    current_time = time.clock()
    try:
        ofile = html_dump_dir + gene + '.elm.html'
        utils_motif.getELMpage(gene, seq, ofile)
        elm_tools.checkForErrors(ofile, elm_dump_dir, gene)
    except:
        pass
    
