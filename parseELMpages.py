#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
"""
This takes tidy resutls of ELM HTML pages
and prints annotation results for protiens.
It requires a list of proteins, the directory
with the tidy results, and an output file.
The output format is:
protein start stop ELMname seq ELM.
"""
import utils_motif, utils_scripting, utils_graph, elm_tools
import sys, os

req_args = ['list of genes to annotate',           
            'elm html dir',
            'elm output file']
examples = ['../../Data/Network/Human/HPRD/hprd.intr.ls',
            '../../Data/ELM/Human/HTML/',
            '../../Data/ELM/Human/human.elm']
utils_scripting.checkStart(sys.argv, req_args, examples, 3, True)
genes_to_annotate_file = sys.argv[1]
html_page_dir = sys.argv[2]
elm_output_file = sys.argv[3]

genes = utils_graph.getNodes(genes_to_annotate_file)

elm_output = open(elm_output_file, 'w')
for gene in genes.keys(): 
    #print 'GENE\t' + gene
    file_name = html_page_dir + gene + '.elm_tidy.html'
    have_file = True
    try:
        f_open = open(file_name)
        f_open.close()
    except:
        have_file = False
    #print gene
    if have_file and not elm_tools.checkForErrors(file_name):
        #print 'yes'
        rows = utils_motif.parseELMpage_v2(gene, file_name)
        for row in rows:
            line = ''
            for item in row:
                line = line + item + '\t'
            elm_output.write(line.strip('\t') + '\n')            
elm_output.close()
