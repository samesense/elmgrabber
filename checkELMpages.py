#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
"""
This is used while getting ELM HTML pages.
Some pages come back with errors.  This finds
those pages and prints a list of proteins that
need to be rescanned with the website.
"""
import utils_scripting, utils_graph, elm_tools
import sys

req_args = ['list of protiens', 'elm html dir']
examples = ['../../Data/Network/Human/HPRD/hprd.intr.ls',
            '../../Data/ELM/Human/HTML/']
utils_scripting.checkStart(sys.argv, req_args, examples, 2, True)
protein_ls_file = sys.argv[1]
elm_html_dir = sys.argv[2]

proteins = utils_graph.getNodes(protein_ls_file)

redo_proteins = {}
for protein in proteins.keys():
    if elm_tools.checkForErrors(elm_html_dir + protein + '.elm_tidy.html'):
        redo_proteins[protein] = True

for protein in redo_proteins.keys():
    print protein
