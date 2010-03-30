#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
"""
This has methods common to getting ELM HTML.
"""
import os

def checkForErrors(file_to_check, elm_html_dir, protein):
    """ Return true if the ELM tool had an error
        or if it crashed (indicted by no summary).
    @param file_to_check: ELM HTML file
    @return True if there was a problem; False otherwise
    """

    have_error = False
    found_summary_line = False
    tidy_file = elm_html_dir + protein + '.elm_tidy.html'

    try:
        # is the file there?
        elm_file = open(file_to_check)
        elm_file.close()
        # if the file is there, try to clean it
        os.system('tidy -f errs -o '
                  + tidy_file + ' '
                  + file_to_check)
        # if tidy is empty it will not find the summary
        f_open = open(file_to_check)
        for line in f_open:
            if line.find('error message') != -1:
                have_error = True
                break
            elif line.find('Summary of features') != -1:
                found_summary_line = True
        f_open.close()            
    except:
        have_error = True
    
    if have_error or not found_summary_line:
        #os.system('rm -f '
        #          + file_to_check)
        os.system('rm -f '
                  + elm_html_dir + protein
                  + '.elm_tidy.html')
        return True
    else:
        return False

