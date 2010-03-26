#---------------------------------------
#
# Author     Perry Evans
#            evansjp@mail.med.upenn.edu
# 2008
#---------------------------------------
"""
This has methods common to getting ELM HTML.
"""
def checkForErrors(file_to_check):
    """ Return true if the ELM tool had an error
        or if it crashed (indicted by no summary).
    @param file_to_check: ELM HTML file
    @return True if there was a problem; False otherwise
    """

    have_error = False
    found_summary_line = False
    try:
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
    return (have_error or not found_summary_line) 
