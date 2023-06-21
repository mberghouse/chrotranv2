"""
Calvin Madsen
6/5/2023

Script should compile all latex files into pdf files if pdflatex is installed
Automatically removes extra 'log' files that come with running pdflatex on
a file. This can be disabled by adding '-e' at the end of running the script.
"""

import os
import sys


def pdflatexConversion(rm_extra_files):

    cwd = os.getcwd()

    for root, dirs, files in os.walk(".", topdown=False):
        for file in files:
            if file.endswith(".tex"):
                name = os.path.splitext(file)[0]

                os.chdir(root)

                # Only convert to pdf if the latex file does not already have a pdf version
                if not os.path.isfile(name + '.pdf'): 
                    os.system('pdflatex ' + file)

                # remove extra file unless specified by user with '-e' (extra)
                if rm_extra_files:
                    files_to_delete = name + '.aux ' + name + '.log ' + name + '.nav ' + name + '.out ' \
                                    + name + '.snm ' + name + '.toc ' + name + '.vrb'
                    os.system('rm ' + files_to_delete + ' 2>&-')

                os.chdir(cwd)


if __name__ == "__main__":
    if(len(sys.argv) > 1):
        if(sys.argv[1] == '-e'):
            rm_extra_files = False
        else:
            rm_extra_files = True
    else:
        rm_extra_files = True

    pdflatexConversion(rm_extra_files)
