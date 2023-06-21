import os
import sys

if(len(sys.argv) > 1):
   if(sys.argv[1] == '-e'):
      rm_extra_files = False
   else:
      rm_extra_files = True
else:
   rm_extra_files = True

cwd = os.getcwd()

for root, dirs, files in os.walk(".", topdown=False):
   for file in files:
      if file.endswith(".tex"):
         path = os.path.join(root, file)
         name = os.path.splitext(file)[0]

         print("The Current working directory is: {0}".format(os.getcwd()))
         print(root)
         print(file)
         print(path)

         os.chdir(root)

         # Only convert to pdf if the latex file does not already have a pdf version
         if(not os.path.isfile(name + '.pdf')): 
            os.system('pdflatex ' + file)

         # remove extra file unless specified by user with '-e' (extra)
         if(rm_extra_files):
            files_to_delete = name + '.aux ' + name + '.log ' + name + '.nav ' + name + '.out ' + name + '.snm ' + name + '.toc ' + name + '.vrb'
            os.system('rm ' + files_to_delete)
         print("The Current working directory now is: {0}".format(os.getcwd()))

         os.chdir(cwd)
         print("The Current working directory now is: {0}".format(os.getcwd()))
