import numpy as np
import os
import sys

if(len(sys.argv)<2):
   print("Syntax: check_nev.py file")



for i in range(1,len(sys.argv)):
    infile=sys.argv[i]
  
    if(not (os.path.exists(infile))):
      print("Error, "+infile+" does not exist!\n")
      exit(1)
    else:
      print("Working on file: "+infile)
    with open(infile,"rb") as fin:
      nevents=np.fromfile(fin,dtype=np.int64,count=1)[0]
      if(i==1):
        nref=nevents
        print("Reference number of events is: "+str(nevents))
      else:
        if(nevents!=nref):
          print("Error, "+infile+" has a different number of events: "+str(nevents)+"\n")
          exit(1)
print("Check completed, all fine!") 
