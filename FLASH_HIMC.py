import re, copy, gzip, os, sys, datetime,subprocess
from subprocess import PIPE 
mate1=sys.argv[1]
mate2=sys.argv[2]
path_ = os.getcwd()
pref_= sys.argv[3]

flash_call = subprocess.Popen('flash '+mate1+' '+mate2 +' -O -o '+pref_ + ' -d '+ path_ +' -m 100 -M 250', shell = True, stdout=PIPE, stderr=PIPE)
