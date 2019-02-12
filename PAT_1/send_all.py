#directory="/lustre/nyx/hades/user/iciepal/Lambda1520_ic/" #oryginal data directory
#fname="list_allfiles" #file with list of expected files
import subprocess
print("run all jobs")
import time
import os

#with open(fname) as f:
 #   content = f.readlines()
#content = [x.strip() for x in content]

os.system("scancel -u knowakow")
print("scancel -u knowakow")

print("./sendScript_01.sh")
os.system("./sendScript_01.sh")

time.sleep(10)

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)


print("./sendScript_02.sh")
os.system("./sendScript_02.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)


print("./sendScript_03.sh")
os.system("./sendScript_03.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)


print("./sendScript_04.sh")
os.system("./sendScript_04.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)

print("./sendScript_05.sh")
os.system("./sendScript_05.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)
     

print("./sendScript_06.sh")
os.system("./sendScript_06.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)

print("./sendScript_07.sh")
os.system("./sendScript_07.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)

print("./sendScript_08sh")
os.system("./sendScript_08.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)

print("./sendScript_09.sh")
os.system("./sendScript_09.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)

print("./sendScript_10.sh")
os.system("./sendScript_10.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)

print("./sendScript_11.sh")
os.system("./sendScript_11.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)

print("./sendScript_00.sh")
os.system("./sendScript_00.sh")

#for k in content:   #take every name from vector content
#   print('try to run files from list{}'.format(k))
while(1 != int(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,))):
    #print('waiting to finish list: {}'.format(k))
    print('still {} jobs to the end'.format(subprocess.check_output('squeue -u knowakow | wc -l',shell=True,)))
    time.sleep(30)
