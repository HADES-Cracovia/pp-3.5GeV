import subprocess
import time
import os


print("Copy and hadd files from shorter lists")
print("type number of file")


value = input()

input_name="x{:02d}".format(value)
output_name="hadron{:02d}.root".format(value)
print("file name: {}".format(input_name))
argument_list=""

with open(input_name) as f:
    content = f.readlines()
    content = [x.strip() for x in content]
    
argument_list=argument_list+" ".join(content)
bashCommand="hadd -f "+output_name+" "+argument_list
#print(bashCommand)

os.system(bashCommand)
