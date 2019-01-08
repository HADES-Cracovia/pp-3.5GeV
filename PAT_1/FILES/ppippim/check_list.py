import subprocess
import time
import os


print("read list and check if files exist")

#print("type number of file")

input_name = "x23"

with open(input_name) as f:
    for line in f:
        print(line)
        if os.path.isfile(line):
            print("file :", line, "exists");
        else:
            print("file :", line ,"does not exist");
