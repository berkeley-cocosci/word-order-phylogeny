#!/usr/bin/env python
import math
import os
import random
import subprocess
import time

CORES=6

def print_eta(initial, start, remaining, now):
    completed = initial - remaining
    timesofar = now - start
    pertask = timesofar / completed
    expected = remaining * pertask
    minutes = int(math.floor(expected / 60.0))
    if minutes > 60:
        hours = int(math.floor(minutes / 60.0))
        minutes = minutes - hours*60
    else:
        hours = 0
    print "Estimated time remaining: %d hours, %d minutes" % (hours, minutes)

commands = []
for treeindex in range(1,101):
    for treeclass, type in enumerate("geographic genetic feature combination".split()):
			commands.append("./inference.exe -b 500 -s 5000 -m -i %d -c %d -o results/common-q/%s/trees_%d" % (treeindex, treeclass, type, treeindex))

initial_length = len(commands)
start_time = time.time()

procs = []
for i in range(0,CORES):
    command = commands.pop()
    proc = subprocess.Popen(command.split(), stdout=open("/dev/null","w"))
    procs.append(proc)

while commands:
    time.sleep(1)
    proc = procs.pop()
    proc.poll()
    if proc.returncode != None:
        # Finished, start a new job
        print "Trees left to analyse: ", len(commands)
        print_eta(initial_length, start_time, len(commands), time.time())
        command = commands.pop()
        proc = subprocess.Popen(command.split(), stdout=open("/dev/null","w"))
    procs.insert(0, proc)

