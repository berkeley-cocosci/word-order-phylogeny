#!/usr/bin/env python
import getopt
import math
import os
import random
import subprocess
import sys
import time

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

def main():

    cores = 1
    burnin = 5000
    samples = 1000
    lag = 100
    treedepth = 100
    opts, args = getopt.gnu_getopt(sys.argv[1:], "c:b:s:l:t:")
    for o, a in opts:
        if o == "-c":
            cores = int(a)
        elif o == "-b":
            burnin = int(a)
        elif o == "-s":
            samples = int(a)
        elif o == "-l":
            lag = int(a)
        elif o == "-t":
            treedepth = int(a)

    commands = []
    for method in range(0, 4):
        # COMMON Q
        # Unsplit, unshuffled
        commands.append("./inference.exe -m -b %d -s %d -l %d -t %d -c %d" % (burnin, samples, lag, treedepth, method))
        # Unsplit, shuffled
        commands.append("./inference.exe -m -b %d -s %d -l %d -t %d -c %d -S" % (burnin, samples, lag, treedepth, method))
        # Split, unshuffled
        commands.append("./inference.exe -m -b %d -s %d -l %d -t %d -c %d -x" % (burnin, samples, lag, treedepth, method))
        # INDIVIDUAL Q
        # Unsplit, unshuffled
        commands.append("./inference.exe -b %d -s %d -l %d -t %d -c %d" % (burnin, samples, lag, treedepth, method))

    initial_length = len(commands)
    start_time = time.time()

    procs = []

    # Start initial procs, one per core
    for i in range(0,cores):
        command = commands.pop()
        print "Starting an initial job..."
        proc = subprocess.Popen(command.split(), stdout=open("/dev/null","w"))
        procs.append(proc)

    # Poll procs to replace when finished
    while procs:
        time.sleep(1)
        proc = procs.pop()
        proc.poll()
        if proc.returncode != None:
            # Finished, start a new job
            print "Trees left to analyse: ", len(commands)
            print_eta(initial_length, start_time, len(commands), time.time())
            if commands:
                command = commands.pop()
                proc = subprocess.Popen(command.split(), stdout=open("/dev/null","w"))
                procs.insert(0, proc)
        else:
            procs.insert(0, proc)

if __name__ == "__main__":
    main()
