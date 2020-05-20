from numpy import mean
import sys
import os
auto = False
have_prefix = False
verbose = False
n_simulations = 4
prefix = "avr_force"
suffix = ".txt"
avg_x_fn = "average_x.txt"
avg_y_fn = "average_y.txt"
def print_help():
    print("TODO: Kirjoita apu.")

for s in sys.argv[1:]:
    # Choose files automatically. Must have prefix, ignore n_simulations.
    if s == "-a" or s == "--auto":
        auto = True
    elif s == "-v" or s == "--verbose":
        verbose = True
    elif s.startswith("-n=") or s.startswith("--nsimulations="):
        n_simulations = int(s.split("=")[-1])
    elif s.startswith("-p=") or s.startswith("--fileprefix="):
        prefix = s.split("=")[-1]
        have_prefix = True
    else:
        print_help()
        sys.exit()

if verbose:
    if have_prefix:
        print "Using prefix = " + prefix + "."
    else:
        print "Using default prefix."
    print "Using n_simulations = " + str(n_simulations) + "."
    if auto:
        print "Using automatic mode."
    else:
        print "Using manual mode."


file_list = []
if auto:
    ls = os.listdir(".")
    file_list = [fn for fn in ls if fn.startswith(prefix)]
else:
    for i in xrange(n_simulations):
        file_list.append(prefix + str(i) + suffix)

if verbose:
    print "Processing the following files:",
    for fn in file_list:
        print fn,
    print
    print "Writing output to " + avg_x_fn + " and " + avg_y_fn + "."

for fn in file_list:
    fp = open(fn, "r")
    x = []
    y = []
    while True:
        line = fp.readline()
        if line == '':
            break
        line = line.split()
        x.append(float(line[0]))
        y.append(float(line[1]))

    fp.close()

    fp_x = open(avg_x_fn, "a")
    fp_y = open(avg_y_fn, "a")

    average_x = mean(x)
    average_y = mean(y)

    fp_x.write(str(average_x) + '\n')
    fp_y.write(str(average_y) + '\n')

    fp_x.close()
    fp_y.close()
