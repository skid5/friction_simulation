from numpy import mean
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
import sys
import os
auto = False
have_prefix = False
verbose = False
n_simulations = 4
prefix = "avr_force"
suffix = ".txt"
averages_fn = "averages.txt"
output_fn = "plots.png"
# TODO: Read force/velocity list from file or command line.
force_list = range(5)
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
    # TODO: Poista purkkaratkaisu
    force_list = range(len(file_list))
else:
    for i in xrange(n_simulations):
        file_list.append(prefix + str(i) + suffix)

if verbose:
    print "Processing the following files:",
    for fn in file_list:
        print fn,
    print
    print "Writing output to " + averages_fn + "."


def gather_data(fn):
    fp = open(fn, "r")
    result = []
    while True:
        line = fp.readline()
        if line == '':
            break
        line = line.split()
        floats = [float(f) for f in line]
        result.append(floats)
    fp.close()
    return result

def write_data(data, fn, sep = "\t"):
    fp = open(fn, "w")
    for datavec in data:
        for d in datavec:
            fp.write(str(d) + sep)
        fp.write("\n")
    fp.close()
    
averages = []
for fn in file_list:
    res = gather_data(fn)
    x, y, z = zip(*res)
    averages.append( (mean(x), mean(y)) )

write_data(averages, averages_fn)

# Make data ready for plotting.
avgs_x, avgs_y = zip(*averages)

# Plotting!
# What do we want?
# We want 2 graphs in one image. (x- and y-components) [pyplot.plot(t, f1(t), t, f2(t), ...)]
# We want the average friction as a function of velocity or force. [pyplot.plot(force, avg_friction_x, force, ...)]
# We want to plot the data points as points(crosses). [pyplot.plot(x,y,"rx") # (red crosses)]
# We might want to fit a curve to data.
# We want to write the graphs to files. [pyplot.savefig(...) # (supports png and pdf at least.)]
# https://matplotlib.org/tutorials/introductory/pyplot.html

# These settings can be modified.
plt.plot(force_list, avgs_x, "bx", label = "x-component")
plt.plot(force_list, avgs_y, "rx", label = "y-component")
plt.xlabel("Normal force")
plt.ylabel("Friction force")
plt.title("Friction force as a function of normal force")
#plt.axis([xmin, xmax, ymin, ymax])
plt.grid(True)
plt.legend(loc = "upper left")
plt.savefig(output_fn)

