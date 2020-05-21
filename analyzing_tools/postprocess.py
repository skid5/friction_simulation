import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from scipy import optimize
from numpy import mean, linspace
import sys
import os

# Function declarations.

def print_help():
    # TODO: Kirjoita apu.
    print "Usage: python " + sys.argv[0] + " [-a] [-v] [-n=number_of_simulations] [-p=file_prefix]" 
    print "-a\t--auto\t\t\tuse automatic mode. Deduce file names and amount of files."
    print "-p=\t--prefix=\t\tuse a custom file prefix. Default is avr_force."
    print "-n=\t--nsimulations=\t\tuse a custom number of simulations. Default is 4. Is omitted if automatic mode is used."
    print "\t\t\t\t Using the linspace option will also override this option."
    print "-v\t--verbose\t\tprint lots of information."
    print "-f=\t\t\t\tset a custom filename for the force list file or create the force list using linspace."
    print "\t\t\t\tUsage -f=filename.txt of -f=low,high,steps. Default action: read data from force_list.txt"
    print "-z\t\t\t\tif set the program will remove zeros from the data. Default action is to not remove zeros."

# TODO: Error checking.
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

# Global variable declarations.

auto = False
have_prefix = False
verbose = False
force_file = True
nozeros = False
n_simulations = 10
prefix = "avr_force"
averages_fn = "averages.txt"
output_fn = "plots.png"
force_list_fn = "force_list.txt"
linspace_min = linspace_max = linspace_steps = -1
force_list = []

for s in sys.argv[1:]:
    # Choose files automatically. Must have prefix, ignore n_simulations.
    if s == "-a" or s == "--auto":
        auto = True
    elif s == "-v" or s == "--verbose":
        verbose = True
    elif s == "-z":
		nozeros = True
    elif s.startswith("-n=") or s.startswith("--nsimulations="):
        n_simulations = int(s.split("=")[-1])
    elif s.startswith("-p=") or s.startswith("--fileprefix="):
        prefix = s.split("=")[-1]
        have_prefix = True
    elif s.startswith("-f="):
		s = s.replace("=", " ").replace(","," ").split()
		if len(s) == 2:
			# Just change default name to a custom name.
			force_list_fn = s[1]
			print 
		elif len(s) == 3:
			# Tell the program that we will not read force data from a file.
			force_file = False
			# Set parameters for linspace, but don't create linspace yet.
			linspace_min = float(s[1])
			linspace_max = float(s[2])
			
		elif len(s) == 4:
			# TODO: Maybe do error checks here for negative values.
			# Same as the previous branch, but also sets linspace steps.
			force_file = False
			linspace_min = float(s[1])
			linspace_max = float(s[2])
			linspace_steps = n_simulations = int(s[3])
    else:
        print_help()
        sys.exit()

if force_file:
    force_list = gather_data(force_list_fn)
    force_list = [a[0] for a in force_list]

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
    print "Reading forces from " + force_list_fn + "."
    if nozeros:
        print "Zeros will be removed in the data."
    else:
        print "Zeros will be left in the data."

file_list = []
if auto:
    ls = os.listdir(".")
    file_list = [fn for fn in ls if fn.startswith(prefix)]
	
	# If linspace steps wasn't given as a parameter, it is the same as the number of files.
    if linspace_steps == -1:
        linspace_steps = len(file_list)
	# If the forces were not read from file only then create them using linspace.
	if not force_file:
	    force_list = linspace(linspace_min, linspace_max, linspace_steps).tolist()
	    if verbose:
		    print "Linspace is using parameters low = " + str(linspace_min) + ", high = " + str(linspace_max) + ", steps = " + str(linspace_steps)
		
else:
    for i in range(n_simulations):
        file_list.append(prefix + str(i) + ".txt")


if verbose:
    print "Processing the following files:",
    for fn in file_list:
        print fn,
    print
    print "Writing output to " + averages_fn + "."


averages = []

# If force list is still empty at this point, print a warning.
if len(force_list) == 0:
	print "Warning: Force list is empty!"
	print "Using a dirty hack."
	force_list = range(len(file_list))

for fn in file_list:
    res = gather_data(fn)
    x, y, z = zip(*res)
	# Remove zeros.
    if nozeros:
        x = [i for i in x if i != 0] 
        y = [i for i in y if i != 0] 
    averages.append( (mean(x), mean(y)) )

write_data(averages, averages_fn)

# Make data ready for plotting.
avgs_x, avgs_y = zip(*averages)

# If the data vectors are of different length at this point, the program must be halted.
if len(avgs_x) != len(avgs_y) or len(avgs_y) != len(force_list):
	print "Warning: data vectors are of different length."
	print "Aborting program."
	sys.exit()

# Plotting!
# What do we want?
# We want 2 graphs in one image. (x- and y-components) [pyplot.plot(t, f1(t), t, f2(t), ...)]
# We want the average friction as a function of velocity or force. [pyplot.plot(force, avg_friction_x, force, ...)]
# We want to plot the data points as points(crosses). [pyplot.plot(x,y,"rx") # (red crosses)]
# We want to fit a curve to data.
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

# Fitting!
# https://towardsdatascience.com/basic-curve-fitting-of-scientific-data-with-python-9592244a2509
# https://scipy-lectures.org/intro/scipy/auto_examples/plot_curve_fit.html
def fitting_function(x, a, b):
	# Line as a placeholder
	return a * x + b

# x_params, x_params_covariance = optimize.curve_fit(fitting_function, force_list, avgs_x)
# y_params, y_params_covariance = optimize.curve_fit(fitting_function, force_list, avgs_y)
