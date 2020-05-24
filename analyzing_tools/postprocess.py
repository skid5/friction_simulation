from __future__ import print_function
import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from scipy import optimize
from numpy import mean, linspace
from input_output import *
import sys
import os

# Function declarations.

def print_help():
    # TODO: Kirjoita apu.
    print("Usage: python " + sys.argv[0] + " [-a] [-v] [-n=number_of_simulations] [-p=file_prefix]")
    print("-a\t--auto\t\t\tuse automatic mode. Deduce file names and amount of files.")
    print("-p=\t--prefix=\t\tuse a custom file prefix. Default is avr_force.")
    print("-n=\t--nsimulations=\t\tuse a custom number of simulations. Default is 4. Is omitted if automatic mode is used.")
    print("\t\t\t\t Using the linspace option will also override this option.")
    print("-v\t--verbose\t\tprint lots of information.")
    print("-f=\t\t\t\tset a custom filename for the force list file or create the force list using linspace.")
    print("\t\t\t\tUsage -f=filename.txt of -f=low,high,steps. Default action: read data from force_list.txt")
    print("-z\t\t\t\tif set the program will remove zeros from the data. Default action is to not remove zeros.")
    print("-t\t\t\t\tmake velocity/force/energy vs time plot if set. force_list should now contain times. Default dt = 4 and SIMTIME = 1e5.")
    print("-i\t\t\t\tset plotting index. Default 0. 0 = Epot, 1 = Ekin, 2 = Etot, 3 = temp. 0 is also x-component of friction.")
    print("\t\t\t\tThis value basically refers to the columns of the data file.")

# Global variable declarations.

auto = False
have_prefix = False
verbose = False
force_file = True
nozeros = False
timeplot = False
n_simulations = 10
prefix = "avr_force"
averages_fn = "averages.txt"
force_list_fn = "force_list.txt"
linspace_min = linspace_max = linspace_steps = -1
force_list = []
# 0 = pot, 1 = kin, 2 = tot, 3 = temp
# 0 = x-component of friction also.
plotting_index = 0

SIMTIME = 1e5
dt = 4

for s in sys.argv[1:]:
    # Choose files automatically. Must have prefix, ignore n_simulations.
    if s == "-a" or s == "--auto":
        auto = True
    elif s == "-v" or s == "--verbose":
        verbose = True
    elif s == "-z":
		nozeros = True
    elif s == "-t":
		timeplot = True
    elif s.startswith("-n=") or s.startswith("--nsimulations="):
        n_simulations = int(s.split("=")[-1])
    elif s.startswith("-i="):
        plotting_index = int(s.split("=")[-1])
    elif s.startswith("-p=") or s.startswith("--fileprefix="):
        prefix = s.split("=")[-1]
        have_prefix = True
    elif s.startswith("-f="):
		s = s.replace("=", " ").replace(","," ").split()
		if len(s) == 2:
			# Just change default name to a custom name.
			force_list_fn = s[1]
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

# Make the file list.
file_list = []
if auto:
	ls = os.listdir(".")
	file_list = [fn for fn in ls if fn.startswith(prefix)]

	# If linspace steps wasn't given as a parameter, it is the same as the number of files.
	if linspace_steps == -1:
		linspace_steps = len(file_list)

else:
    for i in range(n_simulations):
        file_list.append(prefix + str(i) + ".txt")

if verbose:
    if have_prefix:
        print("Using prefix = ", prefix, ".", sep = "")
    else:
        print("Using default prefix.")
    if auto:
        print("Using automatic mode.")
    else:
        print("Using manual mode.")
    print("plotting_index = ", plotting_index)
    if force_list:
		print("Reading forces from ", force_list_fn, ".", sep ="")
    if nozeros:
        print("Zeros will be removed from the data.")
    else:
        print("Zeros will be left in the data.")
    print("Processing the following files: ", end = "")
    for fn in file_list:
        print(fn, " ", end = "")
    print()
    print("Writing output to ", averages_fn, ".", sep = "")

averages = []
if timeplot:
	for i in range(len(file_list)):
		fn = file_list[i]
		res = gather_data(fn)
		res = zip(*res)
		x = res[plotting_index]		
		# TODO: DRY-periaate.
		# TODO: How to calculate time_list using SIMTIME and dt? More precisely how we went from 100000 to 12500? 
		time_list = linspace(0, SIMTIME, len(x))
		plt.plot(time_list, x, "gx", label = "x-component")
		plt.xlabel("Time")
		plt.ylabel("Friction force")
		plt.title("Friction force as a function of time")
		#plt.axis([xmin, xmax, ymin, ymax])
		plt.grid(True)
		plt.legend(loc = "lower right")
		plt.savefig("time" + str(i) + ".svg")
	# End here if we're only plotting times.
	sys.exit()
else:
	for i in range(len(file_list)):
		fn = file_list[i]
		res = gather_data(fn)
		x,y,z = zip(*res)
		# Remove zeros.
		if nozeros:
			x = [j for j in x if j != 0] 
			y = [j for j in y if j != 0]
		averages.append( (mean(x), mean(y)) )
	write_data(averages, averages_fn)

# Make data ready for plotting.
avgs_x, avgs_y = zip(*averages)
# Generate force_list
if force_file:
    force_list = gather_data(force_list_fn)
    force_list = [a[0] for a in force_list]

# If the forces were not read from file only then create them using linspace.
else:
	force_list = linspace(linspace_min, linspace_max, linspace_steps).tolist()
	if verbose:
		print("Linspace is using parameters low = ", linspace_min, ", high = ", linspace_max, ", steps = ", linspace_steps, sep = "")

print(avgs_x)
print(force_list)
# If force list is still empty at this point, print a warning.
if len(force_list) == 0:
	print("Warning: Force list is empty!")
	print("Using a dirty hack.")
	force_list = range(len(file_list))

# If the data vectors are of different length at this point, the program must be halted.
if len(avgs_x) != len(avgs_y) or len(avgs_y) != len(force_list):
	print("Warning: data vectors are of different length.")
	print("Aborting program.")
	sys.exit()

# Fitting!
# https://towardsdatascience.com/basic-curve-fitting-of-scientific-data-with-python-9592244a2509
# https://scipy-lectures.org/intro/scipy/auto_examples/plot_curve_fit.html
def fitting_function(x, a, b):
	return a * x + b

# x_params is what we care about.
x_params, x_params_covariance = optimize.curve_fit(fitting_function, force_list, avgs_x)
y_params, y_params_covariance = optimize.curve_fit(fitting_function, force_list, avgs_y)

xx = linspace(linspace_min, linspace_max, 100)
yy1 = fitting_function(xx, *x_params)
yy2 = fitting_function(xx, *y_params)

# Plotting!
# https://matplotlib.org/tutorials/introductory/pyplot.html

# These settings can be modified.
plt.plot(force_list, avgs_x, "gx", label = "x-component")
plt.plot(force_list, avgs_y, "rx", label = "y-component")
plt.plot(xx, yy1, "k-", label = "linear fit(x)")
plt.plot(xx, yy2, "b-", label = "linear fit(y)")
plt.xlabel("Velocity")
plt.ylabel("Friction force")
plt.title("Friction force as a function of velocity")
#plt.axis([xmin, xmax, ymin, ymax])
plt.grid(True)
plt.legend(loc = "lower right")
plt.savefig("plots.svg")
