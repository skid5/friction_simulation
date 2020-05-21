import random
sep = "\t"
n_files = 20
n_datapoints = 10
len_datavec = 3
for i in range(n_files):
    fn = "avr_force" + str(i) + ".txt"
    fp = open(fn, "w")
    for j in range(n_datapoints):
        for k in range(len_datavec):
            r = random.random()
            if r < 0.5:
                r = 0.0
            fp.write(sep + str(r))
        fp.write("\n")
    fp.close()
