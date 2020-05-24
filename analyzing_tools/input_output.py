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
