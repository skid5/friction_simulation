from numpy import mean
n_simulations = 5
for i in xrange(n_simulations):
    fp = open('avr_force' + str(i) + '.txt', 'r')
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

    fp_x = open('average_x.txt', 'a')
    fp_y = open('average_y.txt', 'a')

    average_x = mean(x)
    average_y = mean(y)

    fp_x.write(str(average_x) + '\n')
    fp_y.write(str(average_y) + '\n')

    fp_x.close()
    fp_y.close()
