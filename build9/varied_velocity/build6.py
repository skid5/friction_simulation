#! /usr/bin/env python

import friction_tools as ft
import numpy as np
import time

t0 = time.time()

SIMTIME = 5000
k = 10**-3

simu = ft.FrictionSimulation()
velo_list = [10*k,9*k,8*k,7*k,6*k,5*k,4*k,3*k,2*k,1*k]
force = -0.1

for i in range(len(velo_list)):
    simu.continue_from_trajectory(filename = "beginning.traj")
    velo = velo_list[i]

    #interactions
    simu.create_interaction(['Fe', 'Fe'], strength = 1.0, equilibrium_distance = 2.39)

    simu.create_interaction(['Ni', 'Ni'], strength = 1.0,  equilibrium_distance = 2.39)

    simu.create_interaction(['Fe', 'Ni'], strength = 0.1, equilibrium_distance = 2.39)
    #indices
    top_slab = simu.get_indices_z_more_than(3.5)
    bottom_slab = simu.get_indices_z_less_than(0.0)
    bottom_indices = simu.get_indices_z_less_than(-3.5)
    top_indices = simu.get_indices_z_more_than(3.5)


    #stats and dynamics
    simu.create_dynamics(dt = 1.0, temperature = 300, coupled_indices = bottom_indices)
    simu.print_stats_during_simulation(interval = 50.0)

    simu.save_trajectory_during_simulation(interval=50.0, filename="simulation"+str(i) + ".traj")

    simu.gather_average_force_during_simulation(indices = top_slab, interval = 1.0, filename = "avr_force" + str(i) + ".txt")
 
    #happenings
    simu.set_temperature(300)
    simu.fix_positions(bottom_indices, [True,True,True])

    simu.fix_velocities(top_indices, [velo, 0, 0], [True,False,False])
    simu.add_constant_force(top_slab, [0, 0, force])

    simu.run_simulation(time = SIMTIME)

    ft.trajectory_to_xyz(filename_in="simulation"+str(i)+".traj", filename_out="simulation"+str(i)+".xyz")

t1=time.time()
print "time taken {ti} s".format(ti=str(int(t1-t0)))

