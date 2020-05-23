#! /usr/bin/env python

import friction_tools as ft
import constants as cts
import time

simu = ft.FrictionSimulation()

SIMTIME = 10**5
velo = 2*10**-5
list = range(10)

t0 = time.time()

for i in list:
    force = (10 - i)*cts.F_min
    
    simu.continue_from_trajectory(filename = "beginning.traj")

    #interactions
    simu.create_interaction(['Fe','Fe'], strength=cts.bond_strength_iron, equilibrium_distance=cts.eq_dist_iron)

    simu.create_interaction(['Cu','Cu'], strength=cts.bond_strength_copper, equilibrium_distance=cts.eq_dist_copper)

    simu.create_interaction(['Fe','Cu'], strength=cts.bond_strength_between, equilibrium_distance=cts.eq_dist_between)

    #indices
    top_indices = simu.get_indices_z_more_than(7)
    bottom_indices = simu.get_indices_z_less_than(-3.5)
    top_slab = simu.get_indices_z_more_than(1)
    bottom_slab = simu.get_indices_z_less_than(1)


    #stats and dynamics
    simu.create_dynamics(dt = 4, temperature = 300, coupled_indices = bottom_slab, strength=cts.thermostat)
    simu.print_stats_during_simulation(interval = 50.0)

    simu.save_trajectory_during_simulation(interval=50.0, filename="simulation"+str(i) + ".traj")

    simu.gather_average_force_during_simulation(indices = top_slab, interval = 8, filename = "avr_force" + str(i) + ".txt")

    simu.gather_energy_and_temperature_during_simulation(interval=50, filename='energy'+str(i)+'.txt')
 
    #happenings
    simu.fix_positions(bottom_indices, [True,True,True])
    simu.fix_positions(top_indices, [False,False,False])
    simu.fix_velocities(top_indices, [velo, 0, 0], [True,True,False])
    simu.add_constant_force(top_indices, [0, 0, force])
    
    simu.run_simulation(steps= SIMTIME)

    ft.trajectory_to_xyz(filename_in="simulation"+str(i)+".traj", filename_out="simulation"+str(i)+".xyz")


t1=time.time()
print "time taken {ti} s".format(ti=str(int(t1-t0)))

