#! /usr/bin/env python

import friction_tools as ft
import constants
import numpy as np

simu = ft.FrictionSimulation()
velo_list = [0.01,0.005,0.001,0.0005,0.0001]
force = -0.5

for i in range(len(velo_list)):
    simu.continue_from_trajectory(filename = "beginning.traj")
    velo = velo_list[i]

    #interactions
    simu.create_interaction(['Fe', 'Fe'], strength = constants.iron_strength, equilibrium_distance = constants.iron_equilibrium_distance)

    simu.create_interaction(['Al', 'Al'], strength = constants.alu_strength, equilibrium_distance = constants.alu_equilibrium_distance)

    simu.create_interaction(['Fe', 'Al'], strength = constants.iron_alu_strength, equilibrium_distance = constants.iron_alu_equilibrium_distance)
    #indices
    top_slab = simu.get_indices_z_more_than(3.5)
    bottom_slab = simu.get_indices_z_less_than(0.0)
    bottom_indices = simu.get_indices_z_less_than(-3.5)
    top_indices = simu.get_indices_z_more_than(7.5)


    #stats and dynamics
    simu.create_dynamics(dt = 1.0, temperature = 300, coupled_indices = bottom_indices)
    simu.print_stats_during_simulation(interval = 50.0)

    simu.save_trajectory_during_simulation(interval=50.0, filename="simulation"+str(i) + ".traj")

    simu.gather_average_force_during_simulation(indices = top_slab, interval = 1.0, filename = "avr_force" + str(i) + ".txt")
 
    #happenings
    simu.fix_positions(bottom_slab, [True,True,True])

    simu.fix_velocities(top_indices, [velo, 0, 0], [True,False,False])
    simu.run_simulation(time = 500)

    simu.add_constant_force(top_slab, [0, 0, force])
    simu.run_simulation(time = 1500)

    ft.trajectory_to_xyz(filename_in="simulation"+str(i)+".traj", filename_out="simulation"+str(i)+".xyz")

