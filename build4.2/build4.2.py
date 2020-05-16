#! /usr/bin/env python

import friction_tools as ft
import constants

simu = ft.FrictionSimulation()
force_list = [-0.5, -0.1, -0.05, -0.01]
for i in range(len(force_list)):
    simu.continue_from_trajectory(filename = "beginning.traj")
    force = force_list[i]

    #interactions
    simu.create_interaction(['Fe', 'Fe'], strength = constants.iron_strength, equilibrium_distance = constants.iron_equilibrium_distance)

    #indices
    top_slab = simu.get_indices_z_more_than(2.0)
    bottom_slab = simu.get_indices_z_less_than(0.0)
    bottom_indices = simu.get_indices_z_less_than(-3.5)


    #stats and dynamics
    simu.create_dynamics(dt = 1.0, temperature = 300, coupled_indices = bottom_indices)
    simu.print_stats_during_simulation(interval = 50.0)

    simu.save_trajectory_during_simulation(interval=50.0, filename="simulation"+str(i) + ".traj")

    simu.gather_average_force_during_simulation(indices = top_slab, interval = 1.0, filename = "avr_force" + str(i) + ".txt")
 
    #happenings
    simu.fix_positions(bottom_slab, [True,True,True])
    simu.add_constant_force(top_slab, [0, 0, force])
    simu.run_simulation(time = 100.0)

    simu.fix_velocities(top_slab, [0.005, 0, 0], [True,False,False])
    simu.run_simulation(time = 500)

    ft.trajectory_to_xyz(filename_in="simulation"+str(i)+".traj", filename_out="simulation"+str(i)+".xyz")

