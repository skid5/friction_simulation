#! /usr/bin/env python

import  friction_tools as ft

simu = ft.FrictionSimulation()
force_list = [-0.5, -0.1, -0.05, -0.01]
for i in range(len(force_list)-1):
    simu.continue_from_trajectory(filename="beginning.traj")


    #interactions
    simu.create_interaction(['Au','Au'], strength=1.0, equilibrium_distance=2.375)
    simu.create_interaction(['Ag','Ag'], strength=1.0, equilibrium_distance=1.781)
    simu.create_interaction(['Au','Ag'], strength=0.1, equilibrium_distance=2.0)

    #indices
    au_indices = simu.get_indices_by_element('Au')
    ag_indices = simu.get_indices_by_element('Ag')
    bottom_indices = simu.get_indices_z_less_than(-3.5)
    #top_indices = simu.get_indices_z_more_than(5.5)


    #stats and dynamics
    if i == 0:
        simu.create_dynamics(dt=1.0, temperature=300, coupled_indices=bottom_indices)
        simu.print_stats_during_simulation(interval=50.0)

    simu.save_trajectory_during_simulation(interval=50.0, filename="simulation"+str(i) + ".traj")

    simu.gather_average_force_during_simulation(indices=ag_indices, filename="avr_force"+str(i)+".txt")
 
    #happenings
    simu.fix_positions(bottom_indices, [True,True,True])
    simu.add_constant_force(ag_indices, [0, 0, force_list[i]])
    simu.run_simulation(time=100.0)

    simu.fix_velocities(ag_indices, [0.005, 0, 0], [True,False,False])
    simu.run_simulation(time=100.0)

    ft.trajectory_to_xyz(filename_in="simulation"+str(i)+".traj", filename_out="simulation"+str(i)+".xyz")


