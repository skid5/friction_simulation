#! /usr/bin/env python

import  friction_tools as ft

simu = ft.FrictionSimulation()

#elements
simu.create_slab(element='Fe',xy_cells=3,z_cells=2,top_z=0.0)
simu.create_slab(element='Fe',xy_cells=3,z_cells=2,bottom_z=2.0)

#interactions
simu.create_interaction(['Fe','Fe'], strength=1.0, equilibrium_distance=2.375)

#indices
bottom_indices = simu.get_indices_z_less_than(-3.5)

#dynamics and stats
simu.create_dynamics(dt=1.0, temperature=500, coupled_indices=bottom_indices)
simu.print_stats_during_simulation(interval=50.0)
simu.save_trajectory_during_simulation(interval=50.0, filename="beginning.traj")
simu.gather_average_force_during_simulation()

#happenings
simu.run_simulation(time=500)
#ft.trajectory_to_xyz(filename_in="beginning.traj", filename_out="beginning.xyz")
