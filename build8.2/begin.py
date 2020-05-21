#! /usr/bin/env python

import  friction_tools as ft

simu = ft.FrictionSimulation()
import constants as cts

#elements
simu.create_slab(element='Fe',xy_cells=3,z_cells=2,top_z=0.0)
simu.create_slab(element='Cu',xy_cells=3,z_cells=2,bottom_z=2.39)

#interactions
simu.create_interaction(['Fe','Fe'], strength=cts.bond_strength, equilibrium_distance=cts.eq_dist)

simu.create_interaction(['Cu','Cu'], strength=cts.bond_strength, equilibrium_distance=cts.eq_dist)

simu.create_interaction(['Fe','Cu'], strength=0.1*cts.bond_strength, equilibrium_distance=cts.eq_dist)

#indices
bottom_indices = simu.get_indices_z_less_than(-3.5)
top_indices = simu.get_indices_z_more_than(7)
bottom_slab = simu.get_indices_z_less_than(0.0)

#dynamics and stats
simu.create_dynamics(dt=2.0, temperature=300, coupled_indices=bottom_slab, strength=cts.thermal_conduct)
simu.print_stats_during_simulation(interval=50.0)
simu.save_trajectory_during_simulation(interval=50.0, filename="beginning.traj")
simu.fix_positions(bottom_indices, [True,True,True])
simu.fix_positions(top_indices, [True,True,True])

#happenings
simu.set_temperature(300)
simu.run_simulation(steps=500)
ft.trajectory_to_xyz(filename_in="beginning.traj", filename_out="beginning.xyz")
