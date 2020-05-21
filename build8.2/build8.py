#! /usr/bin/env python

import friction_tools as ft
import constants as cts
import numpy as np

simu = ft.FrictionSimulation()
velo = 4*10**(-4) 
force = -10**(-3)
simu.continue_from_trajectory(filename = "beginning.traj")

#interactions
simu.create_interaction(['Fe','Fe'], strength=cts.bond_strength, equilibrium_distance=cts.eq_dist)

simu.create_interaction(['Cu','Cu'], strength=cts.bond_strength, equilibrium_distance=cts.eq_dist)

simu.create_interaction(['Fe','Cu'], strength=0.1*cts.bond_strength, equilibrium_distance=cts.eq_dist)

#indices
top_slab = simu.get_indices_z_more_than(3.5)
bottom_slab = simu.get_indices_z_less_than(0.0)
bottom_indices = simu.get_indices_z_less_than(-3.5)
top_indices = simu.get_indices_z_more_than(7)


#stats and dynamics
simu.create_dynamics(dt = 2.0, temperature = 300, coupled_indices = bottom_slab, strength=cts.thermal_conduct)
simu.print_stats_during_simulation(interval = 50.0)

simu.save_trajectory_during_simulation(interval= 50.0)

simu.gather_average_force_during_simulation(indices = top_slab, interval = 50.0)
 
#happenings
simu.set_temperature(300)
simu.fix_positions(bottom_indices, [True,True,True])
simu.fix_positions(top_indices, [False,False,False])

simu.fix_velocities(top_indices, [velo, 0, 0], [True,True,False])

simu.add_constant_force(top_slab, [0, 0, force])
simu.run_simulation(steps= 10**4)

ft.trajectory_to_xyz()

