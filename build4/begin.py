#! /usr/bin/env python

import  friction_tools as ft
import constants

simu = ft.FrictionSimulation()

#elements
simu.create_slab(element='Fe',xy_cells=3,z_cells=2,top_z=0.0)
simu.create_slab(element='Fe',xy_cells=4,z_cells=2,bottom_z=2.0)

#interactions
simu.create_interaction(['Fe', 'Fe'], strength = constants.iron_strength, equilibrium_distance = constants.iron_equilibrium_distance)

#indices
top_slab = simu.get_indices_z_more_than(2.0)
bottom_slab = simu.get_indices_z_less_than(0.0)
bottom_indices = simu.get_indices_z_less_than(-3.5)

#dynamics and stats
simu.create_dynamics(dt = 1.0, temperature = 500, coupled_indices = bottom_indices)
simu.print_stats_during_simulation(interval = 50.0)
simu.save_trajectory_during_simulation(interval = 50.0, filename = "beginning.traj")
simu.gather_average_force_during_simulation()

#happenings
simu.set_temperature(500)

simu.run_simulation(time=2000.0)
