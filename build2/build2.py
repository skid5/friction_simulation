#! /usr/bin/env python

import  friction_tools as ft

simu = ft.FrictionSimulation()

#elements
simu.create_slab(element='Au',xy_cells=3,z_cells=2,top_z=0.0)
simu.create_slab(element='Ag',xy_cells=4,z_cells=2,bottom_z=3.0)

#interactions
simu.create_interaction(['Au','Au'], strength=1.0, equilibrium_distance=2.375)
simu.create_interaction(['Ag','Ag'], strength=1.0, equilibrium_distance=1.781)
simu.create_interaction(['Au','Ag'], strength=0.1, equilibrium_distance=2.0)

list_atoms()

#indices
au_indices = simu.get_indices_by_element('Au')
ag_indices = simu.get_indices_by_element('Ag')
bottom_indices = simu.get_indices_z_less_than(-3.5)
top_indices = simu.get_indices_z_more_than(8.0)

#dynamics and stats
simu.create_dynamics(dt=1.0, temperature=500)
simu.print_stats_during_simulation(interval=50.0)
simu.save_trajectory_during_simulation(interval=50.0)

simu.gather_average_force_during_simulation(interval=1.0, indices=[1], filename='force.txt')

#happenings
simu.fix_positions(au_indices, [True,True,True])
simu.add_constant_force(ag_indices, [0, 0, -0.02])
simu.fix_velocities(ag_indices, [0.005, 0, 0], [True,False,False])
simu.set_temperature(500)

simu.run_simulation(time=2000.0)
ft.trajectory_to_xyz()
