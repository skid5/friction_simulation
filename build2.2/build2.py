#! /usr/bin/env python

import  friction_tools as ft

simu = ft.FrictionSimulation()
simu.continue_from_trajectory(filename="beginning.traj")


#interactions
simu.create_interaction(['Au','Au'], strength=1.0, equilibrium_distance=2.375)
simu.create_interaction(['Ag','Ag'], strength=1.0, equilibrium_distance=1.781)
simu.create_interaction(['Au','Ag'], strength=0.1, equilibrium_distance=2.0)

#indices
au_indices = simu.get_indices_by_element('Au')
ag_indices = simu.get_indices_by_element('Ag')
bottom_indices = simu.get_indices_z_less_than(-3.5)


#dynamics and stats
simu.create_dynamics(dt=1.0, temperature=500, coupled_indices=bottom_indices)
simu.print_stats_during_simulation(interval=50.0)
simu.save_trajectory_during_simulation(interval=50.0)

simu.gather_average_force_during_simulation(indices=ag_indices)

#happenings
simu.fix_positions(bottom_indices, [True,True,True])
simu.add_constant_force(ag_indices, [0, 0, -0.1])
simu.run_simulation(time=2000.0)

simu.fix_velocities(ag_indices, [0.005, 0, 0], [True,False,False])
simu.run_simulation(time=5000.0)

ft.trajectory_to_xyz()
