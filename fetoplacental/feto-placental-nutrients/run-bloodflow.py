#!/usr/bin/env python

# This routine reads in an arterial tree and solves Pressure-Resistance-Flow equations within this tree.
import os
import sys
import shutil

from reprosim.diagnostics import set_diagnostics_level
from reprosim.indices import perfusion_indices, get_ne_radius
from reprosim.geometry import append_units, define_node_geometry, define_1d_element_placenta, define_rad_from_geom, \
    add_matching_mesh, define_capillary_model, update_1d_elem_field, update_radius_by_order
from reprosim.repro_exports import export_1d_elem_geometry, export_node_geometry, export_1d_elem_field, \
    export_node_field, export_terminal_perfusion
from reprosim.pressure_resistance_flow import evaluate_prq, calculate_stats


def main(args=None):
    if args is None:
        args = sys.argv[1:]

    blood_flow_ml = float(args[0])
    remove_frac = float(args[1])
    blood_flow_model_units = blood_flow_ml / 60. * 1000.
    print(blood_flow_ml, blood_flow_model_units, remove_frac)

    ## Model parameterisation
    set_diagnostics_level(0)  # level 0 - no diagnostics; level 1 - only prints subroutine names (default); level 2 - prints subroutine names and contents of variables

    # define model geometry and indices
    anast_elem = 3
    anast_radius = 1.0

    perfusion_indices()
    define_node_geometry('../two_umb_arteries/sample_geometry/FullTree.ipnode')
    define_1d_element_placenta('../two_umb_arteries/sample_geometry/FullTree.ipelem',anast_elem)
    export_directory = 'output_flowbc_' + str(blood_flow_ml) + 'ml_min_' + str(remove_frac)

    if not os.path.exists(export_directory):
        os.makedirs(export_directory)
    mesh_type = 'full_plus_tube'
    # mesh_type: can be 'simple_tree' or 'full_plus_tube'. Simple_tree is the input
    ## arterial tree without any special features at the terminal level
    # 'full_plus_tube' creates a matching venous mesh and has arteries and
    ## veins connected by capillary units (capillaries are just tubes represented by an element)
    # define terminal units (this subroutine always needs to be called regardless of mesh_type
    append_units()

    # creates a mesh that converges (a venous mesh)
    umbilical_elem_option = 'single_umbilical_vein'
    umbilical_elements = [1, 2, 3, 4, 5]
    add_matching_mesh(umbilical_elem_option, umbilical_elements)

    # define radius by Strahler order in diverging (arterial mesh)
    s_ratio = 1.38  # rate of decrease in radius at each order of the arterial tree  1.38
    inlet_rad = 1.3  # inlet radius
    order_system = 'strahler'
    order_options = 'arterial'
    name = 'inlet'
    define_rad_from_geom(order_system, s_ratio, name, inlet_rad, order_options, '')
    # defines radius by STrahler order in converging (venous mesh)
    s_ratio_ven = 1.46  # rate of decrease in radius at each order of the venous tree 1.46
    inlet_rad_ven = 2.7  # inlet radius
    order_system = 'strahler'
    order_options = 'venous'
    first_ven_no = ''  # number of elements read in plus one
    last_ven_no = ''  # 2x the original number of elements + number of connections
    define_rad_from_geom(order_system, s_ratio_ven, first_ven_no, inlet_rad_ven, order_options, last_ven_no)

    ne_radius = get_ne_radius()
    update_1d_elem_field(ne_radius, anast_elem, anast_radius)

    num_convolutes = 10  # number of terminal convolute connections per intermediate villous
    num_generations = 3  # number of generations of symmetric intermediate villous trees
    num_parallel = 6  # number of parallel "capillary units" in a convolute
    capillary_model = 'erlich_resistance'
    define_capillary_model(num_convolutes, num_generations, num_parallel, capillary_model)
    # Call solve
    bc_type = 'flow'  # 'pressure' or 'flow'
    if bc_type == 'pressure':
        inlet_pressure = 6650  # Pa (~50mmHg)
        outlet_pressure = 2660  # Pa (~20mmHg)
        inlet_flow = 0  # set to 0 for bc_type = pressure;

    if bc_type == 'flow':
        inlet_pressure = 0
        outlet_pressure = 2660
        inlet_flow = blood_flow_model_units/2.0   # mm3/s (flow is divided between two inlets)

    rheology_type = 'pries_vessel'
    vessel_type = 'elastic'

    if remove_frac > 0.0:
        print("Occluding order 6 vessels")
        update_radius_by_order(6, 0.01, 'random', remove_frac)
    evaluate_prq(mesh_type, bc_type, rheology_type, vessel_type, inlet_flow, inlet_pressure, outlet_pressure)

    # this parameter is used to calculate the volume of vessels in the reconstructed tree that can't be
    # resolved in the images
    # set to 0 if not using
    image_voxel_size = 0
    calculate_stats(export_directory + '/terminal_flow_per_generation.csv', image_voxel_size, 1)

    ##export geometry
    group_name = 'perf_model'
    export_1d_elem_geometry(export_directory + '/full_tree.exelem', group_name)
    export_node_geometry(export_directory + '/full_tree.exnode', group_name)

    # export element field for radius
    field_name = 'radius_perf'
    ne_radius = get_ne_radius()
    export_1d_elem_field(ne_radius, export_directory + '/radius_perf.exelem', group_name, field_name)
    # export flow in each element
    field_name = 'flow'
    export_1d_elem_field(7, export_directory + '/flow_perf.exelem', group_name, field_name)
    # export node field for pressure
    field_name = 'pressure_perf'
    export_node_field(1, export_directory + '/pressure_perf.exnode', group_name, field_name)
    # Export terminal solution
    export_terminal_perfusion(export_directory + '/terminal.exnode', 'terminal_soln')
    shutil.move('micro_flow_results.out', export_directory + '/micro_flow_results.out')


if __name__ == '__main__':
    main()
