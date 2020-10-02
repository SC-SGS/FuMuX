#!/usr/bin/env python3

import argparse
import collections
import json
import logging
import math
import os
import shutil
import time
import xml.etree.ElementTree as ET  # dependency to parse ParaView pvd file

import meshio  # dependency for easy reading of vtu files
import numpy as np
import precice_future as precice

"""
Run the script with 
	python3 fumux.py case1/case1_single_tracer_fracture.pvd case1/precice-config.xml
"""

TIME_EPS = 1e-10
SLEEP_TIME = 0.5
DataSet = collections.namedtuple("DataSet", ["t", "path"])


class DumuxMesh:
    """"""

    def __init__(self, points=None):

        if points is not None:
            self.points = points
        else:
            self.points = []

        self.vertex_ids = []

        self.mesh_id = -1

        self.mesh_name = "unnamed"

    def __str__(self):
        numbers = ""
        if len(self.points) > 4:
            numbers = "{} {} ... {} {}".format(
                self.points[0], self.points[1], self.points[-2], self.points[-1]
            )
        else:
            for p in self.points:
                numbers += "{}".format(p)

        return "DuMuX mesh stats:\n  {} points\n  Coordinates: {}\n".format(
            len(self.points), numbers
        )


def create_file_list(path_to_pvd):
    base_path = path_to_pvd.rsplit("/", maxsplit=1)[0]
    if base_path == path_to_pvd:
        base_path = ""
    else:
        base_path += "/"

    print('Path to VTK files: "{}"'.format(base_path))

    tree = ET.parse(path_to_pvd)

    root = tree.getroot()
    file_list = []
    for dataset in root.iter("DataSet"):
        print(dataset.attrib["file"])
        file_list.append(
            DataSet(
                float(dataset.attrib["timestep"]), base_path + dataset.attrib["file"]
            )
        )

    print(file_list)
    return file_list


def extract_mesh(file_list, cell_type="quad"):

    print("Extracting mesh from first vtu file!")

    paraview_mesh = meshio.read(file_list[0].path)
    print(paraview_mesh)

    vis_points = np.zeros((len(paraview_mesh.cells[cell_type]), 3))
    for idx, point_list in zip(
        range(0, len(vis_points)), paraview_mesh.cells[cell_type]
    ):
        tmp = np.zeros((1, 3))
        for cell in point_list:
            for point in paraview_mesh.points[cell]:
                tmp = np.add(tmp, point)
        vis_points[idx] = tmp / 4.0

    print(vis_points)

    return DumuxMesh(vis_points)


def extract_cell_data(file_name, data_label, cell_type="quad"):
    vtu_data = meshio.read(file_name.path)

    print(
        "{} has shape {} ".format(
            data_label, np.shape(vtu_data.cell_data[cell_type][data_label])
        )
    )
    # print( "Extracting mesh from first vtu file!")
    return vtu_data.cell_data[cell_type][data_label]


def extract_data_from_vtu(vtu_data, data_label, cell_type="quad"):
    old_shape = np.shape(vtu_data.cell_data[cell_type][data_label])
    return np.reshape(vtu_data.cell_data[cell_type][data_label], (old_shape[0],))


def main():
    args = parse_args()
    vtu_file_list = create_file_list(args.pvd_filename)
    print(vtu_file_list)
    my_mesh = extract_mesh(vtu_file_list)
    my_mesh.mesh_name = "DumuxMesh"
    print(my_mesh)

    print("Starting visualization coupling...")

    configuration_file_name = args.precice_config
    participant_name = "Dumux"

    ### Create preCICE interfacemesh_id
    interface = precice.Interface(participant_name, 0, 1)
    interface.configure(configuration_file_name)

    # dimensions = interface.get_dimensions()
    dimensions = 3

    ### Get mesh ID for preCICE
    my_mesh.mesh_id = interface.get_mesh_id(my_mesh.mesh_name)

    ### Set mesh for preCICE
    my_mesh.vertex_ids = interface.set_mesh_vertices(my_mesh.mesh_id, my_mesh.points)
    print(my_mesh.vertex_ids)

    ### Get ID of data
    pressure_id = interface.get_data_id("Pressure", my_mesh.mesh_id)
    concentration_id = interface.get_data_id("Concentration", my_mesh.mesh_id)

    fileNumber = 0

    print("preCICE initialized. Begin coupling iterations...")

    t = 0
    dt = interface.initialize()

    # while interface.is_coupling_ongoing():
    for vtu_file in vtu_file_list:
        if not interface.is_coupling_ongoing():
            break

        vtu_data = meshio.read(vtu_file.path)
        print("VTU data: {}".format(vtu_file))
        print(
            "Time difference: abs(t - vtu_file.t) / t = {} ".format(
                abs(t - vtu_file.t) / (t + TIME_EPS)
            )
        )

        assert (
            abs(t - vtu_file.t) / (t + TIME_EPS) < TIME_EPS
        ), "Time does not fit with time from data set!\n    Coupling time: {}\n    Data set time: {}".format(
            t, vtu_data.t
        )

        pressure = extract_data_from_vtu(vtu_data, "rho")
        concentration = extract_data_from_vtu(vtu_data, "X^tracer_0")

        # print("Concentration: ", concentration[my_mesh.vertex_ids])
        print("fileNumber: ", fileNumber)

        interface.write_block_scalar_data(pressure_id, my_mesh.vertex_ids, pressure)
        interface.write_block_scalar_data(
            concentration_id, my_mesh.vertex_ids, concentration
        )

        dt = interface.advance(dt)
        print("Sleeping for a {}s".format(SLEEP_TIME))
        time.sleep(SLEEP_TIME)

        fileNumber += 1

        t += dt
        print("Simulation time:  {}".format(t))

    print(
        "This process does only terminate automatically after the final simulation time has been reached AND the visualization program is closed."
    )
    print(
        "No more new files to process! If the final simulation time has not been reached, you have to kill the process now."
    )

    interface.finalize()
    print("Closing visualization...")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Read meshes, partition them and write them out in internal format."
    )
    # parser.add_argument("in_meshname", nargs="+", help="The meshes used as input")
    parser.add_argument(
        "pvd_filename", type=str, help="Path to ParaView PVD file to parse"
    )
    parser.add_argument(
        "precice_config", type=str, help="Path to preCICE XML configuration"
    )
    parser.add_argument(
        "--log",
        "-l",
        dest="logging",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Set the log level. Default is INFO",
    )
    return parser.parse_args()


if __name__ == "__main__":
    main()
