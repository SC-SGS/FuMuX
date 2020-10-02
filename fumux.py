#!/usr/bin/env python3

import argparse, logging, math
import numpy as np
import meshio  # dependency for easy reading of vtu files
import precice_future as precice
import xml.etree.ElementTree as ET  # dependency to parse ParaView pvd file
import collections
import time

"""
Run the script with 
	python3 fumux.py case1/case1_single_tracer_fracture.pvd case1/precice-config.xml
"""

TIME_EPS = 1e-10
COMPARE_EPS = 1e-12
SLEEP_TIME = 0.5
DataSet = collections.namedtuple("DataSet", ["t", "path"])

GridOptions = collections.namedtuple("GridOptions", "id, x0, x1")


gridOptions = [
    GridOptions(
        id=0, x0=np.array([0.500, 0.000, 0.000]), x1=np.array([0.500, 1.000, 1.000])
    ),
    GridOptions(
        id=1, x0=np.array([0.000, 0.500, 0.000]), x1=np.array([1.000, 0.500, 1.000])
    ),
    GridOptions(
        id=2, x0=np.array([0.000, 0.000, 0.500]), x1=np.array([1.000, 1.000, 0.500])
    ),
    GridOptions(
        id=3, x0=np.array([0.750, 0.500, 0.500]), x1=np.array([0.750, 1.000, 1.000])
    ),
    GridOptions(
        id=5, x0=np.array([0.500, 0.750, 0.500]), x1=np.array([1.000, 0.750, 1.000])
    ),
    GridOptions(
        id=4, x0=np.array([0.500, 0.500, 0.750]), x1=np.array([1.000, 1.000, 0.750])
    ),
    GridOptions(
        id=7, x0=np.array([0.625, 0.500, 0.500]), x1=np.array([0.625, 0.750, 0.750])
    ),
    GridOptions(
        id=6, x0=np.array([0.500, 0.625, 0.500]), x1=np.array([0.750, 0.625, 0.750])
    ),
    GridOptions(
        id=8, x0=np.array([0.500, 0.500, 0.625]), x1=np.array([0.750, 0.750, 0.625])
    ),
]


class DumuxMesh:
    """"""

    def __init__(self, grid_option, points=None):

        self.grid_option = grid_option
        self.gridDimensions = -1

        if points is not None:
            self.points = points
        else:
            self.points = []

        self.vertex_ids = []

        self.mesh_id = -1
        self.data_id = -1

        self.global_to_mesh_id = []

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

        mapped_ids = ""
        if len(self.global_to_mesh_id) > 4:
            mapped_ids = "{} {} ... {} {}".format(
                self.global_to_mesh_id[0],
                self.global_to_mesh_id[1],
                self.global_to_mesh_id[-2],
                self.global_to_mesh_id[-1],
            )
        else:
            for p in self.points:
                mapped_ids += "{}".global_to_mesh_id(p)

        return "DuMuX mesh stats:\n  Mesh name: {}\n  Mesh id: {}\n  Data id: {}\n  {} points\n  Coordinates: {}\n  Mapped ids: {}\n".format(
            self.mesh_name,
            self.mesh_id,
            self.data_id,
            len(self.points),
            numbers,
            mapped_ids,
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


def extract_meshes(file_list, cell_type="quad"):

    print("Extracting meshes from first vtu file: {}".format(file_list[0].path))

    paraview_mesh = meshio.read(file_list[0].path)
    print(paraview_mesh)

    meshes = []
    for grid_option in gridOptions:
        meshes.append(DumuxMesh(grid_option))
    print(meshes)

    for idx, cell_point_list in zip(
        range(0, len(paraview_mesh.cells[cell_type])), paraview_mesh.cells[cell_type]
    ):
        tmp = np.zeros(3)
        for cell_point_id in cell_point_list:

            # Check if all points of cell are on surface
            for mesh in meshes:
                p0 = mesh.grid_option.x0
                p1 = mesh.grid_option.x1
                cell_is_in_grid = True
                for point in paraview_mesh.points[cell_point_id]:
                    point = np.reshape(point[0], (3))

                    # Is cell in current mesh?
                    if not (
                        (
                            min(p0[0], p1[0]) - COMPARE_EPS
                            <= point[0]
                            <= max(p0[0], p1[0]) + COMPARE_EPS
                        )
                        and (
                            min(p0[1], p1[1]) - COMPARE_EPS
                            <= point[1]
                            <= max(p0[1], p1[1]) + COMPARE_EPS
                        )
                        and (
                            min(p0[2], p1[2]) - COMPARE_EPS
                            <= point[2]
                            <= max(p0[2], p1[2]) + COMPARE_EPS
                        )
                    ):
                        cell_is_in_grid = False
                        break

                if cell_is_in_grid == True:
                    for point in paraview_mesh.points[cell_point_id]:
                        point = np.reshape(point[0], (3))
                        tmp = np.add(tmp, point)
                    tmp /= len(paraview_mesh.points[cell_point_id])
                    mesh.points.append(tmp)
                    mesh.global_to_mesh_id.append(idx)
                    break

    sum_of_points = 0
    for mesh in meshes:
        sum_of_points += len(mesh.points)

        mesh.points = np.array(mesh.points)
        mesh.global_to_mesh_id = np.array(mesh.global_to_mesh_id)

        mesh.mesh_name = "DumuxMesh{}".format(mesh.grid_option.id)

    print("Total number of points: ", sum_of_points)

    return meshes


def extract_cell_data(file_name, data_label, cell_type="quad"):
    vtu_data = meshio.read(file_name.path)
    return vtu_data.cell_data[cell_type][data_label]


def extract_data_from_vtu(vtu_data, data_label, cell_type="quad"):
    old_shape = np.shape(vtu_data.cell_data[cell_type][data_label])
    return np.reshape(vtu_data.cell_data[cell_type][data_label], (old_shape[0],))


def main():
    args = parse_args()
    vtu_file_list = create_file_list(args.pvd_filename)
    print(vtu_file_list)
    my_meshes = extract_meshes(vtu_file_list, "triangle")

    print("Starting visualization coupling...")

    configuration_file_name = args.precice_config  # "precice-config.xml"
    participant_name = "Dumux"

    ### Create preCICE interfacemesh_id
    interface = precice.Interface(participant_name, 0, 1)
    interface.configure(configuration_file_name)

    dimensions = 3

    ### Get mesh IDs and data IDs from preCICE and announce vertices
    for mesh in my_meshes:
        mesh.mesh_id = interface.get_mesh_id(mesh.mesh_name)
        mesh.vertex_ids = interface.set_mesh_vertices(mesh.mesh_id, mesh.points)

        ### Get ID of data
        mesh.data_id = interface.get_data_id(
            "Concentration{}".format(mesh.grid_option.id), mesh.mesh_id
        )

        print(mesh)

    fileNumber = 0

    print("preCICE initialized. Begin coupling iterations...")

    t = 0
    dt = interface.initialize()

    for vtu_file in vtu_file_list:
        if not interface.is_coupling_ongoing():
            break

        print("File no. {}".format(fileNumber))

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

        concentration = extract_data_from_vtu(vtu_data, "X^tracer_0", "triangle")
        for mesh in my_meshes:
            # print("Concentration: ", concentration[mesh.global_to_mesh_id] )

            interface.write_block_scalar_data(
                mesh.data_id, mesh.vertex_ids, concentration[mesh.global_to_mesh_id]
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
