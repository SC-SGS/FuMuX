import argparse, logging, math, os
import numpy as np
import shutil
from ctypes import *
import mesh_io
import json
import precice

'''
Run the script with 
	python3 visual.py data-0.vtk


'''

class Mesh:
	"""
	A Mesh consists of:
		- Points: A list of tuples of floats representing coordinates of points
		- Cells: A list of tuples of ints representing mesh elements
		- Pointdata: A list of floats representing data values at the respective point
	"""
	def __init__(self, points = None, cells = None, cell_types = None, pointdata = None):
		if points is not None:
			self.points = points
		else:
			self.points = []
		if cells is not None:
			assert(cell_types is not None)
			self.cells = cells
			self.cell_types = cell_types
		else:
			self.cells = []
			self.cell_types = []
		if pointdata is not None:
			self.pointdata = pointdata
		else:
			self.pointdata = []

		def __str__(self):
			return "Mesh with {} Points and {} Cells ({} Cell Types)".format(len(self.points), len(self.cells), len(self.cell_types))




def read_mesh(filename, tag, datadim=1):
	points, cells, cell_types, pointdata = mesh_io.read_mesh(filename, tag)
	#print("Points: ", len(points))
	#print("Point data: ", pointdata)
	return Mesh(points, cells, cell_types, pointdata)



def main():
	print("Starting visualization coupling...")

	configuration_file_name = "precice-config.xml"
	participant_name = "visualizer"
	mesh_name = "visMesh"
	
	### Create preCICE interface
	#interface = precice.Interface(participant_name, configuration_file_name, 0, 1)

	#dimensions = interface.get_dimensions()
	dimensions = 3

	### Get mesh ID for preCICE
	#mesh_id = interface.get_mesh_id(mesh_name)


	args = parse_args()
	logging.basicConfig(level=getattr(logging, args.logging))
	if len(args.in_meshname) > 1 and args.out_meshname:
		logging.warn("--out ignored")
	mesh_names = args.in_meshname
	for mesh_name in mesh_names:
		assert os.path.isfile(mesh_name), ("Invalid filename: "  + mesh_name)
	mesh = read_mesh(mesh_name, args.tag, args.datadim)
	print("Points: ", len(mesh.points))
	#print("Points: ", mesh.points)
	#print("Point data: ", mesh.celldata)
	n = len(mesh.points)

	vertices_mesh = np.zeros((n, dimensions))
	pressure = []
	concentration = []
	
	for i in range(0,n):
		for j in range(0,dimensions):
			vertices_mesh[i][j] = mesh.points[i][j]
		pressure.append(mesh.pointdata[i])
		concentration.append(mesh.pointdata[i])


	# This is the vertices for coupling and also the correct array sizes for the data. Data gets overwritten
	print("Vertices_mesh: ", vertices_mesh)
	#print("Pressure: ", pressure)
	#print("concentration: ", concentration)
	
	### Set mesh for preCICE
	#data_indices = interface.set_mesh_vertices(mesh_id, vertices_mesh)

	### Get ID of data
	#pressure_id = interface.get_data_id("pressure",mesh_id)
	#concentration_id = interface.get_data_id("concentration",mesh_id)

	#dt = interface.initialize()
	fileNumber = 0

	print("preCICE initialized. Begin coupling iterations...")

	#while interface.is_coupling_ongoing():
	while fileNumber < 3:
		
		mesh_name = "case-" + str(fileNumber) + ".vtu"
		mesh = read_mesh(mesh_name, args.tag, args.datadim)

		for i in range(0,n):
			pressure[i] = mesh.pointdata[i]
			concentration[i] = mesh.pointdata[i]

		print("Pressure: ", pressure)
		print("fileNumber: ", fileNumber)
	
		#interface.write_block_scalar_data(pressure_id, data_indices, pressure)
    	#interface.write_block_scalar_data(concentration_id, data_indices, concentration)

		#dt = interface.advance(dt)

		fileNumber += 1

	#interface.finalize()
	print("Closing visualization...")


def parse_args():
	parser = argparse.ArgumentParser(description=
                                     "Read meshes, partition them and write them out in internal format.")
	parser.add_argument("in_meshname", nargs="+", help="The meshes used as input")
	parser.add_argument("--tag", "-t", dest="tag", default=None,
		help="The PointData tag for vtk meshes")
	parser.add_argument("--log", "-l", dest="logging", default="INFO", 
		choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
		help="Set the log level. Default is INFO")
	parser.add_argument("--datadim", "-d", dest="datadim", default=1, type=int, help="Dimensions of the function. Default is 1 (Scalar function.")
	return parser.parse_args()

if __name__ == "__main__":
	main()


