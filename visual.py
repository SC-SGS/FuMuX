#!/usr/bin/env python3

import argparse, logging, math, os
import numpy as np
import shutil
#from ctypes import *
#import mesh_io
import meshio #dependency for easy reading of vtu files
import json
import precice_future as precice
import xml.etree.ElementTree as ET #dependency to parse ParaView pvd file 
import collections
import time

'''
Run the script with 
	python3 visual.py case1/case1_single_tracer_fracture.pvd case1/precice-config.xml
'''

TIME_EPS = 1e-10
SLEEP_TIME = 1.5
DataSet = collections.namedtuple("DataSet", [ "t", "path" ])

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

class DumuxMesh:
	"""
	"""
	def __init__(self, points = None ):


		if points is not None:
			self.points = points
		else:
			self.points = []
		
		self.vertex_ids = []

		self.mesh_id = -1

		self.mesh_name = "unnamed"

	def __str__(self):
		numbers = ""
		if len( self.points) > 4:
			numbers = "{} {} ... {} {}".format( self.points[0], self.points[1], self.points[-2], self.points[-1] )
		else:
			for p in self.points:
				numbers += "{}".format( p )

		return "DuMuX mesh stats:\n  {} points\n  Coordinates: {}\n".format( len( self.points), numbers ) 





def create_file_list( path_to_pvd ):
	base_path = path_to_pvd.rsplit("/", maxsplit=1)[0]
	if ( base_path == path_to_pvd ):
		base_path = ""
	else:
		base_path += "/"

	print( "Path to VTK files: \"{}\"".format(base_path) )

	tree = ET.parse( path_to_pvd )

	root = tree.getroot()
	file_list = []
	for dataset in root.iter("DataSet"):
		print( dataset.attrib["file"] )
		#file_list.append( base_path + dataset.attrib["file"] )
		file_list.append( DataSet( float(dataset.attrib["timestep"]), base_path + dataset.attrib["file"]) )

	print( file_list )
	return file_list


def extract_mesh( file_list, cell_type = "quad" ):

	print( "Extracting mesh from first vtu file!")

	paraview_mesh = meshio.read( file_list[0].path )
	print( paraview_mesh )

	vis_points = np.zeros(( len(paraview_mesh.cells[cell_type]), 3)) 
	for idx, point_list in zip( range(0, len(vis_points)), paraview_mesh.cells[cell_type] ):
		tmp = np.zeros((1,3))
		for cell in point_list:
			for point in paraview_mesh.points[ cell ]:
		#print( dataset[0] )
			#print( cell )
				#print( point )
				#print( paraview_mesh.points[ cell ] )
				#print( paraview_mesh.points[ point_id ] )
				tmp = np.add(tmp, point )
		vis_points[idx] = tmp
	
	return DumuxMesh( vis_points )

def extract_cell_data( file_name, data_label, cell_type = "quad" ):

	#print( "Extracting mesh from first vtu file!")
	vtu_data = meshio.read( file_name.path )
	#print( paraview_mesh )
	#print( vtu_data.cell_data[cell_type][data_label] )

	print( "{} has shape {} ".format( data_label, np.shape( vtu_data.cell_data[cell_type][data_label] ) ) )
	return vtu_data.cell_data[cell_type][data_label]

def extract_data_from_vtu( vtu_data, data_label, cell_type = "quad" ):
	#print( vtu_data.cell_data[cell_type][data_label] )
	#print( "{} has shape {} ".format( data_label, np.shape( vtu_data.cell_data[cell_type][data_label] ) ) )
	old_shape = np.shape( vtu_data.cell_data[cell_type][data_label] ) 
	return np.reshape( vtu_data.cell_data[cell_type][data_label], (old_shape[0],) ) 

def main():
	args = parse_args()
	vtu_file_list = create_file_list( args.pvd_filename )
	print( vtu_file_list )
	my_mesh = extract_mesh( vtu_file_list )
	my_mesh.mesh_name = "DumuxMesh"
	print( my_mesh )

	#print( extract_cell_data( vtu_file_list[0], "x^tracer_0" ) )


	print("Starting visualization coupling...")

	configuration_file_name = args.precice_config #"precice-config.xml"
	participant_name = "Dumux"
	#mesh_name = "DumuxMesh"
	
	### Create preCICE interfacemesh_id
	interface = precice.Interface(participant_name, 0, 1)
	interface.configure(configuration_file_name)

	#dimensions = interface.get_dimensions()
	dimensions = 3

	### Get mesh ID for preCICE
	my_mesh.mesh_id = interface.get_mesh_id( my_mesh.mesh_name )


	#logging.basicConfig(level=getattr(logging, args.logging))
	#if len(args.in_meshname) > 1 and args.out_meshname:
	#	logging.warn("--out ignored")
	#mesh_names = args.in_meshname
	#for mesh_name in mesh_names:vtu_data = meshio.read( file_name.path )
	#	assert os.path.isfile(mesh_name), ("Invalid filename: "  + mesh_name)
	#mesh = read_mesh(mesh_name, args.tag, args.datadim)
	#print("Points: ", len(mesh.points))
	#print("Points: ", mesh.points)
	#print("Point data: ", mesh.celldata)
	#n = len(mesh.points)

	#vertices_mesh = np.zeros( )
	#pressure = []
	#concentration = []
	#for i in range(0,n):
	#	for j in range(0,dimensions):
	#		vertices_mesh[i][j] = mesh.points[i][j]
	#	pressure.append(mesh.pointdata[i])
	#	concentration.append(mesh.pointdata[i])


	# This is the vertices for coupling and also the correct array sizes for the data. Data gets overwritten
	#print("Vertices_mesh: ", vertices_mesh)
	#print("Pressure: ", pressure)
	#print("concentration: ", concentration)
	
	### Set mesh for preCICE
	my_mesh.vertex_ids = interface.set_mesh_vertices( my_mesh.mesh_id, my_mesh.points )
	print( my_mesh.vertex_ids )

	### Get ID of data
	pressure_id = interface.get_data_id("Pressure", my_mesh.mesh_id)
	concentration_id = interface.get_data_id("Concentration", my_mesh.mesh_id)

	#dt = interface.initialize()
	fileNumber = 0

	print("preCICE initialized. Begin coupling iterations...")

	t = 0
	dt = interface.initialize()

	#while interface.is_coupling_ongoing():
	for vtu_file in vtu_file_list:

		vtu_data = meshio.read( vtu_file.path )
		print( "VTU data: {}".format(vtu_file) )
		print( "Time difference: abs(t - vtu_file.t) / t = {} ".format(abs(t - vtu_file.t) / (t +TIME_EPS )))

		assert (abs(t - vtu_file.t) / (t+TIME_EPS) < TIME_EPS), "Time does not fit with time from data set!\n    Coupling time: {}\n    Data set time: {}".format( t, vtu_data.t ) 

		pressure = extract_data_from_vtu( vtu_data, "rho" )
		concentration = extract_data_from_vtu( vtu_data, "x^tracer_0" )
#		for i in range(0,n):
#			pressure[i] = mesh.pointdata[i]
#			concentration[i] = mesh.pointdata[i]

		#print("Pressure: ", pressure)
		#print("Concentration: ", concentration)
		#print("Vertex Ids shape: ", np.shape( my_mesh.vertex_ids ) )
		#print("Vertex Ids: ", my_mesh.vertex_ids)
		print("fileNumber: ", fileNumber)
		
		interface.write_block_scalar_data(pressure_id, my_mesh.vertex_ids, pressure)
		interface.write_block_scalar_data(concentration_id, my_mesh.vertex_ids, concentration)

		dt = interface.advance(dt)
		print( "Sleeping for a {}s".format(SLEEP_TIME) )
		time.sleep(SLEEP_TIME)

		fileNumber += 1

		t += dt
		print( "Simulation time:  {}".format(t) )

	print( "No more new files to process! If the final simulation time has not been reached, you have to kill the process now." )

	interface.finalize()
	print("Closing visualization...")


def parse_args():
	parser = argparse.ArgumentParser(description=
                                     "Read meshes, partition them and write them out in internal format.")
	#parser.add_argument("in_meshname", nargs="+", help="The meshes used as input")
	parser.add_argument("pvd_filename", type=str, help="Path to ParaView PVD file to parse")
	parser.add_argument("precice_config", type=str, help="Path to preCICE XML configuration")
	parser.add_argument("--tag", "-t", dest="tag", default=None,
		help="The CellData tag for vtk meshes")
	parser.add_argument("--log", "-l", dest="logging", default="INFO", 
		choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], 
		help="Set the log level. Default is INFO")
	#parser.add_argument("--datadim", "-d", dest="datadim", default=1, type=int, help="Dimensions of the function. Default is 1 (scalar function.)")
	return parser.parse_args()

if __name__ == "__main__":
	main()


