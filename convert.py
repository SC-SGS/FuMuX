import meshio
import vtk

'''
mesh = meshio.read(
    "case1-1.vtu",  # string, os.PathLike, or a buffer/open file
    file_format="vtu",  # optional if filename is a path; inferred from extension
)
# mesh.points, mesh.cells, mesh.cells_dict, ...

# mesh.vtk.read() is also possible

print(mesh.cells)


#meshio.write(
#    "foo.vtk",  # str, os.PathLike, or buffer/ open file
#    mesh,
#)

meshio.write_points_cells(
    "foo.vtk",
    mesh.points,
    mesh.cells,
    file_format="vtk-ascii",
)

#mesh.vtk.write() #is also possible
'''

#m = meshio.Mesh.read("case1-1.vtu", "vtu")  # same arguments as meshio.read
#print(m.pointdata)
#m.write("foo.vtk")  # same arguments as meshio.write, besides `mesh`

reader = vtk.vtkXMLUnstructuredGridReader()
reader.SetFileName("case1-3.vtu")
reader.Update()
vtkmesh = reader.GetOutput()
print(vtkmesh.GetCellData)

writer = vtk.vtkUnstructuredGridWriter()
writer.SetFileName("data-3.vtk")
writer.SetInputData(vtkmesh)
writer.Write()



