# FuMuX a fake DuMuX "solver"

The python script reads in DuMuX output files (in Paraview's VTU format) and transfers it via preCICE to the visualization component. The python script mimics the true behavior of the real DuMuX solver. Slight deviations in the solution are possible due to round-off and mapping errors when the data is read in from the VTU files. 

## Important notes

- One needs appropriate Paraview output from DuMuX, i.e. a `pvd`-file and the corresponding `vtu`-files that contain the solution data. I placed the data in the directories `case1/` and `case2/`.
- There are scripts for two test cases
    - `case1`: **Single fracture** case residing in the branch `case1-precice-v1.X`.
        - One may start the `fumux` solver as: 
        ```
        fumux.py case1/case1_single_tracer_fracture.pvd case1/precice-config.xml
        ``` 
    - `case2`: Test case containing **9 fractures**. Data is transferred and received seperately on these 9 meshes. The code is in the branch `case2-precice-v1.X`.
        - One may start the `fumux` solver as: 
        ```
        ./fumux.py case2/case2_regular_tracer_fracture.pvd case2/precice-config.xml 
        ```
- The script has been tested with preCICE v1.6.X. Newer releases of preCICE (v2.X.X) are using a new API. If that is needed for testing, please let me know and I will update the scripts!

## Dependencies

- `preCICE` 1.6.X including preCICE python bindings `precice_future`
- Python 3
- Several Python 3 packages
    - `numpy`
    - `meshio`
    - `xml`