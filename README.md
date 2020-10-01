# FuMuX a fake DuMuX "solver"

The python script reads in DuMuX output files (in Paraview's VTU format) and transfers it via preCICE to the visualization component. The python script mimics the true behavior of the real DuMuX solver. Slight deviations in the solution are possible due to round-off and mapping errors when the data is read in from the VTU files. 

## Important notes

- There are scripts for two test cases
    - `case1`: Single fracture case residing in the branch `case1-precice-v1.X`
        - One may start the `fumux` solver as: `python3 fumux.py case1/case1_single_tracer_fracture.pvd case1/precice-config.xml` for example.
    - `case2`: Test case containing 9 fractures. Data is transferred and received seperately on these 9 meshes.
- The script has been tested with preCICE v1.6.X. Newer releases of preCICE (v2.X.X) are using a new API. If that is needed for testing, please let me know and I will update the scripts!
