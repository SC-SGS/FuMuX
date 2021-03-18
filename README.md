# FuMuX (=fake DuMuX)

The python script reads in DuMuX output files (in Paraview's VTU format) and transfers it via [preCICE](https://www.precice.org/) to the visualization component [inpreCICE](https://github.com/SteScheller/inpreCICE/). The python script mimics the behavior of a real [DuMuX](https://dumux.org/) solver. Slight deviations in the solution are possible due to round-off and mapping errors when the data is read in from the VTU files instead of accessing the data directly from DuMuX. This script is meant to simplify the development process of [inpreCICE](https://github.com/SteScheller/inpreCICE/).

Currently supported data sets originate from the publication *Verification benchmarks forsingle-phase flow in three-dimensional fracturedporous media*. You can more information in the [Call for Participation](https://arxiv.org/pdf/1809.06926.pdf) and the [final paper](https://doi.org/10.1016/j.advwatres.2020.103759). The supported data sets are *case 1* and *case 2* of the verification benchmarks.

## Important notes / How to

- One needs appropriate Paraview output from DuMuX, i.e. a `pvd`-file and the corresponding `vtu`-files that contain the solution data. I recommend to place the data in the subdirectories `case1/` and `case2/`.
- Datasets for case 1  (`case1-data.tar.gz `) and case 2 (`case2-data.tar.gz `) can be found on the DaRUS of the University of Stuttgart: [https://doi.org/10.18419/darus-1773](https://doi.org/10.18419/darus-1773)
- Unpack the data set into the corresponding subdirectory `case1/` and `case2/`. The directories already contain a appropriate `precice-config.xml` for the corresponding test case.
- There are scripts for two test cases residing in *branches*. Currently, full functionality for the two supported cases is **only** provided by these branches:
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
- FuMuX has been tested with preCICE v1.6.X. Newer releases of preCICE (v2.X.X) are using a new API and would require an update of inpreCICE as well.

## Dependencies

- [preCICE](https://www.precice.org/) 1.6.X including preCICE python bindings `precice_future`
- Python 3
- The following Python 3 packages:
    - `numpy`
    - `meshio`
    - `xml`
- [inpreCICE](https://github.com/SteScheller/inpreCICE/) is not a direct dependency for running FuMuX, but nothing will happen/FuMuX will hang if the visualization tool is not running at the same time.

## Acknowledgements

FuMuX has been developed in the [SFB1313](https://www.sfb1313.uni-stuttgart.de/) which has been funded by the [DFG](https://www.dfg.de/) under research grant [327154368](https://gepris.dfg.de/gepris/projekt/327154368).
## License

This work is licensed under the [GPLv3](./LICENSE).