# T-matrix
For the calculation of T-matrices, with a simple MPI implementation

## Usage 
* Compile in root folder using make (please check the corresponding library locations of your environment).
* Create a scatterer using the Python script provided
* Run using the command
```
./T_matrix -m mesh.h5 -T T.h5 -a 1.d-7 -lambda 550d-9
```
or for 4-node parallel MPI run
```
mpirun -np 4 T_matrix -m mesh.h5 -T T.h5 -a 1.d-7 -lambda 550d-9
```
and so on.

### Notes
The geometry meshes are to be `tetgen`-compatible. Geometry generation routines are available, and can be run e.g. in `tetgen` or `quartet` mode. Former tends to generate unevenly sized tetrahedra inside the geometry while the latter is much more optimal on the inside, though the surface can be poor. Optimal geometry generation thus depends on the choice of initial surface refinement level, tetrahedralization refinement level and the tetrahedralization engine.
