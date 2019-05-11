### Create a file with coordinates and main eigen vectors ###
* Open write_fabForVTK.py
* Adapt the path to ground_truth.txt
* Run `$ python write_fabForVTK.py`
* A file "fabForVTK.dat" is created.

### Convert fabForVTK.dat to a VTK file ###
* Compile TensorImageView
```
$ cd TensorImageView
$ mkdir build
$ cd build
$ cmake ..
$ make
```
* Convert fabForVTK.dat to a VTK file
```
$ \path\to\TensorImageView \path\to\F16_R_stance_3p0_mask.mhd fabForVTK.dat ground_truth_fab.vtk
```
### Visualize with Paraview ###
* Download Paraview from [https://www.paraview.org/download/](https://www.paraview.org/download/)
* Open Paraview `$ paraview`
* Open images in paraview: 
  - LowRes image
    - Files > Open > F16_R_stance_3p0_mask.mhd > Apply
    - In the drop-down menu above the image view, change the view from `Outline` to `Slice`
    - In the tab `Properties` (the lower left panel), change from `XY plane` to `XZ plane`
    - Change from slice `15` to `7`
  - Fabric image
    - Files > Open > ground_truth_fab.vtk > Apply
    - In the menu bar, select `Clip` 
    - In the tab `Properties`
        - Plane parameters: select `Y Normal` with y value at 20, click `Apply`
        - Coloring: edit color to black
