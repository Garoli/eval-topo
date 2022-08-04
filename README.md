#  Quantification of topological errors due to image transformations

**Context** : Third year project \
**Year** : 2021-2022

**Students** :

* Sarah Brood *sarah.brood@ecole.ensicaen.fr*
* Heithem Dridi *heithem.dridi@ecole.ensicaen.fr*

**Tutors**:

* Yukiko Kenmochi *yukiko.kenmochi@unicaen.fr*
* Sebastien Fourey *sebastien.fourey@unicaen.fr*

Description
-----
Sources and examples to compute different topological characteristics, that aims to quantify the topolgical errors due to image transformations: 

Tools available for 2 and 3 dimensions: 
- Transformations
- Interpolation
- Computation of betti-numbers
- Persistent homology

Display tools : 
- Boards 2D for topological characteristics 
- 3D viewer for objects 
- plot of the persistent homology

*More informations can be found in the report in the document folder*

Installation & Usage
-----
Be sure to install : 
* [DGtal](https://dgtal.org/download/) 
with following dependencies : ```-DWITH_GMP=true -DWITH_CAIRO=true -DWITH_QGLVIEWER=true -DWITH_MAGICK=true ```
* [Gudhi](https://gudhi.inria.fr/doc/latest/installation.html) with CGAL


Examples can be found in the source folder : 
- evalTranspoTrans2D for transformation, interpolation, topological characteristics, and homological persistent for 2D 
- evalTranspoTrans3D for transformation, interpolation, topological characteristics, and homological persistent for 3D
- viewer3D to use to get an application viewer for 3D objects

```shell
mkdir build && cd build
cmake ..
make

./evalTopoTrans2D -h
This file is part of the project : Topological evaluation of image transformations
The aim of this file is to show an example of use of our project for 2D
Usage: ./evalTopoTrans2D [OPTIONS] 1 [2] [3] [4]

Positionals:
  1 TEXT:FILE REQUIRED                  input file in .pgm / .ppm
  2 TEXT                                output file in .pgm / .ppm. Default to output/output.pgm
  3 TEXT:{NN,BIL,BIC}                   set interpolation: NN: Nearest Neighbour, BIL: Bilinear, BIC: Bicubic
  4 BOOLEAN                             if true compute homology. Default to true

Options:
  -h,--help                             Print this help message and exit
  -i,--input TEXT:FILE REQUIRED         input file in .pgm / .ppm
  -o,--output TEXT                      output file in .pgm / .ppm. Default to output/output.pgm
  -I,--interpolation TEXT:{NN,BIL,BIC}  set interpolation: NN: Nearest Neighbour, BIL: Bilinear, BIC: Bicubic
  -H,--homology BOOLEAN                 if true compute homology. Default to true

```

File format that work fine are PGM for 2D and .vol for 3D.  
Doxygen documentation is available.

Sources
----
Documentation and courses we used are available in the directory named "Documents/documentations".

Images sources :
- PGM files : [PGMA Files](https://people.sc.fsu.edu/~jburkardt/data/pgma/pgma.html )
- Volumes : [VolGallery](https://github.com/dcoeurjo/VolGallery)
