# optiMesh

OpenFOAM mesh smoothing. Place the settings for optiMesh in a system/optiMeshDict file.

See the tutorials/omesh case for an example.

# Installation

run the Allwmake script

# Other tools included

* removeCells
* collapseCells

# Alternatives to laplacian smoothing
An objective optimizer is included in optiMesh (eg. optimizing the orthogonality), but I haven't yet seen it produce better results than laplacian smoothing (and it's much slower). So this is certainly work in progress!
