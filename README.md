## Slicer
This python script slices 3D models from STL files to produce gcode with infill and supports.

# Credits
This algorithm was written by Aaron Schaer and Michelle Ross, with help from Steven Hernandez, for Digital Fabrication.

# Usage
Our script takes exactly three inputs:

1. Name of STL file to slice
2. The desired thickness of each layer in millimeters
3. The desired percentage of infill, between 0.0 and 1.0

Our script runs best on Python 2.7.

Example usage:
`python slicing.py cube.stl 0.08 0.2`

## Github
We used git's version control to develop our algorithm. For our actual repository, please visit https://github.com/Nehri/slicing_algorithm 
