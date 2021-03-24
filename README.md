# lignumvtk.py: 

Initial proposition to use vtk library
to visualize Lignum with ParaView. The 'lignumvtk'
program creates 'vtp' xml file that can be imported
to ParaView.

lignumvtk is a python3 program that requires numpy 
and vtk python packages. These and the ParaView 
can be installed with MacPorts. Another
possibility is to create python3 virtual environment
and install nympy and vtk with 'pip'.

## Usage:

LUKE160201L:vtk$ python lignumvtk.py -h

Usage: lignumvtk.py [options]

Options:
 + -h, --help        show this help message and exit
 + -i F1, --fxml=F1  Read Lignum xml file
 + -o F2, --fvtp=F2  VTP output file
 + -c, --cylinder    Use segment base radius as segment top radius (pure cylinder)
LUKE160201L:vtk$ 
   
# coloredlines.py:

Demonstrates the use of vtk to visualize three
cylindrical segments with tube filter.



