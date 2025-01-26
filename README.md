# Technical Task No. 2

## Description
A convex polygon (closed polyline) is given in 2D by the list of nodes’ coordinates.
Implement an algorithm that searches all axes of symmetry for this polygon or reports that
the polygon is non-symmetric.

## Input-output
The algorithm should be framed as a console application:
1. Input (command line options): name of text file containing coordinates of nodes;
2. Output: pairs of points defining the axis of symmetry.

### Samples


| Example 1.                                                                                                            | Example 2.                                                                 | Example 3.                                                                         |
|-----------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------|------------------------------------------------------------------------------------|
| Nodes:<br/> 0 0<br/>1 0<br/>1 1<br/>0 1<br/>Output: <br/>0 0 - 1 1<br/>1 0 - 0 1<br/> 0.5 0 - 0.5 1<br/>0 0.5 - 1 0.5 | Nodes:<br/> 0 0<br/> 2 1 <br/> 0 3 <br/> -2 1<br/> Output: <br/> 0 0 - 0 3 | Nodes:<br/>0.1 1.0 <br/>-1.0 0.0<br/>-1.0 0.0<br/>1.0 -0.5<br/>2.0 1.0<br/>Output:<br/>non-symmetric |


Requirements
1. The algorithm should be implemented in C++ 
2. All classes and methods should be precisely documented. 
3. The implementation should not involve external libraries like OpenCascade, CGAL, boost, etc.
4. Code accuracy is the main criteria for the successful result


