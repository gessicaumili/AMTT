# AMTT
Automated Markland's Test Tool
The code refers to an automatic procedure based on the GIS environment working principles and developed it in Matlab language. 
Main discontinuity sets orientation and relative friction angles, along with slope and aspect data representing the rockface orientation of the considered outcrop, are the input data. 
The slope and aspect data are in GeoTIFF format. 
The Matlab code performs Markland's tests for planar and wedge sliding and flexural toppling, considering all the possible sets or intersections of sets.
The outputs are a series of GeoTIFF raster files describing the result for each kinematism separately and globally, which can be imported directly into GIS, with the same extent and georeferencing of the input data.
The global results can be also used to map source areas for 3D rockfall numerical simulations. 

How to cite the code:
Taboni, B.; Tagliaferri, I.D.; Umili, G. A Tool for Performing Automatic Kinematic Analysis on Rock Outcrops. Geosciences 2022, 12, 435. https://doi.org/10.3390/geosciences12120435 
