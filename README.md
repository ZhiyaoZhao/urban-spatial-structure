# codeshare

This repository is used for maintaining the FAIR ( find-able, accessible, inter-operable, and reusable) data project

#Notice
1. To read and write shapefile，we used the Pyshp(shapefile.py) that is contributed by karimbahgat and his team. You can find the source code from https://github.com/GeospatialPython/pyshp.

2. If there is some error such as "No module named 'new_point_rank'",when you are running the code. You can add the full path of file(new_point_rank.py) to the directories list. 
## Example
    >>>import sys
    >>>sys.path
    >>>sys.path.append("..\\")#the full path of the folder

#Introduce

1.You can use shape_sim.py to calculate the shape similarity between two polygons.

2.The Self_organizing_clustering_process.py can make you to cluster the buildings by the setted similarity threshold.

3.The main orientation of polygons can be measured by the use of  shape_moment.py.

4.In the paper that we will submit to IJGIS, we proposed five VAU-based spatial metrics, three of them(HBO,CTI,and EOV) can be calculated by the spatial_metrics.py. the other spatial metrics are calcualated by the field calculator of attribute table in Arcgis10.3 software.

5.The three kinds of social function likelihood can be calculated by social_function_likelihood_formula.py 

6.The difference minimum spanning tree and the difference maximum spanning tree can be constructured by DMIST_and_DMAST.py.
