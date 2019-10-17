# codeshare
this repository is used for maintaining the FAIR ( find-able, accessible, inter-operable, and reusable) data project

#Notice
1. TO read and write shapefile，we used the Pyshp(shapefile.py) that is contributed by karimbahgat and his team. You can fine the source code from https://github.com/GeospatialPython/pyshp.
2. If there is some error such as "No module named 'new_point_rank'",when you are running the code. You can add the full path of file(new_point_rank.py) to the directories list. 
>>>import sys
>>>sys.path
>>>sys.path.append("..\\")#the full path of the folder

#Introduce
1.You can use shape_sim.py to calculate the shape similarity between two polygons.
2.The Self_organizing_clustering_process.py can make you to cluster the buildings by the setted similarity threshold.
