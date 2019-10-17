# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 14:48:07 2019

@author: anonymity  
"""
import shapefile

import shape_sim

import shape_moment

import numpy as np

import networkx as nx

import gdal

import math

from scipy.spatial import Delaunay

from scipy.spatial import KDTree

from scipy.spatial import ConvexHull

import matplotlib.pyplot as plt

from scipy.spatial import Voronoi, voronoi_plot_2d


#getting a full path of shapefile
myshp = open("..\\test_region2.shp", "rb")

mydbf = open("..\\test_region2.dbf", "rb")

shapes = shapefile.Reader(shp=myshp, dbf=mydbf)

shapep = shapes.shapes()

recordp = shapes.records()
 
shape_list = []

for i in range(len(shapep)):

    dic = {'osmid':[],'code':[],'fclass':[],'name':[],'typee':[],'shapetype':[],'point':[],'bbox':[],'parts':[],}
    #recode
    dic['osmid'] = recordp[i][0]

    #shape
    dic['shapetype'] = shapep[i].shapeType
    for j in range(len(shapep[i].points)):
        dic.setdefault('point',[]).append(shapep[i].points[j])
    dic['bbox'] = shapep[i].bbox
    dic['parts'] = shapep[i].parts
    
    if (len(shapep[i].points)) > 2:
               
        shape_list.append(dic);
            
for j in range(len(shape_list)):
    
    point1 = shape_list[j]['point']
    
    delet_point_index_list = []
    
    part_list = list(shapep[j].parts)
    
    part_list.append(len(point1)-1)
    
    delet_count = 0
    #print(shapep[j].parts)

#getting a full path of buildind center point shapefile 
loop_test=nx.read_shp("..\\text_region2_point.shp")

#identifying all the adjacent relationship of every building footprints
def search_neighbor(triangulation_net,shape):
    
    shape_number= len(shape)
    
    #recording all the triangles of the triangulation_net
    all_triangle = triangulation_net.simplices
    
    #recording all the edges of the triangulation_net
    edge_linear =  tri_into_node_of_edge(all_triangle)
    
    #print(len(edge_linear))
    
    edge_linear.sort()
    
    edge_linear = delet_repeat(edge_linear)
    
    #print(len(edge_linear))
            
    edge_linear = sorted(edge_linear,key=lambda x: (x[0], -x[1]))
            
    edge_linear = delet_repeat(edge_linear)
    
    #print(len(edge_linear))
    
    #initializing the set of connections
    connect_set = np.zeros((shape_number,1))
    
    connect_set = connect_set.tolist()
    
    #print(edge_linear)
    
    for i_tuple in edge_linear:
        
        #print(i_tuple[0])
        
        if isinstance(connect_set[(i_tuple[0])][0],float):
            
            #temp = []
            
            #temp.append(i_tuple[1])
            
            del connect_set[i_tuple[0]][0]
            
            connect_set[i_tuple[0]].append(i_tuple[1])
            
        elif isinstance(connect_set[(i_tuple[0])][0],int):
            
            connect_set[(i_tuple[0])].append(i_tuple[1])
            
        if isinstance(connect_set[(i_tuple[1])][0],float):
            
            #temp = []
            
            #temp.append(i_tuple[1])
            
            del connect_set[i_tuple[1]][0]
            
            connect_set[i_tuple[1]].append(i_tuple[0])
            
        elif isinstance(connect_set[(i_tuple[1])][0],int):
            
            connect_set[(i_tuple[1])].append(i_tuple[0])
            
    for i_connect in range(len(connect_set)):
        
        connect_set[i_connect] = list(set(connect_set[i_connect]))
        
    return(connect_set);

#finding clustered buildings
def find_pattern(present_index,head_list,leaf_list,jumplist):
    
    present_leaf = leaf_list.copy()
    
    for i_component in present_leaf:
        
        #print(i_component)
        
        if i_component in jumplist:
            
            continue;
            
        print(present_index,i_component)
        
        #In this study, the similarity threshold was set to 0.9,
        if shape_sim.SN(shape_list[present_index]['point'],shape_list[i_component]['point'])>0.9:
            
            moment_present_index = shape_moment.moment_shape(shape_list[present_index]['point'])
            
            moment_i_component = shape_moment.moment_shape(shape_list[i_component]['point'])
            
            area_present_index = abs(moment_present_index[0])
            
            area_i_component = abs(moment_i_component[0])
            
            if area_present_index >area_i_component:
                
                area_present_index,area_i_component = area_i_component,area_present_index
            
            if abs((area_present_index/area_i_component)-1)<0.3:
                
                center1 = shape_moment.grivate_point(shape_list[present_index]['point'])
                
                center2 = shape_moment.grivate_point(shape_list[i_component]['point'])
                
                #if the distance between two buildings is more than 200m,they will not be identified as adjacent builidngs 
                if calculatelength(center1[0],center1[1],center2[0],center2[1]) <200:
                    
                    #print('@@@',i_component)
            
                    head_list.append(i_component)
            
                    #leaf_list.extend(connect_set[i_component])
            
                    leaf_list = list(set(leaf_list).union(set(connect_set[i_component])))
            
                    if present_index in leaf_list:
            
                        leaf_list.remove(present_index)
                
                    if i_component in leaf_list:
            
                        leaf_list.remove(i_component)
            
                    jumplist.append(i_component)
            
    return(head_list,leaf_list,jumplist);

#this function is used to traverse all building footprints including the scenario. meanwhile, all clustered building footprint will be found.
def traverse_all_shape(connec_set):
    
    #this list is used for storing shape index of compared building footprints 
    jump_list = []

    i_present= 0
    
    class_list = []
    
    while len(jump_list) < len(connec_set):
        
        #print('##')
        
        if i_present in jump_list:
            
            i_present += 1
            
            continue;
            
        else :
                
            #print(i_present)
        
            head_list = []
        
            leaf_list = []
        
            head_list.append(i_present)
        
            leaf_list.extend(connec_set[i_present])
        
            list_new = []
        
            list_old = []
        
            list_new.extend(connec_set[i_present])
            
            jump_list.append(i_present)
        
            #If no new leaves are added, the loop stops
            while (len(list_new)-len(list_old)) > 0:
            
                list_old = leaf_list.copy()
                
                #print(list_old)
                
                present_index = i_present
            
                head_list,leaf_list,jump_list = find_pattern(present_index,head_list,leaf_list,jump_list)
            
                #print(head_list,leaf_list,jump_list)
            
                list_new = leaf_list.copy()
                
                #print(list_new)
                
                #print(list_old)
            
            i_present += 1;
            
            class_list.append(head_list)
    
    return(class_list)

#the triangulation network is dispersed into vertex of each  triangle 
def tri_into_node_of_edge(A):
    
    list_of_node_of_edge = []
    
    for i in range(len(A)):
        
        temp = A[i,:]
        
        temp_list = temp.tolist()
        
        edge_1 = [temp_list[0],temp_list[1]]
        edge_1 = tuple(edge_1)
        edge_2 = [temp_list[0],temp_list[2]]
        edge_2 = tuple(edge_2)
        edge_3 = [temp_list[1],temp_list[2]]
        edge_3 = tuple(edge_3)
        list_of_node_of_edge.append(edge_1)
        list_of_node_of_edge.append(edge_2)
        list_of_node_of_edge.append(edge_3)
        
    return(list_of_node_of_edge);

#the repetitive 
def delet_repeat(list_of_node_of_edge_loop):
    
    k = 0
    
    for i in range(len(list_of_node_of_edge_loop)-1):
        
        if list_of_node_of_edge_loop[i-k] == list_of_node_of_edge_loop[i+1-k]:
            
            del list_of_node_of_edge_loop[i-k]
            
            k = k + 1
            
    return(list_of_node_of_edge_loop);

def vector_angle(x1,y1,x2,y2,x3,y3):
    #origial vector
    vector_a = [x2-x1,y2-y1]
    #new vector
    vector_b = [x3-x2,y3-y2]
    
    #calculating the module value of vector
    mo_vector_a = math.sqrt(math.pow(vector_a[0],2)+math.pow(vector_a[1],2))
    mo_vector_b = math.sqrt(math.pow(vector_b[0],2)+math.pow(vector_b[1],2))
    
    if mo_vector_a*mo_vector_b == 0:
        
        return(0)
    
    cosin_angle = ((vector_a[0]*vector_b[0]+vector_a[1]*vector_b[1]))/(mo_vector_a*mo_vector_b)
    
    if cosin_angle>1:
        
        cosin_angle = 1
        
    if cosin_angle<-1:
        
        cosin_angle = -1
    
    angle = math.acos(cosin_angle)
    
    #return radian type
    return(angle);

def return_temp_adjacent(index_list,indicator,list_of_node_of_edge_linear):
    
    #the list 'temp' stores all the connections between the current index builidng to adjacent buildings, the connection is recorded in the form of a tuple（index，near）
    temp = []
        
    for i in range(len(list_of_node_of_edge_linear)):
                
        if index_list[indicator] in list_of_node_of_edge_linear[i]:
                    
            if index_list[indicator] == list_of_node_of_edge_linear[i][1]:
                        
                replace = (list_of_node_of_edge_linear[i][1],list_of_node_of_edge_linear[i][0])
                        
                temp.append(replace)
                        
            else:
                
                temp.append(list_of_node_of_edge_linear[i])
    #deleting the adjacent and same elements      
    temp1 = delet_repeat(temp)
            
    b = sorted(temp1,key=lambda x: (x[0], -x[1]))
            
    temp1 = delet_repeat(b)
    
    return(temp1);

def calculatelength(X1,Y1,X2,Y2):
    
    s = math.sqrt(math.pow((X2 - X1),2) + math.pow((Y2 - Y1),2))
    
    return s

def calculate_adjacent_distance(temp1):
    
    adjacent_distance = []
    
    for i_adjacent in range(len(temp1)):
        
        temp_distance_index = [0,0]
        
        X1 = (loop_test_m[(temp1[i_adjacent][0])])[0]
                        
        Y1 = (loop_test_m[(temp1[i_adjacent][0])])[1]
                        
        X2 = (loop_test_m[(temp1[i_adjacent][1])])[0]
                        
        Y2 = (loop_test_m[(temp1[i_adjacent][1])])[1]
        
        temp_distance_index = [calculatelength(X1,Y1,X2,Y2),temp1[i_adjacent][1]]
        
        adjacent_distance.append(temp_distance_index)
        
    adjacent_distance.sort()
    
    return(adjacent_distance);

def print_deplicate_shape(shape_lists):
    
    lists = []
    
    for i_shape in range(len(shape_lists)):
        
        if shape_list[i_shape]['point'] in lists:
            
            print(i_shape)
        
        lists.append(shape_list[i_shape]['point'])
        
    return(lists)
    
##Finding the vertices of all graphics in each pattern, saving the vertices, and finding its convex hull
def hull_of_pattern(class_list_si):
    
    hull_point_list = []
    
    id_of_component = []
    
    for i_pattern in class_list_si:
        
        if len(i_pattern)<2:
            
            continue;
            
        else:
            
            id_of_component.append(i_pattern)
            
            temp_point_list = []
            
            for j_shape in i_pattern:
                
                temp_point_list.extend(shape_list[j_shape]['point'])
                
        hull_point_list.append(temp_point_list)
    
    #vertexs of each clusted building footprints are used to generate 
    ConvexHull_list = []
    
    for i_patter_hull in hull_point_list:
        
        vertex_point = []
        
        i_patter_hull = np.array(i_patter_hull)
        
        temp_hull = ConvexHull(i_patter_hull)
        
        temp_hull_point = temp_hull.vertices
        
        #temp_hull_point = sp.spatial.qhull.Delaunay(temp_hull).convex_hull
        
        temp_hull_point = temp_hull_point.tolist()
        
        for j_vertex in temp_hull_point:
            
            vertex_point.append(i_patter_hull[j_vertex].tolist())
            
        #vertex_point.append(vertex_point[0])
        
        ConvexHull_list.append(vertex_point)
    
    
    return(ConvexHull_list,id_of_component)

#preserve these generated minimum convex hull into shapefile
def plot_polygon(hull_polygoni):
    
    w = shapefile.Writer()
      
    w.field("field2")
            
    for i_hull in range(len(hull_polygoni)):
        
        if len(hull_polygoni[i_hull])<=2:
            
            continue;
            
        else:

            w.record(i_hull)
            
            polygon_point_list = []
 
            polygon_point_list.append(hull_polygoni[i_hull])
            
            #print(polygon_point_list)
            
            w.poly(parts=polygon_point_list)
    #save path of the generated VAU        
    w.save("..\\VAU_of_test_region2.shp")
    
    return('true');
    
def calculatelength(X1,Y1,X2,Y2):
    
    s = math.sqrt(math.pow((X2 - X1),2) + math.pow((Y2 - Y1),2))#计算两点的欧氏距离
    
    return s
