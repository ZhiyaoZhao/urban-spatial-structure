# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:52:11 2019

@author: anonymity
"""
import math

import shapefile

#import new_point_rank
import xlrd

import numpy as np

import networkx as nx

#import gdal

import scipy

from scipy.spatial import Delaunay

#from scipy.spatial import KDTree

import Self_organizing_clustering_process as SO

import shape_moment

import cv2

import shape_sim

from pandas import Series


##################center of building footprints
                
loop_test=nx.read_shp("..\\text_region2_point.shp")

########## loading building footprints data
myshp = open("..\\test_region2.shp", "rb")

mydbf = open("..\\test_region2.dbf", "rb")

shapes = shapefile.Reader(shp=myshp, dbf=mydbf)

shapep = shapes.shapes()

recordp = shapes.records()
 
shape_list = []

for i in range(len(shapep)):

    dic = {'osmid':[],'code':[],'fclass':[],'name':[],'typee':[],'shapetype':[],'point':[],'bbox':[],'parts':[],}
    #recode操作
    dic['osmid'] = recordp[i][0]
    #dic['code'] = recordp[i][1]
    #dic['fclass'] = recordp[i][2]
    #dic['name'] = recordp[i][3]
    #dic['typee'] = recordp[i][4]
               
    #shape操作
    dic['shapetype'] = shapep[i].shapeType
    for j in range(len(shapep[i].points)):
        dic.setdefault('point',[]).append(shapep[i].points[j])
    dic['bbox'] = shapep[i].bbox
    dic['parts'] = shapep[i].parts
    
    if (len(shapep[i].points)) > 2:
               
        shape_list.append(dic);
    
for j in range(len(shape_list)):
    
    point1 = shape_list[j]['point']
    
    #print(len(point1))

    if point1[0][0] == point1[len(point1)-1][0] and point1[0][1] == point1[len(point1)-1][1]:
        
        del point1[len(point1)-1]
        
    if abs(shape_sim.vector_angle(point1[len(point1)-1][0],point1[len(point1)-1][1],point1[0][0],point1[0][1],point1[1][0],point1[1][1])*180/math.pi-180)<1:
    
        del point1[0]
       
    if abs(shape_sim.vector_angle(point1[len(point1)-2][0],point1[len(point1)-2][1],point1[len(point1)-1][0],point1[len(point1)-1][1],point1[0][0],point1[0][1])*180/math.pi-180)<1:
        
        del point1[len(point1)-1]
        
    num_of_shape_point = len(point1)-1
        
    for xx in range(1,num_of_shape_point):
        
        pointer_delete = 0 
    
        for v in range(1,len(point1)-1):
    
            v =v-pointer_delete
        
            x1 = point1[v-1][0]
            y1 = point1[v-1][1]
            x2 = point1[v][0]
            y2 = point1[v][1]
            x3 = point1[v+1][0]
            y3 = point1[v+1][1]

            angle = shape_sim.vector_angle(x1,y1,x2,y2,x3,y3)
    
            if (abs(angle*180/math.pi-180))<10:
            
                del point1[v]
            
                pointer_delete = pointer_delete + 1

######################data processing
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

node_list =list(loop_test.nodes)

loop_test_m = np.array(node_list)

loop_tri_region2 = Delaunay(loop_test_m )

connect_set = SO.search_neighbor(loop_tri_region2,node_list)

class_list_s = traverse_all_shape(connect_set)

########################## HBO ###########################

#calculation of euclidean distance between two building footprints
def calculatelength(X1,Y1,X2,Y2):
    
    s = math.sqrt(math.pow((X2 - X1),2) + math.pow((Y2 - Y1),2))#计算两点的欧氏距离
    
    return s

#weight matrix of HBO metric
def IDW_matrix(shape_set,center_list,hull_id):
    
    number_of_shape = len(shape_set[hull_id])
    
    #print(shape_set[hull_id])
    
    temp_matrix = np.zeros((number_of_shape,number_of_shape))
    
    for i_row in range(number_of_shape):
        
        for j_line in range(number_of_shape):
            
            index_row = shape_set[hull_id][i_row]

            index_line = shape_set[hull_id][j_line]
            
            X1 = center_list[index_row][0]
            X2 = center_list[index_line][0]
            Y1 = center_list[index_row][1]
            Y2 = center_list[index_line][1]
    
            temp_matrix[i_row][j_line] = calculatelength(X1,Y1,X2,Y2)
        
    inverse_matrix = np.power(temp_matrix,-1)
    
    IDW_matrix = np.zeros((number_of_shape,number_of_shape))
    
    sum_list = []
    
    for ii_row in range(number_of_shape):
        
        temp_sum = 0
        
        for jj_line in range(number_of_shape):
            
            if ii_row == jj_line:
                
                continue;
                
            else:
                
                temp_sum = temp_sum + inverse_matrix[ii_row][jj_line]
                
        sum_list.append(temp_sum)
        
        
    for iii_row in range(number_of_shape):
        
        for jjj_line in range(number_of_shape):
            
            if iii_row == jjj_line:
                
                continue;
                
            else:
                
                IDW_matrix[iii_row][jjj_line] = inverse_matrix[iii_row][jjj_line]/sum_list[iii_row]
                    
    return(IDW_matrix)

#Calculate the Angle between the main direction of building footprints and the positive direction 
def calculate_main_direction(n,shape_listt):
    
    #print(shape_list[j]['point'])
    #point_gravity = shape_moment.grivate_point(shape_listt[n]['point'])
    
    point_temp = shape_listt[n]['point']
        
    point_temp_np = np.array(point_temp, dtype = int)

    rect = cv2.minAreaRect(point_temp_np)
    
    angle = math.radians(rect[2])    
    
    width = rect[1][0]
        
    height = rect[1][1]
    
    if width >= height:
        
        angle = (-1)*angle
        
    else:
            
        angle = (-1)*angle+(0.5*math.pi)
        

    angle = 180 - angle*180/math.pi
        
    return(angle)

def Moran_index(shape_set,center_list,hull_id,shape_listt):
    
    fenmu = 0
    
    fenzi = 0
    
    w = IDW_matrix(shape_set,center_list,hull_id)
    
    number_of_shape = len(shape_set[hull_id])
    
    average_angle = 0
    
    ang = []
    
    for i_shape in range(number_of_shape):
        
        n = shape_set[hull_id][i_shape]
        
        ang.append(calculate_main_direction(n,shape_listt))
    
    average_angle = sum(ang)/len(ang)
    
    #print(ang,average_angle)
    
    for i_row in range(number_of_shape):
        
        for j_line in range(number_of_shape):
            
            if i_row == j_line:
                
                continue;
            
            else:
                
                fenzi = fenzi + w[i_row][j_line]*(abs(ang[i_row]-average_angle))*(abs(ang[j_line]-average_angle))
                
        
        fenmu = fenmu + (ang[i_row]-average_angle)*(ang[i_row]-average_angle)
        
    if fenmu == 0 :
        
        moran_i = 1
        
    else:
        
        moran_i = fenzi/fenmu
        
    return(moran_i)   
    
#the resulted 'moran_i' is the HBO of i-th convex hull  
    
########################## CTI ###########################
    
hull_polygon,component_set = SO.hull_of_pattern(class_list_s)

#the ratio of the sum area of involved buildings to the area of VAU
def area_building_hull_ratio(hull_polygon1,component_set1):
    
    hull_area = []
    
    ratio = []
    
    for i_hull in range(len(hull_polygon1)):
        
        moment_hull = shape_moment.moment_shape(hull_polygon1[i_hull])
        
        area_hull = abs(moment_hull[0])
        
        hull_area.append(area_hull)
        
        sum_area = 0
        
        for j_building in range(len(component_set1[i_hull])):
            
            building_id = component_set1[i_hull][j_building]
            
            moment_building_area = shape_moment.moment_shape(shape_list[building_id]['point'])
            
            building_area = abs(moment_building_area[0])
            
            sum_area = sum_area + building_area
            
        ratio.append(sum_area/area_hull)

    return(hull_area,ratio)
    
#the resulted 'ratio' is the CTI
    
########################## EOV ###################################

#The sample data including average building area of the convex hull 
data = xlrd.open_workbook("..\\Average_building_area_of_each_hull.xlsx")

table = data.sheets()[0]

#row
nrows = table.nrows

#col
ncols = table.ncols

def read_cell_of_excel(table,nrows):
    
    list_of_nonrec = []
    
    #length of the list is equal to the number of street blocks. in this study,we set the length to 371
    hull_area = [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],
                [],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
    
    print(len(hull_area))
    
    for row in range(1,nrows):
        
        region = int(table.cell_value(row,1))
        
        #print(region)
        
        mean_area = table.cell_value(row,4)
        
        if mean_area == 0:
            
            continue;
            
        else:
        
            hull_area[region].append(table.cell_value(row,2))
        
    return(hull_area)
    
#calcualteing the entropy of a series of data
def ent(data):
    p_data= data.value_counts()/len(data) # calculates the probabilities
    
    entropy=scipy.stats.entropy(p_data)  # input probabilities to get the entropy 
    return entropy

#计算图形在设定完成的采样间隔下的斜率信息熵
def entropy_of_hull_area(hull_area):
    
    entropy_list = []
    
    for i in range(len(hull_area)):
        
        ser = Series(hull_area[i])
        
        entropy = ent(ser)
        
        #print(slope)
        entropy_list.append(entropy)
        
    return(entropy_list);
    
    