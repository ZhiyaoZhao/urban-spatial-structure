# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 11:56:59 2019

@author: anonymity
"""
from scipy.spatial import Delaunay

import numpy as np

import matplotlib.pyplot as plt

import networkx as nx

import xlrd

import xlwt

import math

import shapefile

import Self_organizing_clustering_process as SO

####################loading center of urban blocks data

loop_test=nx.read_shp('..\\street_block_point.shp')

#####################load social function likelihood
data = xlrd.open_workbook('..\\social function ikelihood.xls')

table = data.sheets()[0]

nrows = table.nrows

ncols = table.ncols

####################data processing

center_point_list =list(loop_test.nodes)


#center_point_list = construe_Delaunay_tri(shape_list)
loop_test_m = np.array(center_point_list)
loop_tri_region2 = Delaunay(loop_test_m )


all_the_triangle = loop_tri_region2.simplices

list_of_node_of_edge_linear =  SO.tri_into_node_of_edge(all_the_triangle)

list_of_node_of_edge_linear.sort()

list_of_node_of_edge_linear = SO.delet_repeat(list_of_node_of_edge_linear)

########################################

def vector_angle(x1,y1,x2,y2,x3,y3):
    
    vector_a = [x2-x1,y2-y1]
    
    vector_b = [x3-x2,y3-y2]
    
    
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
    
    return(angle);

def construe_Delaunay_tri(ploint_list):
    
    gravity_point_list = []
    
    for i in range(len(ploint_list)):
        
        gravity_point_list.append(gravity_list(ploint_list[i]['point']))
        
    return(gravity_point_list)

def gravity_list(point1):
    
    sum_x =0
    
    row = len(point1)
    
    column = len(point1)
    
    for i in range(row):
        
        sum_x = sum_x + point1[i][0]
    
    point_x = sum_x/(row)
    
    sum_y =0
    
    for j in range(column):
        
        sum_y = sum_y + point1[j][1]
    
    point_y = sum_y/(column)
    
    point_gravity = [0,0]
    
    point_gravity[0] = point_x
    
    point_gravity[1] = point_y

    return(point_gravity);

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

def delet_repeat(list_of_node_of_edge_loop):
    
    k = 0
    
    for i in range(len(list_of_node_of_edge_loop)-1):
        
        if list_of_node_of_edge_loop[i-k] == list_of_node_of_edge_loop[i+1-k]:
            
            del list_of_node_of_edge_loop[i-k]
            
            k = k + 1
            
    return(list_of_node_of_edge_loop);

#for your own data, please adjust the row number !!!!!!
def read_likelihood_from_excel(excel_data):
    
    ind_list = []
    com_list = []
    resid_list = []
    
    for row in range(1,nrows):
        #the row number is depend on the real excel
        industry_temp = table.cell_value(row,3)
        ind_list.append(industry_temp)

        residential_temp = table.cell_value(row,4)
        resid_list.append(residential_temp)
        
        commercial_temp = table.cell_value(row,5)
        com_list.append(commercial_temp)
    
    return(ind_list,com_list,resid_list)

#generating the wight of the edges
def generate_weighted_edges(list_of_node_of_edge_linear_list,likelihood):
    
    edge_list = []
    
    for i_edge in range(len(list_of_node_of_edge_linear_list)):
        
        likelihood_temp = abs(likelihood[list_of_node_of_edge_linear_list[i_edge][0]] - likelihood[list_of_node_of_edge_linear_list[i_edge][1]])
        
        temp_tuple = (list_of_node_of_edge_linear_list[i_edge][0],list_of_node_of_edge_linear_list[i_edge][1],likelihood_temp)
        
        edge_list.append(temp_tuple)
        
    return(edge_list);

def calculatelength(X1,Y1,X2,Y2):
    
    s = math.sqrt(math.pow((X2 - X1),2) + math.pow((Y2 - Y1),2))
    
    return (s)

def normal_list(dis):
    
    temp = []
    
    for i in range(len(dis)):
        
        temp.append((max(dis)+min(dis)) - dis[i])
    
    normal_temp = []
        
    for j in range(len(temp)):
        
        normal_temp.append(temp[j]/sum(temp))
    
    return(normal_temp)

def normal_list_new(dis):
    
    normal_temp = []
        
    for j in range(len(dis)):
        
        normal_temp.append(dis[j]/sum(dis))
    
    return(normal_temp)
    
#recoreding all the edge length of adjacent edges 
def adjacent_edge_length(adjacent_edge):
    
    length_list = []
    
    for i in range(len(adjacent_edge)):
    
        x1 = center_point_list[adjacent_edge[i][0]][0]
        y1 = center_point_list[adjacent_edge[i][0]][1]
        x2 = center_point_list[adjacent_edge[i][1]][0]
        y2 = center_point_list[adjacent_edge[i][1]][1]
        
        length_list.append(calculatelength(x1,y1,x2,y2))
        
    return(length_list)


#calculating the difference of likelihood between connected blocks
def local_difference(tree_edge,likelihood):
    
    value_set = []
    
    local_difference_list = []
    
    x1_list = []
    
    # traverse all the edges of the tree
    for i in range(len(tree_edge)):
        
        #print(i)
        temp_u = tree_edge[i][0]

        temp_v = tree_edge[i][1]
        
        temp_adjacent_edge = []
        
        temp_value = []
        
        for j in range(len(tree_edge)):
            
            if j == i:
                
                continue;
            
            if temp_u in tree_edge[j] or temp_v in tree_edge[j]:
                
                temp_adjacent_edge.append(tree_edge[j])
                
                temp_value.append(abs(likelihood[tree_edge[j][0]]-likelihood[tree_edge[j][1]]))
        
        value_set.append(temp_value)
        
        temp_wight = normal_list_new(temp_value)
        
        #Xi
        Xi = abs(likelihood[tree_edge[i][0]]-likelihood[tree_edge[i][1]])
        #print(Xi)
        
        x1_list.append(Xi)
        
        #X_mean
        local_mean_value = (sum(temp_value) + Xi)/(len(temp_value)+1)
        #print(local_mean_value)
        
        temp_sum = 0
        
        for ij in range(len(temp_wight)):
            
            temp_sum = temp_sum + temp_wight[ij] * (temp_value[ij] - local_mean_value)
            
        local_difference_up = temp_sum * (Xi-local_mean_value)
        
        temp_sum_j = 0
        
        for ji in range(len(temp_value)):
            
            temp_sum_j = temp_sum_j + pow((temp_value[ji] - local_mean_value),2)
        
        local_difference_i_down = temp_sum_j / len(temp_value) - pow(local_mean_value,2)
        
        local_difference = local_difference_up/local_difference_i_down
        
        if local_difference >1 :
            
            local_difference = 1
            
        if local_difference <-1:
        
            local_difference = -1
        
        local_difference_list.append(local_difference) 
        
    return(x1_list)


############################################
#preserve the edge of spanning tree into 
def plot_line(linear_list):
    
    w = shapefile.Writer()
    
    aa = []
    
    for i_number_line in range(len(linear_list)):
        
        if len(linear_list[i_number_line])<=2:
            
            continue;
            
        else:
    
            aa.append(1)
    
    w.field("field1")
            
    
    
    for i in range(len(linear_list)):
        
        #print(i)
            
        line = []
            
        w.record(i)
            
        line_point_list = []
            
        for j in range(len(linear_list[i])):
                
            point = [0,0]
                
            point[0] = center_point_list[(linear_list[i][j])][0]
                
            point[1] = center_point_list[(linear_list[i][j])][1]
                
            line_point_list.append(point)
                
        line.append(line_point_list)
            
        w.line(parts=line)
    
    #the save path for shape file   
    w.save("..\\commercial_MAST")
    #w.save("..\\residential_MIST")
    #w.save(".\\industry_MAST")
    
    return('true');
    
    
def write_excel(list1):
    
    #x_list,y_list,y_interpolate_list = direction_distance_descriptor(shape_list)
    y_interpolate_list = list1
    
    #创建workbook
    workbook = xlwt.Workbook(encoding = 'ascii')
    
    #创建表
    worksheet = workbook.add_sheet('My Worksheet')
    
    ii = 1
    
    nn = len(y_interpolate_list)
    
    print(nn)
    
    worksheet.write(0, 0, label = 'fid')
        
    worksheet.write(0, 1, label = 'local_moran')
    
    for jj in range(nn):
        
        osmid = y_interpolate_list[jj]
            
        worksheet.write(ii, 1, label = osmid)
        
        worksheet.write(ii, 0, label = jj)
            
        ii =ii + 1
    
    # the wight of tree edges
    workbook.save("..\\commercial_in_DMAST.xls")
    #workbook.save("\\residential_in_DMIST.xls")
    #workbook.save("\\industry_in_DMAST.xls")
    
    return('true')
    
############################################
#reading three kinds of social function likelihood from excel
industrial,commercial,residential = read_likelihood_from_excel(data)

#calculating the edge weight
industrial_network = generate_weighted_edges(list_of_node_of_edge_linear,industrial)
#initialization of industrial graph
ind_graph = nx.Graph()
#generating the graph with weight
ind_graph.add_weighted_edges_from(industrial_network)

residential_network = generate_weighted_edges(list_of_node_of_edge_linear,residential)

res_graph = nx.Graph()

res_graph.add_weighted_edges_from(residential_network)


commercial_network = generate_weighted_edges(list_of_node_of_edge_linear,commercial)

com_graph = nx.Graph()

com_graph.add_weighted_edges_from(commercial_network)


ind_DMIST = nx.minimum_spanning_tree(ind_graph)

ind_DMAST = nx.maximum_spanning_tree(ind_graph)


res_DMIST = nx.minimum_spanning_tree(res_graph)

res_DMAST = nx.maximum_spanning_tree(res_graph)


com_DMIST = nx.minimum_spanning_tree(com_graph)

com_DMAST = nx.maximum_spanning_tree(com_graph)

##############################################
#extracting the edge and  the likelihood difference of spanning tree
def capture_edge_of_tree(tree,correspondind_likelihood):
    
    temp_edges = tree.edges
    
    edges = list(temp_edges)
    
    #print(type(sss1[0]))

    #ssnode = np.array(center_point_list)
    
    #print(len(sss1))
    
    ##this part code can be used to show your result
    '''
    plt.figure(figsize=(10,10))
    for simplex in edges:
    #print(simplex)
    #print(center_point_list[simplex[0]])
    #print('+++++')
    #print(center_point_list[simplex[1]])
        plt.plot([center_point_list[simplex[0]][0],center_point_list[simplex[1]][0]], [center_point_list[simplex[0]][1],center_point_list[simplex[1]][1]],c = 'r')
    plt.show()
    '''
    difference_list = local_difference(edges,residential)
    
    return(edges,difference_list)
    
    
    