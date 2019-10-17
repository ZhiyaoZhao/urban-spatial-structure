# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 10:50:52 2019

@author: Zhao Zhiyao
"""

import cv2

import numpy as np

import shape_moment

import matplotlib.pyplot as plt

import numpy as np

import math

import new_point_rank

def include_angle(x1,y1,x2,y2):
    
    vector_a = [x2-x1,y2-y1]

    vector_east = [1,0]

    mo_vector_a = math.sqrt(math.pow(vector_a[0],2)+math.pow(vector_a[1],2))
    mo_vector_east = math.sqrt(math.pow(vector_east[0],2)+math.pow(vector_east[1],2))
    
    #print(mo_vector_a)
    
    if mo_vector_a == 0:
        
        mo_vector_a = 0.001
    
    cosin_angle = ((vector_a[0]*vector_east[0]+vector_a[1]*vector_east[1]))/(mo_vector_a*mo_vector_east)
    
    angle = math.acos(cosin_angle)
    
    if vector_a[1] < 0:
        angle = math.pi * 2 - angle
    

    return(angle);

def vector_angle(x1,y1,x2,y2,x3,y3):

    vector_a = [x2-x1,y2-y1]

    vector_b = [x3-x2,y3-y2]
    #vector_b = [x2-x3,y2-y3]
    

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

def tangent_angle(x1,y1,x2,y2):
    
    #original vector
    vector_a = [x2-x1,y2-y1]
    #new vector
    vector_b = [1,0]

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
    
    if (y2-y1) < 0:
        
        angle = (2*math.pi) - angle
    
    #return radian type result
    return(angle);

#Determining whether there is a clockwise and counterclockwise transition at the Angle
def judge_turningorient(x1,y1,x2,y2,x3,y3,x4,y4):
    
    vector_judge = [x4-x3,y4-y3]
    
    vector_ref   = [x2-x1,y2-y1]
    
    vector_middle= [x3-x2,y3-y2]
    
    mo_vector_judge  = math.sqrt(math.pow(vector_judge[0],2)+math.pow(vector_judge[1],2))
    #print(mo_vector_judge,1)
    mo_vector_ref    = math.sqrt(math.pow(vector_ref[0],2)+math.pow(vector_ref[1],2))
    #print(mo_vector_ref,2)
    mo_vector_middle = math.sqrt(math.pow(vector_middle[0],2)+math.pow(vector_middle[1],2))
    #print(mo_vector_middle,3)
    
    cosin_angle_j_m = ((vector_judge[0]*vector_middle[0]+vector_judge[1]*vector_middle[1]))/(mo_vector_judge*mo_vector_middle)
    
    #print(cosin_angle_j_m)
    
    cosin_angle_j_m = round(cosin_angle_j_m,4)
    
    angle_j_m = math.acos(cosin_angle_j_m)
    
    cosin_angle_j_r = ((vector_judge[0]*vector_ref[0]+vector_judge[1]*vector_ref[1]))/(mo_vector_judge*mo_vector_ref)
    
    #print(cosin_angle_j_r)
    
    cosin_angle_j_r = round(cosin_angle_j_r,4)
    
    angle_j_r = math.acos(cosin_angle_j_r)
    
    if angle_j_m > angle_j_r:
        
        return(1);
    else:
        
        return(0);

def calculatelength(X1,Y1,X2,Y2):
    
    s = math.sqrt(math.pow((X2 - X1),2) + math.pow((Y2 - Y1),2))#计算两点的欧氏距离
    
    return(s)

#calicualting the perimeter of the polygon
def Z(points):
    
    num = len(points)
    
    Z = 0
    
    for k in range(len(points) - 1):
        
        X1 = points[k][0]
        X2 = points[k+1][0]
        Y1 = points[k][1]
        Y2 = points[k+1][1]
        
        Z = Z + calculatelength(X1,Y1,X2,Y2)
    
    FX1 = points[0][0]
    FX2 = points[num-1][0]
    FY1 = points[0][1]
    FY2 = points[num-1][1]
    
    final = Z + calculatelength(FX1,FY1,FX2,FY2)
    
    return (final);

#calculating the route distance between the k-th and original point
def Sk(points,k):
    
    if k == 0:
        
        return(0);
        
    else :
    
        Sk = 0
    
        for j in range(k):
        
            X1 = points[j][0]
            X2 = points[j+1][0]
            Y1 = points[j][1]
            Y2 = points[j+1][1]
        
            Sk = Sk + calculatelength(X1,Y1,X2,Y2)
    
        return Sk;
    
#generating turning——function from vertexs of polygon   
def  turning_function(points):
    
    #get perimeter of polygon
    length = Z(points)
    
    arc_length = []
    
    for i in range(len(points)-1):
        
        arcl = Sk(points,i+1)
        
        arc_length.append(arcl)
    
    arc_length.append(length)
    
    #intializing the strip
    ratio = []
    
    #Calculating the ratio of the edge to the perimeter
    for j in range(len(arc_length)-1):
        
        portion = arc_length[j]/length
        
        ratio.append(portion)
        
    ratio.append(1)
    
    #initializing the list of angle
    clockwise_angle_list = []
    
    angle_list = []
    
    manber = 0
    
    judge_number = 0
    
    angle_menber = 0
    
    if len(points) == 3:
        
        for k in range(len(points)):
            
            if k == 0:
                
                angle_menber = angle_menber + include_angle(points[0][0],points[0][1],points[1][0],points[1][1])
                
                angle_list.append(angle_menber)
                
            elif k ==(len(points)-1):
                
                angle_menber = angle_menber + vector_angle(points[k-1][0],points[k-1][1],points[k][0],points[k][1],points[0][0],points[0][1])
                
                angle_list.append(angle_menber)
                
            else :

                angle_menber = angle_menber + vector_angle(points[0][0],points[0][1],points[1][0],points[1][1],points[2][0],points[2][1])

                angle_list.append(angle_menber)
                
    else :
        
        for k in range(len(points)):
            
            if k == 0:
                
                angle_menber = angle_menber + include_angle(points[0][0],points[0][1],points[1][0],points[1][1])
                
                angle_list.append(angle_menber)
                
            elif k ==1:
                
                judge = rotate_matrix(points[0][0],points[0][1],points[k][0],points[k][1],points[k+1][0],points[k+1][1])
                
                angle_menber = angle_menber + judge*vector_angle(points[0][0],points[0][1],points[k][0],points[k][1],points[k+1][0],points[k+1][1])
                
                angle_list.append(angle_menber)
                
            elif k ==(len(points)-1):
                
                
                judge = rotate_matrix(points[len(points)-2][0],points[len(points)-2][1],points[len(points)-1][0],points[len(points)-1][1],points[0][0],points[0][1])
                
                angle_menber = angle_menber + judge*vector_angle(points[k-1][0],points[k-1][1],points[k][0],points[k][1],points[0][0],points[0][1])
                
                angle_list.append(angle_menber)
                
            else :
                    
                judge = rotate_matrix(points[k-1][0],points[k-1][1],points[k][0],points[k][1],points[k+1][0],points[k+1][1])
                
                angle_menber = angle_menber + judge*vector_angle(points[k-1][0],points[k-1][1],points[k][0],points[k][1],points[k+1][0],points[k+1][1])
                
                angle_list.append(angle_menber)
   
    turning_f = {'strip':ratio,'angle':angle_list}

    
    return(turning_f);

def integration_f(turnning_f):
    
    length = len(turnning_f['strip'])
    
    #print(length)
    
    sum_jifen = 0
    
    for i in range(length):
        if i ==0:
            sum_jifen = sum_jifen + turnning_f['strip'][i]*turnning_f['angle'][i]
        else:
            sum_jifen = sum_jifen + (turnning_f['strip'][i]-turnning_f['strip'][i-1])*turnning_f['angle'][i]
    return(sum_jifen);


def aerfa(turnning_f,turnning_g):
    integra_f = integration_f(turnning_f)
    integra_g = integration_f(turnning_g)
    
    interim_4 = integra_g - integra_f
    
    return(interim_4);

def new_list(strip1,strip2):
    
    #strip1.pop()
    
    newl = []
    
    new1 = strip1 + strip2
    
    #print(strip1)
    #print(strip2)
    
    new1.sort()
    
    #print(newl)
    
    #newlist =strip1
    
    contral_number = 0
    
    for i in range(len(new1)-1):
        
        if new1[i - contral_number] == new1[i - contral_number + 1]:
            
            del new1[i - contral_number + 1]
            
            contral_number = contral_number + 1
            
    return(new1);

#according the strip and angle stored in TF, the angle can be look up by the index
def return_turnningfunction(x,turnning_f):
    
    return_angle = 0

    for j in range(len(turnning_f['strip'])):
        
            if x <= turnning_f['strip'][j]:
            
                return_angle = turnning_f['angle'][j]
            
                break
    return(return_angle);

def area_ratio_of_shape_with_MAR(shape_point1):
    
    frame_ID_list_np = np.array(shape_point1, dtype = float)
        
    frame_ID_list_np =frame_ID_list_np*100
    
    #print(int(frame_ID_list_np[0][0]/1000000))
    
    longitude_three = int(frame_ID_list_np[0][0]/1000000)
    
    #print(int(frame_ID_list_np[0][1]/1000000))
    
    latitude_three = int(frame_ID_list_np[0][1]/1000000)
    
    #print(type(frame_ID_list_np))
    
    #print(frame_ID_list_np.shape)
    
    longitude = frame_ID_list_np[:,0]
    
    latitude = frame_ID_list_np[:,1]
    
    longitude = longitude - longitude_three*1000000
    
    latitude = latitude - latitude_three*1000000

    #print(longitude,latitude)
    
    frame_ID_list_np = np.vstack((longitude,latitude))
    
    frame_ID_list_np = frame_ID_list_np.T
    
    #frame_ID_list_np[0] = frame_ID_list_np[0] -(int(frame_ID_list_np[0][0]/10000000))*10000000
    
    #print(frame_ID_list_np)
        
    frame_ID_list_np = frame_ID_list_np.astype(np.int64)
        
    #print(frame_ID_list_np)

    rect = cv2.minAreaRect(frame_ID_list_np)
    
    #print(rect[0],rect[1],rect[2])

    points1 = cv2.boxPoints(rect)
    
    #print(points.tolist())
        
    points1 = points1/100
    
    points1[:,0] = points1[:,0] + longitude_three*10000
    
    points1[:,1] = points1[:,1] + latitude_three*10000
        
    points1 =points1.tolist()
    
    moment_MAR = shape_moment.moment_shape(points1)
    
    #print(points,moment_MAR)
    
    area_MAR = moment_MAR[0]
    
    moment_shape = shape_moment.moment_shape(shape_point1)
    
    area_shape = moment_shape[0]
    
    hull_vertices = HULL_shape(shape_point1)
    
    moment_hull = shape_moment.moment_shape(hull_vertices)
    
    area_hull = moment_hull[0]
    
    if area_MAR ==0:
        
        RATIO = 0
        
    else:
    
        RATIO = area_shape/area_MAR
        
        #ratio_1 = area_shape/area_hull
    
    if RATIO < 0:
        
        RATIO = (-1*RATIO)
        
    #if ratio_1 < 0:
        
        #ratio_1 = (-1*ratio_1)
        
    #print(area_shape,area_MAR)
    
    #return(RATIO,ratio_1,points,hull_vertices);
    
    return(RATIO,points1,hull_vertices);

#calculating shape similarity between shapeA and shapeB	
def SN(point_A,point_B):
    
    D_A_B = (distance6(point_A,point_B)+distance6(point_B,point_A))/2
    
    area1,area2 = envelope_area(point_A,point_B)
    
    SN_A_B = 1 - D_A_B/((area1+area2)/2)
    
    return(SN_A_B)

def plot_turning_function(turning_function):
    
    angle = turning_function['angle']
    
    strip = turning_function['strip']
    
    print(angle[0])
    
    temp_angle = []
    temp_strip = []
    
    for i_angle in range(len(angle)):
        
        temp_angle.append(angle[i_angle])
        temp_angle.append(angle[i_angle])
    
    for i_strip in range(len(strip)):
        
        if i_strip == 0:
        
            temp_strip.append(0)
            temp_strip.append(strip[i_strip])
            
        else :
            
            temp_strip.append(strip[i_strip-1])
            temp_strip.append(strip[i_strip])
    
    #print(len(temp_strip),len(temp_angle))

    x = np.array(temp_strip)
    
    y = np.array(temp_angle)
    
    plt.plot(x, y, marker="*", linewidth=3, linestyle="--", color="orange")
    
    plt.legend(["Y","Z"], loc="upper right")
    
    plt.grid(True)
    
    plt.show()
    
    return(0);

def rotate_matrix(x1,y1,x2,y2,x3,y3):
    
    rotate_angle = include_angle(x1,y1,x2,y2)
    
    #print(rotate_angle)
    
    new_x3 = (x3-x2)*math.cos(rotate_angle) + (y3-y2)*math.sin(rotate_angle)
    
    #print(new_x3)
    
    new_y3 = (y3-y2)*math.cos(rotate_angle) - (x3-x2)*math.sin(rotate_angle)
    
    #print(new_y3)
    
    if new_y3 < 0:
        
        return(-1)
    
    else :
        
        return(1)
		
def HULL_shape(point): 
    
    d = np.array(point)
    
    hull = ConvexHull(d)
        
    list_hull = (hull.simplices).tolist()
    
    list_index = hull.vertices
    
    hull_vertices_list = []
    
    for i_vertices in list_index:
        
        temp = [0,0]
        
        temp[0] = point[i_vertices][0]
        
        temp[1] = point[i_vertices][1]
        
        hull_vertices_list.append(temp)
        
    ranked_hull = []
        
    for i_list_hull in range(len(list_hull)):
    
        temp_sort = [0,0]
        element1 = list_hull[i_list_hull][0]
        element2 = list_hull[i_list_hull][1]
    
        if element1>element2:
        
            element1,element2=element2,element1
        
        temp_sort[0] = element1
    
        temp_sort[1] = element2
    
        ranked_hull.append(temp_sort)
        
    return(hull_vertices_list);

def calculate_SN(shape_list):
    
    result = []
    
    for i_shape in range(len(shape_list)):
        
        result_1 = []
        
        point1 = new_point_rank.re_rank(shape_list[i_shape]['point'])
        
        for j_shape in range(len(shape_list)):
            
            point2 = new_point_rank.re_rank(shape_list[j_shape]['point'])
            
            temp = SN(point1,point2)
            
            result_1.append(temp)
            
        result.append(result_1)
        
    return(result)

#calculating the shape distance by turning function
def distance6(point1,point2):
    
    point1 = new_point_rank.re_rank(point1)
    
    #t_turnning_f = turning_function(point1)
    
    #variance_t = t_turnning_f['strip']
    
    point2 = new_point_rank.re_rank(point2)
    
    result = []
    
    for t in range(len(point1)):
        
        temp = []
        
        if t == 0:
            
            temp = point1
            
        elif t == len(point1)-1:
            
            temp.append(point1[len(point1)-1])
            
            temp.extend(point1[0:t])
            
        else: 
            
            temp.extend(point1[t:len(point1)])
        
            temp.extend(point1[0:t])
            
        #print(temp)
        
        turnning_f = turning_function(temp)
        
        #plot_turning_function(turnning_f)
        
        turnning_g = turning_function(point2)
        
        #plot_turning_function(turnning_g)
        
        #print(turnning_f,turnning_g)
        
        ds = new_list(turnning_f['strip'],turnning_g['strip'])
        
        #print(ds)
        
        normal_fx = []
        
        normal_gx = []
        
        normal_fx_temp = []
        
        normal_gx_temp = []
        
        for i in range(len(ds)-1):
            
            vari = ds[i]
            
            fx = return_turnningfunction(vari,turnning_f)
            
            #print(fx)
            
            normal_fx_temp.append(fx)

            gx = return_turnningfunction(vari,turnning_g)
            
            #print(gx)
            
            normal_gx_temp.append(gx)
        
        temp_min_fx = min(normal_fx_temp)
        
        temp_min_gx = min(normal_gx_temp)
        
        #print(temp_min_fx,temp_min_gx)
        
        for i_fx in range(len(normal_fx_temp)):
            
            #print(normal_fx_temp[i_fx])
            
            #print(normal_fx_temp[i_fx]-temp_min_fx)
            
            normal_fx.append(normal_fx_temp[i_fx]-temp_min_fx)
            
        #print(normal_fx_temp[0],normal_fx_temp[1])
            
        for i_gx in range(len(normal_gx_temp)):
            
            #print(normal_gx_temp[i_gx]-temp_min_gx)
            
            normal_gx.append(normal_gx_temp[i_gx]-temp_min_gx)
            
        #print(normal_fx)
        #print(normal_gx)
        #print(ds)
        
        sumi = 0
        
        for i_interpolation in range(len(ds)-1):
            
            
            if i_interpolation == 0:
                
                sumi = sumi + abs(ds[i_interpolation] * (normal_fx[i_interpolation] - normal_gx[i_interpolation]))
                
                #print(normal_fx[i_interpolation] - normal_gx[i_interpolation])
            
            else :
                
                sumi = sumi + abs((ds[i_interpolation]-ds[i_interpolation-1])*(normal_fx[i_interpolation] - normal_gx[i_interpolation]))
                
                #print(normal_fx[i_interpolation] - normal_gx[i_interpolation])
        
        result.append(sumi)

    return(min(result));

#calculating the envelop area of polygonA and polygonB
def envelope_area(point1,point2):
    
    point1 = new_point_rank.re_rank(point1)
    
    point2 = new_point_rank.re_rank(point2)
    
    turnning_f = turning_function(point1)
        
    #plot_turning_function(turnning_f)
        
    turnning_g = turning_function(point2)
    
    ds = new_list(turnning_f['strip'],turnning_g['strip'])
        
    #print(ds)
        
    normal_fx = []
        
    normal_gx = []
        
    normal_fx_temp = []
        
    normal_gx_temp = []
    
    for i in range(len(ds)-1):
            
        vari = ds[i]
            
        fx = return_turnningfunction(vari,turnning_f)
            
        #print(fx)
            
        normal_fx_temp.append(fx)

        gx = return_turnningfunction(vari,turnning_g)
            
        #print(gx)
            
        normal_gx_temp.append(gx)
        
    temp_min_fx = min(normal_fx_temp)
        
    temp_min_gx = min(normal_gx_temp)
        
    #print(temp_min_fx,temp_min_gx)
        
    for i_fx in range(len(normal_fx_temp)):
            
        normal_fx.append(normal_fx_temp[i_fx]-temp_min_fx)
            
    for i_gx in range(len(normal_gx_temp)):
            
        normal_gx.append(normal_gx_temp[i_gx]-temp_min_gx)
    
    A_envelope = 0
    
    B_envelope = 0
    
    for i_area in range(len(ds)-1):
        
        #print(ds[i_area+1],ds[i_area],normal_gx[i_area])
        
        if i_area == 0:
        
            A_envelope = A_envelope + (ds[i_area])*normal_fx[i_area]
        
            B_envelope = B_envelope + (ds[i_area])*normal_gx[i_area]
            
        else :
            
            A_envelope = A_envelope + (ds[i_area]-ds[i_area-1])*normal_fx[i_area]
        
            B_envelope = B_envelope + (ds[i_area]-ds[i_area-1])*normal_gx[i_area]
    
    return(A_envelope,B_envelope)