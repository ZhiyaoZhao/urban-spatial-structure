# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 16:50:58 2019

@author: anonymity
"""
import xlrd

import shapefile

import xlwt

data_street_block = xlrd.open_workbook("..\\spatial metrics.xlsx")

table_street_block = data_street_block.sheets()[0]

#读取excel的总行数
nrows_street_block = table_street_block.nrows

#读取excel的总列数
ncols_street_block = table_street_block.ncols


#COM
COV1 = table_street_block.col_values(4,1,nrows_street_block)
#BD
BD1 = table_street_block.col_values(5,1,nrows_street_block)
#RI
RI1 = table_street_block.col_values(6,1,nrows_street_block)
#EOM
EOV1 = table_street_block.col_values(7,1,nrows_street_block)
#COB
COB1 = table_street_block.col_values(8,1,nrows_street_block)
#CI
CI1 = table_street_block.col_values(9,1,nrows_street_block)
#HOB
HOB1 = table_street_block.col_values(10,1,nrows_street_block)
#BS
BS1 = table_street_block.col_values(11,1,nrows_street_block)


#########redidential social function likelihood
#def residential_index(ri_col,doc_col,bd_col):
def residential_index(COV,RI,HOB,BD,EOV):
    
    max_COV = max(COV)
    min_COV = min(COV)
    
    max_RI = max(RI)
    min_RI = min(RI)
    
    max_HOB = max(HOB)
    min_HOB = min(HOB)
    
    max_BD = max(BD)
    min_BD = min(BD)
    
    max_EOV = max(EOV)
    min_EOV = min(EOV)
    
    list_COV = []
    list_RI = []
    list_HOB = []
    list_BD = []
    list_EOV = []
    
    for i_block in range(len(COV)):
        
        list_COV.append((COV[i_block]- min_COV)/(max_COV - min_COV))
        list_RI.append((RI[i_block] - min_RI)/(max_RI - min_RI))
        list_HOB.append((HOB[i_block] - min_HOB)/(max_HOB - min_HOB))
        list_BD.append((BD[i_block] - min_BD)/(max_BD - min_BD))
        list_EOV.append((EOV[i_block] - min_EOV)/(max_EOV - min_EOV))
        
    #calcluating the weight of corresponding metric
    
    sort_COV = COV.copy()
    sort_COV.sort()
    sort_RI = RI.copy()
    sort_RI.sort()
    sort_HOB = HOB.copy()
    sort_HOB.sort()
    sort_BD = BD.copy()
    sort_BD.sort()
    sort_EOV = EOV.copy()
    sort_EOV.sort()
    
    residential_index_value = []
    
    for i_w in range(len(COV)):
        
        W1 = sort_COV.index(COV[i_w])
        W2 = sort_RI.index(RI[i_w])
        W3 = sort_HOB.index(HOB[i_w])
        W4 = sort_BD.index(BD[i_w])
        W5 = sort_EOV.index(EOV[i_w])
        
        if (W1+W2+W3+W4+W5) == 0:
            
            residential_index_value.append(0)
            
        else :
        
            w1 = W1/(W1+W2+W3+W4+W5)
            w2 = W2/(W1+W2+W3+W4+W5)
            w3 = W3/(W1+W2+W3+W4+W5)
            w4 = W3/(W1+W2+W3+W4+W5)
            w5 = W3/(W1+W2+W3+W4+W5)
        
            residential_index_value.append(w1*list_COV[i_w] + w2*list_RI[i_w] + w3*list_HOB[i_w] + w4*list_BD[i_w]- w5*list_EOV[i_w])
    
        
    return(residential_index_value);
    

#########industrial social function likelihood
#calcualting industrial likelihood
def industry_index(COV,HOB,COB,BS,BD):
    
    max_COV = max(COV)
    min_COV = min(COV)
    
    max_HOB = max(HOB)
    min_HOB = min(HOB)
    
    max_COB = max(COB)
    min_COB = min(COB)
    
    max_BS = max(BS)
    min_BS = min(BS)
    
    max_BD = max(BD)
    min_BD = min(BD)
    
    list_COV = []
    list_HOB = []
    list_COB = []
    list_BS = []
    list_BD = []
    
    for i_block in range(len(COV)):
        
        list_COV.append((COV[i_block]- min_COV)/(max_COV - min_COV))
        list_HOB.append((HOB[i_block]- min_HOB)/(max_HOB - min_HOB))
        list_COB.append((COB[i_block]- min_COB)/(max_COB - min_COB))
        list_BS.append((BS[i_block]- min_BS)/(max_BS - min_BS))
        list_BD.append((BD[i_block]- min_BD)/(max_BD - min_BD))
        
    
    sort_COV = COV.copy()
    sort_COV.sort()
    sort_HOB = HOB.copy()
    sort_HOB.sort()
    sort_COB = COB.copy()
    sort_COB.sort()
    sort_BS = BS.copy()
    sort_BS.sort()
    sort_BD = BD.copy()
    sort_BD.sort()
    
    industry_index_value = []
    
    for i_w in range(len(COV)):
        
        W1 = sort_COV.index(COM[i_w])+1
        W2 = sort_HOB.index(HOB[i_w])+1
        W3 = sort_COB.index(COB[i_w])+1
        W4 = sort_BS.index(BS[i_w])+1
        W5 = sort_BD.index(BD[i_w])+1
        
        
        if (W1+W2+W3+W4+W5) == 0:
            
            industry_index_value.append(0)
            
        else :
        
            w1 = W1/(W1+W2+W3+W4+W5)
            w2 = W2/(W1+W2+W3+W4+W5)
            w3 = W3/(W1+W2+W3+W4+W5)
            w4 = W4/(W1+W2+W3+W4+W5)
            w5 = W4/(W1+W2+W3+W4+W5)

            industry_index_value.append(w1*list_COV[i_w] + w2*list_HOB[i_w] + w3*list_COB[i_w] + w4*list_BS[i_w]- w5*list_BD[i_w])
            
    return(industry_index_value);
    
#########commercial social function likelihood
def commercial_index(EOV,CI,COB,BS,RI):
    
    max_EOV = max(EOV)
    min_EOV = min(EOV)
    
    max_CI = max(CI)
    min_CI = min(CI)
    
    max_COB = max(COB)
    min_COB = min(COB)
    
    max_BS = max(BS)
    min_BS = min(BS)
    
    max_RI = max(RI)
    min_RI = min(RI)
    
    list_EOV  = []
    list_CI = []
    list_COB = []
    list_BS = []
    list_RI = []
    
    for i_block in range(len(EOV)):
        
        list_EOV.append((EOV[i_block]- min_EOV )/(max_EOV  - min_EOV ))
        list_CI.append(1 - (CI[i_block]- min_CI)/(max_CI - min_CI))
        list_COB.append((COB[i_block]- min_COB)/(max_COB - min_COB))
        list_BS.append((BS[i_block]- min_BS)/(max_BS - min_BS))
        list_RI.append((RI[i_block]- min_RI)/(max_RI - min_RI))
    
    sort_EOV  = EOV.copy()
    sort_EOV.sort()
    sort_CI = CI.copy()
    sort_CI.sort()
    sort_COB = COB.copy()
    sort_COB.sort()
    sort_BS = BS.copy()
    sort_BS.sort()
    sort_RI = RI.copy()
    sort_RI.sort()
    
    commercial_index_value = []
    
    for i_w in range(len(EOV)):
        
        W1 = sort_EOV.index(EOV[i_w])+1
        W2 = sort_CI.index(CI[i_w])+1
        W3 = sort_COB.index(COB[i_w])+1
        W4 = sort_BS.index(BS[i_w])+1
        W5 = sort_RI.index(RI[i_w])+1
        
        if (W1+W2+W3+W4+W5) == 0:
            
            commercial_index_value.append(0)
            
        else :
        
            w1 = W1/(W1+W2+W3+W4+W5)
            w2 = W2/(W1+W2+W3+W4+W5)
            w3 = W3/(W1+W2+W3+W4+W5)
            w4 = W4/(W1+W2+W3+W4+W5)
            w5 = W4/(W1+W2+W3+W4+W5)
        
            commercial_index_value.append(w1*list_EOV[i_w] + w2*list_CI[i_w] + w3*list_COB[i_w] + w4*list_BS[i_w] - w5*list_RI[i_w])
    
    return(commercial_index_value);

def write_LIKELIHOOD_list_excel(component_set):
    
    y_interpolate_list = component_set
    
    workbook = xlwt.Workbook(encoding = 'ascii')
    
    worksheet = workbook.add_sheet('My Worksheet')
    
    ii = 0
    
    nn = len(y_interpolate_list)
    
    print(nn)
    
    for jj in range(nn):
                
        worksheet.write(ii, 0, label = jj)
        
        worksheet.write(ii, 1, label = component_set[jj])
            
        ii =ii + 1
        
    #workbook.save("..\\residential_likelihood.xls")
    #workbook.save("..\\industrial_likelihood.xls")
    workbook.save("..\\commercial_likelihood.xls")
    
    return('true')