# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 14:48:05 2021

@author: lenovo
"""


import xml.etree.ElementTree as ET
import numpy as np
import pickle
import os

def extract_trace_grps(inkml_file_abs_path):
    trace_grps = []
    tree = ET.parse(inkml_file_abs_path)
    root = tree.getroot()
    doc_namespace = "{http://www.w3.org/2003/InkML}"
    for i in root.findall(doc_namespace + 'annotation'):
        if i.attrib['type'] == 'truth':
            latex = i.text

    latex = latex[1:]
    latex = latex[:-1]
    traceGrpWrapper = root.findall(doc_namespace + 'traceGroup')[0]
    traceGroups = traceGrpWrapper.findall(doc_namespace + 'traceGroup')
    for traceGrp in traceGroups:
        latex_class = traceGrp.findall(doc_namespace + 'annotation')[0].text
        traceViews = traceGrp.findall(doc_namespace + 'traceView')
            # Get traceid of traces that refer to latex_class extracted above
        id_traces = [traceView.get('traceDataRef') for traceView in traceViews]
        # Construct pattern object
        trace_grp = {'label': latex_class, 'traces': [], 'traces_id': [int(x) for x in id_traces]}
    
            # Find traces with referenced by latex_class
        traces = [trace for trace in root.findall(doc_namespace + 'trace') if trace.get('id') in id_traces]
    
        for idx, trace in enumerate(traces):
            coords = []
            for coord in trace.text.replace('\n', '').split(','):
                    # Remove empty strings from coord list (e.g. ['', '-238', '-91'] -> [-238', '-91'])
                coord = list(filter(None, coord.split(' ')))
                    # Unpack coordinates
                x, y = coord[:2]
                    # print('{}, {}'.format(x, y))
                if not float(x).is_integer():
                        # Count decimal places of x coordinate
                    d_places = len(x.split('.')[-1])
                        # ! Get rid of decimal places (e.g. '13.5662' -> '135662')
                        # x = float(x) * (10 ** len(x.split('.')[-1]) + 1)
                    x = float(x) * 10000
                else:
                    x = float(x)
                if not float(y).is_integer():
                        # Count decimal places of y coordinate
                    d_places = len(y.split('.')[-1])
                        # ! Get rid of decimal places (e.g. '13.5662' -> '135662')
                        # y = float(y) * (10 ** len(y.split('.')[-1]) + 1)
                    y = float(y) * 10000
                else:
                    y = float(y)
    
                    # Cast x & y coords to integer
                x, y = round(x), round(y)
                coords.append([x, y])
            trace_grp['traces'].append(coords)
        trace_grps.append(trace_grp)
    return trace_grps, latex

def get_seq(sample):
    traces_seq = {}
    for x in sample:
        traces, traces_id = x['traces'], x['traces_id']
        for i in range(len(traces)):
            traces_seq[traces_id[i]] = traces[i]
    
    traces_seq = sorted(traces_seq.items(), key=lambda x:x[0])

    seq = []
    for id_, traces in traces_seq:
        tmp = []
        for point in traces:
            tmp.append(point + [id_])
        seq += tmp
    return seq

def outlier_filter(raw_seq, dis, cos):
    res = []
    def dist(a,b):
        return ((a[0]-b[0])**2 + (a[1]-b[1])**2)**0.5
    def dot(a,b):
        return a[0]*b[0]+a[1]*b[1]
    for i in range(len(raw_seq)):
        if i == 0 or i == len(raw_seq) - 1:
            res.append(raw_seq[i])
        
        else:
            if raw_seq[i][2] == raw_seq[i+1][2] and  raw_seq[i][2] == raw_seq[i-1][2]:
                delta_x, delta_y = raw_seq[i][0] - raw_seq[i-1][0], raw_seq[i][1] - raw_seq[i-1][1]
                tmp_dist = dist(raw_seq[i],raw_seq[i-1])
                v1 = [raw_seq[i+1][0]-raw_seq[i][0], raw_seq[i+1][1]-raw_seq[i+1][1]]
                v2 = [raw_seq[i][0]-raw_seq[i-1][0], raw_seq[i][1]-raw_seq[i-1][1]]
                tmp_cos = dot(v1,v2)/(dist(raw_seq[i],raw_seq[i-1])*dist(raw_seq[i+1],raw_seq[i])+1e-10)
                
                if tmp_cos >= cos or tmp_dist <= dis:
                    #print(tmp_cos, )
                    continue
                else:
                    res.append(raw_seq[i])
            else:
                res.append(raw_seq[i])
   # print(len(res))
    return res
    

def normalization(raw_seq):
    sum_px_l = 0
    sum_py_l = 0
    sum_len = 0
    for i in range(len(raw_seq)-1):
        if raw_seq[i][2] == raw_seq[i+1][2]:
            x1, x2, y1, y2 = raw_seq[i][0], raw_seq[i+1][0], raw_seq[i][1], raw_seq[i+1][1]
            len_ = ((x2-x1)**2 + (y2-y1)**2)**0.5
            px_l = 0.5 * len_ * (x1+x2)
            py_l = 0.5 * len_ * (y1+y2)
            sum_px_l += px_l
            sum_py_l += py_l
            sum_len += len_
    mu_x, mu_y = sum_px_l/sum_len, sum_py_l/sum_len
    sum_delta = 0
    for i in range(len(raw_seq)-1):
        x1, x2, y1, y2 = raw_seq[i][0], raw_seq[i+1][0], raw_seq[i][1], raw_seq[i+1][1]
        len_ = ((x2-x1)**2 + (y2-y1)**2)**0.5
        tmp = (x2-mu_x)**2 + (x1-mu_x)**2  + (x1-mu_x)*(x2-mu_x)
        delta = 1/3 * len_ * tmp
        sum_delta += delta
    
    delta_x = (sum_delta/sum_len) ** 0.5
    for i in range(len(raw_seq)-1):
        raw_seq[i][0] = (raw_seq[i][0] - mu_x) / delta_x
        raw_seq[i][1] = (raw_seq[i][1] - mu_y) / delta_x
    raw_seq[-1][0] = (raw_seq[-1][0] - mu_x) / delta_x
    raw_seq[-1][1] = (raw_seq[-1][1] - mu_y) / delta_x
         
    return raw_seq

def feature_exactor(raw_seq):
    features =[]
    for i in range(len(raw_seq)):        
        x, y = raw_seq[i][0], raw_seq[i][1]
        
        
        if i == len(raw_seq) - 1:
            delta_x, delta_y = 0, 0
            s1, s2 = 0, 1
        else:
            x1, y1 = raw_seq[i+1][0], raw_seq[i+1][1]
            delta_x, delta_y = x1 - x, y1 - y
            s1, s2 = int(raw_seq[i][2]==raw_seq[i+1][2]), int(raw_seq[i][2]!=raw_seq[i+1][2])
        if i >= len(raw_seq) - 2:
            delta_x1, delta_y1 = 0, 0
        else:
            x2, y2 = raw_seq[i+2][0], raw_seq[i+2][1]
            delta_x1, delta_y1 = x2 - x, y2 - y
        
        features.append([x,y,delta_x,delta_y,delta_x1,delta_y1,s1,s2])
    return features
        
trace_groups, latex = extract_trace_grps('UN_101_em_15.inkml')
  
seq = get_seq(trace_groups)
seq = outlier_filter(seq, 0.0001, 1)
# normalization(seq) 
# features = feature_exactor(seq)
# with open('data_2016_features/' + 'test' + '.ascii', 'w') as fw:
#     for row in features:
#         row = str(row)
#         row = row.replace(',', '')
#         row = row.replace('[', '')
#         row = row.replace(']', '')
#         fw.write(row + '\n')
    

align = []
for i in trace_groups:
    id_ = str(i['traces_id'])
    id_ = id_.replace(',', '')
    id_ = id_.replace('[', '')
    id_ = id_.replace(']', '')
    align.append([i['label'], id_])



# with open('data_2016_align/' + 'test' + '.align', 'w') as fw:
#     for row in align:
#         row = str(row)
#         row = row.replace(',', ' ')
#         row = row.replace('\'', '')
#         row = row.replace('[', '')
#         row = row.replace(']', '')
#         fw.write(row + '\n')

# tree = ET.parse('UN_101_em_15.inkml')
# root = tree.getroot()
# doc_namespace = "{http://www.w3.org/2003/InkML}"
# traceGrpWrapper = root.findall(doc_namespace + 'traceGroup')[0]
# traceGroups = traceGrpWrapper.findall(doc_namespace + 'traceGroup')
# latex = root.findall(doc_namespace + 'annotation')[]

# tmp_caption = latex
# tmp_caption_true = []
# i = 0
# while i < len(tmp_caption): 
#     if tmp_caption[i] == '\\':  
#         j = i+1
#         for j in range(i+1, len(tmp_caption)):
#             asc = ord(tmp_caption[j])
#             if (asc >= 97 and asc <= 122) or (asc >= 65 and asc <= 90):
#                 j += 1
#             else:
#                 tmp_caption_true.append(tmp_caption[i:j])
#                 i = j
#                 break
#     elif tmp_caption[i] == ' ':
#         i += 1
#     else:
#         tmp_caption_true.append(tmp_caption[i])    
#         i += 1
# tmp_caption_true = ' '.join(tmp_caption_true)    
# with open('test.txt', 'w') as fw:
#     fw.write(tmp_caption_true)


def read_datafile(file_dir):
    caption = []
    size = 0
    for training_inkml in os.listdir(file_dir):
        size += 1
        name = training_inkml.replace('.inkml', '')
        #print(training_inkml)
        trace_groups, tmp_caption = extract_trace_grps(file_dir + training_inkml)
        
        tmp_caption_true = []
        i = 0
        while i < len(tmp_caption): 
            if tmp_caption[i] == '\\':  
                j = i+1
                while True:
                    if j == len(tmp_caption):
                        tmp_caption_true.append(tmp_caption[i:j])
                        i = j
                        break
                    asc = ord(tmp_caption[j])
                    if (asc >= 97 and asc <= 122) or (asc >= 65 and asc <= 90):
                        j += 1
                    else:
                        if j == i+1:
                            tmp_caption_true.append(tmp_caption[i:j+1])
                            i = j+1
                        else:
                            tmp_caption_true.append(tmp_caption[i:j])
                            i = j
                        break
                    
            elif tmp_caption[i] == ' ':
                i += 1
            else:
                tmp_caption_true.append(tmp_caption[i])    
                i += 1
        tmp_caption_true = ' '.join(tmp_caption_true) 
        caption.append([name, tmp_caption_true])
        
        raw_seq = get_seq(trace_groups)
        #raw_seq = outlier_filter(raw_seq, 0.0001, 1)
        normalization(raw_seq)
        features = feature_exactor(raw_seq)
            
        #train.append(trace_groups)
    
        
        align = []
        for i in trace_groups:
            id_ = str(i['traces_id'])
            id_ = id_.replace(',', '')
            id_ = id_.replace('[', '')
            id_ = id_.replace(']', '')
            align.append([i['label'], id_])
            
        with open('data_2016_features/' + name + '.ascii', 'w') as fw:
            for row in features:
                for x in row:
                    fw.write(str(x) + ' ')
                fw.write('\n')
        
        
        with open('data_2016_align/' + name + '.align', 'w') as fw:
            for row in align:
                fw.write(row[0] + '\t' + row[1] + '\n')
    
    print('Total size: ' + str(size)) 
           
    with open('test_caption2016.txt', 'w') as fw:
        for i in caption:
            fw.write(i[0]+'\t'+i[1]+'\n')

read_datafile('TEST2016_INKML_GT/')



   
        
        

