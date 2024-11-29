#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: C. Martinez-Sebastian

GOAL: Creating a program to read all and create pkl from all given models in a 
grid. 
"""
import os
from class_FASTWIND import *
import time

try: os.chdir("/net/nas/proyectos/hots/AAA_WORKING_2030/grids/fastwind/FW10.6/iac_iacobgrid_HHeCNO_Z100/iac_iacobgrid_HHeCNO_Z100_He10_v106")
except: os.chdir("/media/cms/Elements/PhD/DATA/FASTWIND/iac_iacobgrid_HHeCNO_Z100_He10_v106")

folders = os.listdir('.')
models = []
for m in folders:
    if 'CNO' in m and m[0]=='T':
        if float(m[1:4])>=250:
            if m[:22] not in models:
                models.append(m[:22])
elems = ['C','N','O','H_He']

model_lines_created = os.listdir('../iac_iacobgrid_HHeCNO_Z100_He10_v106_pkl/original')
for model in models:  
    print('Reading model '+model)
    for el in elems:
        if model[:-3]+el+'_lines.pkl' not in model_lines_created:
            start_fin = time.time()
            A = OUT_to_pkl_lines(model, el)
            end_fin = time.time()
            print(el,end_fin - start_fin)
        else:
            print('Model '+model+' already exits')
