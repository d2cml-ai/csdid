# -*- coding: utf-8 -*-
"""
Created on Wed Jun  7 20:50:27 2023

@author: Carlos Guevara
"""
import os
import pickle

folder_path = r"E:\Google drive\Repositorio\csdid\did\aggte" 
os.chdir(folder_path)

from aggte import *


pickle_file = "att_gt_object.pkl"
with open(pickle_file, "rb") as f:
    out = pickle.load(f)
    
out =  aggte(out)

