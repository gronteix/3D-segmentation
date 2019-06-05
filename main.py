#!/usr/bin/env python

# -*- coding: utf-8 -*-



# ==================================================================

# SpheroidSegment

#  A spheroid segmentation tool custom-made for image analysis
# and graph construction.

#

#   Copyright 2019 Gustave Ronteix

#   MIT License

#

#   Source Repository: https://github.com/gronteix/3D-segmentation

# ==================================================================

import multiprocessing as mp
import os

import process

###### PARAMETERS ######

path = r'D:\Gustave\27052019_TL'
zRatio = 1/4
rNoyau = 12
dCells = 80

###### START OF STUDY ######

output = mp.Queue()

processes = [mp.Process(target=_makeSpheroidClass, args=(path, key, zRatio, rNoyau, dCells))
    for key in tqdm(os.listdir(path))]

for p in tqdm(processes):
    p.start()

for p in tqdm(processes):
    p.join()
