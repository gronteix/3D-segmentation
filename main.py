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

livePosition = 0
deadPosition = 1

channels = [livePosition, deadPosition]

zRatio = 1/5
rNoyau = 9
dCells = 70

path = r'PATH NAME'

###### START OF STUDY ######

process._sortFiles(path)
process._makeSpheroidClass(path, zRatio, rNoyau, dCells)

output = mp.Queue()

processes = [mp.Process(target=_makeSphParrallel, args=(path, key, zRatio,
    rNoyau, dCells)) for key in glob.glob(path + r'\\**\\**')]

for p in tqdm(processes):
    p.start()

for p in tqdm(processes):
    p.join()
