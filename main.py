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
import glob
import tqdm

import process

###### PARAMETERS ######

livePosition = 1
deadPosition = 0

channels = [livePosition, deadPosition]

zRatio = 1/5
rNoyau = 10
dCells = 70

path = r'X:\Gustave\Experiments\Nuclei Segmentation\01072019\tif'

###### START OF STUDY ######

process._sortFiles(path)

if __name__ == '__main__':

    output = mp.Queue()

    processes = [mp.Process(target=process._makeSphParrallel, args=(path, key, zRatio,
        rNoyau, dCells, channels)) for key in glob.glob(path + r'\\**\\**')]

    for p in processes:
        p.start()

    for p in tqdm(processes):
        p.join()
