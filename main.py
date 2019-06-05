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
from spheroid import spheroid

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



def _makeSpheroidClass(path, spheroidFolder, zRatio, rNoyau, dCells):

        spheroidPath = path + '/' + spheroidFolder

        if os.path.isdir(spheroidPath):

            for timeFolder in os.listdir(spheroidPath):

                timePath = spheroidPath + r'\\' + timeFolder

                if os.path.isdir(timePath):

                    print('prep image: ' + spheroidFolder + ' folder and time ' + timeFolder)

                    Sph = spheroid(timePath, spheroidFolder, timeFolder, zRatio, rNoyau, dCells)

                    try:

                        # Pay attention to the channel order in the experiments

                        Sph._loadImageNuclei(1)
                        Sph._loadImageDead(0)
                        Sph._getNuclei()
                        Sph._makeSpheroid()

                        Sph._initializeDead()

                        with open(path + '\spheroid_' + spheroidFolder + r'_' +  timeFolder + '.json', 'w') as fp:

                            json.dump(Sph.Spheroid, fp, default = default)

                        Sph._verifySegmentation()

                    except: print('Error on: '+ spheroidFolder + ' folder and time ' + timeFolder)
