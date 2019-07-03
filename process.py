from __future__ import (absolute_import, division, print_function,
    unicode_literals)

import os
import json
from tqdm import tqdm_notebook as tqdm

from spheroid import spheroid

def default(o):
    if isinstance(o, np.int64): return int(o)
    raise TypeError

def _sortFiles(path):

    print('Please verify that the filename order is $xy$ then $z$ then $t$')

    for fileName in tqdm(os.listdir(path)):

        if not os.path.isdir(path + r'\\' + fileName):

            try:

                fileName, ending = fileName.split('.')

                if not 't' in fileName:

                    _, position = fileName.split('xy')
                    position, _ = position.split('z')

                    time = '00'

                elif 't' in fileName:

                    _, position = fileName.split('xy')
                    position, time = position.split('z')
                    _, time = time.split('t')


                if not os.path.exists(path + r'\\' + position):
                    os.mkdir(path + r'\\' + position)

                if not os.path.exists(path + r'\\' + position + r'\\' + time):
                    os.mkdir(path + r'\\' + position + r'\\' + time)

                os.rename(path + r'\\' + fileName, path + r'\\' + position + r'\\'
                    + time + r'\\' + fileName)

                return print('job done')

            except: print("check the file name structure")


def _saveSpheroid(sph, path):

    print(path)

    with open(path, 'w') as fp:

        json.dump(sph, fp, default = default)


def _makeSingleSpheroidClass(path, spheroidFolder, timeFolder, zRatio, rNoyau,
    dCells, channels):

    print('prep image: ' + spheroidFolder + ' folder and time ' + timeFolder)

    filePath =  path + r'\\' + spheroidFolder + r'\\' + timeFolder

    Sph = spheroid(filePath, spheroidFolder, timeFolder, zRatio, rNoyau, dCells)
    # Initialize spheroid

    try:

        if len(channels) == 2: # Improve dependancy on channel number...

            Sph._loadImage(channels[0], 'NucImage') # Load live cells
            Sph._loadImage(channels[1], 'DeadImage') # Load dead cells

        elif len(channels) == 1:

            Sph._loadImage(channels[0], 'NucImage')

        print('image made, starting nuclei ID')

        Sph._getNuclei() # identification of nuclei positions

        print('nuclei gotten, make spheroid')

        Sph._makeSpheroid() # creation of dict object

        if len(channels) == 2:

            Sph._initializeDead()

        if not os.path.exists(path + r'\\' + 'Spheroids'):
            os.mkdir(path + r'\\' + 'Spheroids')

        _saveSpheroid(Sph.Spheroid, path + r'\\' + 'Spheroids' +
            '\\spheroid_' + spheroidFolder + r'_' +  timeFolder + '.json')

        Sph._verifySegmentation()

    except: print('Error on: '+ spheroidFolder + ' folder and time ' + timeFolder)


def _makeSpheroidClass(path, zRatio, rNoyau, dCells, channels):

    """
    ====== COMMENT =======

    Function to be optimized for parrallization.
    """

    for spheroidFolder in tqdm(os.listdir(path)):

        spheroidPath = path + '/' + spheroidFolder

        if os.path.isdir(spheroidPath):

            for timeFolder in os.listdir(spheroidPath):

                timePath = spheroidPath + r'\\' + timeFolder

                if os.path.isdir(timePath):

                    _makeSingleSpheroidClass(path, spheroidPath,
                        timePath, zRatio, dCells, channels)

    return print('Spheroids made')

def _makeSphParrallel(path, key, zRatio, rNoyau, dCells, channels):

    """
    ====== COMMENT =======

    Function to be optimized for parrallization.
    """

    _, spheroidFolder = key.split(path)
    _, spheroidFolder, timeFolder = spheroidFolder.split('\\')

    _makeSingleSpheroidClass(path, spheroidFolder, timeFolder, zRatio, rNoyau,
        dCells, channels)

    return print('Spheroids made')
