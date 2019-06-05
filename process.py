
import os
import json
from spheroid import spheroid

def _makeSpheroidClass(path, spheroidFolder, zRatio, rNoyau, dCells):

        spheroidPath = path + '/' + spheroidFolder

        if os.path.isdir(spheroidPath):

            for timeFolder in os.listdir(spheroidPath):

                timePath = spheroidPath + r'\\' + timeFolder

                if os.path.isdir(timePath):

                    print('prep image: ' + spheroidFolder + ' folder and time '
                        + timeFolder)

                    Sph = spheroid(timePath, spheroidFolder, timeFolder, zRatio,
                        rNoyau, dCells)

                    try:

                        # Pay attention to the channel order in the experiments

                        Sph._loadImageNuclei(1)
                        Sph._loadImageDead(0)
                        Sph._getNuclei()
                        Sph._makeSpheroid()

                        Sph._initializeDead()

                        with open(path + '\spheroid_' + spheroidFolder + r'_' +
                            timeFolder + '.json', 'w') as fp:

                            json.dump(Sph.Spheroid, fp, default = default)

                        Sph._verifySegmentation()

                    except: print('Error on: '+ spheroidFolder +
                        ' folder and time ' + timeFolder)
