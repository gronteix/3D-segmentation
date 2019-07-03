from __future__ import (absolute_import, division, print_function,
    unicode_literals)

import os
import numpy as np
import pandas
import os
import skimage
from scipy import ndimage, signal
from skimage import exposure
import matplotlib.pyplot as plt
from skimage import io
import scipy
from scipy import signal
from skimage import feature
import json
from matplotlib_scalebar.scalebar import ScaleBar
import networkx as nx
import timeit


class spheroid:

    """ Spheroid class containing the necessary info to build the spheroid.

    ====== COMMENT ======

     - All variables starting with a capital letter refer to class variables.
     - All variables starting with '_' refer to a dict


    ====== PARAMETERS ======

    path: string object, path to the folder storing images
    position: string object, well ID
    time: string object, time of experiment"""

    def __init__(self, path, position, time, zRatio, rNoyau, dCells):

        self.Path = path
        self.Position = position
        self.Time = time
        self.ZRatio = zRatio
        self.RNoyau = rNoyau
        self.DCells = dCells
        self.NucImage = []
        self.DeadImage = []
        self.BorderCrop = 300 # pixels cropped on border
        self.Thresh = 300 # thresh for dead cell detection
        self.ThreshCell = 150 # thresh for live cell detection
        self.Percentile = 30 # you dump pixels below this relative threshold


    def _loadImage(self, channel, type):

        """ Function to load the images corresponding to a given channel

        ====== COMMENT ======

        The function needs to be improved so as to add new channels without
        requiring to manually add new channels by hand.

        """

        image_list = []
        for filename in sorted(os.listdir(self.Path)): #assuming tif

            if '.tif' in filename:

                im = io.imread(self.Path + '/' + filename)
                image_list.append(im[channel])

        if type == 'NucImage':

            self.NucImage = np.reshape(image_list, (len(image_list),
                np.shape(image_list[0])[0], np.shape(image_list[0])[1]))[:,
                self.BorderCrop:-self.BorderCrop, self.BorderCrop:-self.BorderCrop]

        if type == 'DeadImage':

            self.DeadImage = np.reshape(image_list, (len(image_list),
                np.shape(image_list[0])[0], np.shape(image_list[0])[1]))[:,
                self.BorderCrop:-self.BorderCrop, self.BorderCrop:-self.BorderCrop]


    def _getNuclei(self):

        """
        Creates the dataframe containing all the cells of the Spheroid.
        The duplicata clean function is eliminated. Indeed, the local maximum
        makes it impossible for any cell to be segmented twice along the z-axis.
        """

        self._getMaximaFrame(self._getMaskImage(), self.ZRatio, self.RNoyau)


    def _makeSpheroid(self):

        """Generates the spheroid dict containing all the essential information about the spheroid.

        ====== PARAMETERS ======

        df: DataFrame containing all the positional information of the cells
        dCells: maximum distance for two cells to be considered as neighbours
        zRatio: pixel ratio between z and xy dimensions
        Image: original, multichannel, 3D image
        state: True/False variable stating if we seek the dead-alive information

        ====== RETURNS ======

        _Spheroid: dict object"""

        _Spheroid = {}

        _Spheroid['spheroid position'] = self.Position
        _Spheroid['time'] = self.Time
        _Spheroid['cells'] = self._generateCells()

        self.Spheroid = _Spheroid


    def _initializeDead(self):

        X = np.arange(0, 40)
        Y = np.arange(0, 40)
        Z = np.arange(0, 40)
        X, Y, Z = np.meshgrid(X, Y, Z)

        mask = np.sqrt((X-20)**2 + (Y-20)**2 + (Z-20)**2/self.ZRatio**2) < self.RNoyau
        mask = np.transpose(mask, (2,1,0)).astype(np.int)

        deadConv = scipy.signal.fftconvolve(self.DeadImage, mask, mode='same')


        for cellLabel in self.Spheroid['cells'].keys():

            x = int(self.Spheroid['cells'][cellLabel]['x'])
            y = int(self.Spheroid['cells'][cellLabel]['y'])
            z = int(self.Spheroid['cells'][cellLabel]['z'])

            zlen, _, _ = np.nonzero(mask)

            # Test dimension order to verify coherence
            if deadConv[z,x,y]/len(zlen) > self.Thresh:

                print('Dead cell')

                self.Spheroid['cells'][cellLabel]['state'] = 'Dead'


    def _getMaskImage(self):

        blurred = ndimage.gaussian_filter(self.NucImage, sigma=1)
        mask = (blurred > np.percentile(blurred, self.Percentile)).astype(np.float)
        mask += 0.1

        binary_img = mask > 0.5
        binary_img = skimage.morphology.binary_closing(ndimage.binary_dilation(
            ndimage.binary_erosion(binary_img)).astype(np.int))

        return np.multiply(blurred, binary_img)

    def _getMaximaFrame(self, Image, zRatio, rNoyau):

        """Generates the spheroid dict containing all the essential information about the spheroid.

        ====== PARAMETERS ======

        Image: original, multichannel, 3D image
        zRatio: pixel ratio between z and xy dimensions
        rNoyau: radius of the prior (i.e. expected radius of the nuclei)

        ====== RETURNS ======

        _Spheroid: dict object"""

        z, x, y = np.nonzero(Image)
        binary_img = Image > 0 # Regarder si le seuil est le bon ou pas

        mask_image_crop = Image[min(z):max(z), min(x):max(x), min(y):max(y)]

        # Making of the convolution mask

        X = np.arange(0, 40)
        Y = np.arange(0, 40)
        Z = np.arange(0, 40)
        X, Y, Z = np.meshgrid(X, Y, Z)

        mask = np.sqrt((X-20)**2 + (Y-20)**2 + (Z-20)**2/zRatio**2) < rNoyau
        mask = np.transpose(mask, (2,1,0))
        zlen, _, _ = np.nonzero(mask)

        start = timeit.default_timer()

        # We divide by len(zlen) so has to have the average value over the mask

        conv = scipy.signal.fftconvolve(ndimage.gaussian_filter(self.NucImage, sigma=1),
            mask/len(zlen), mode='same')

        mask_conv = ndimage.gaussian_filter(np.multiply(conv, binary_img), sigma=1)

        stop = timeit.default_timer()

        print('Convolution Time: ', stop - start)

        # Distance minimum entre deux noyaux repérés par l'algorithme de convolution
        # On prend exprès 10% de marge supplémentaire pour éviter les contacts inopportuns.
        # On fait aussi attention au liveau minimal d'intensite des pics.

        r = int(self.RNoyau)

        start = timeit.default_timer()

        coordinates = np.asarray(skimage.feature.peak_local_max(mask_conv, min_distance = r,
            threshold_abs = self.ThreshCell))

        stop = timeit.default_timer()

        print('Corrdinates ID Time: ', stop - start)
        print(str(len(coordinates)) + ' cells ID')

        a, b = np.shape(coordinates)

        coordinates = np.hstack((coordinates, np.asarray([np.arange(a)]).T))

        df = pandas.DataFrame(coordinates, columns = ['z', 'x', 'y', 'label'])

        self.NucFrame = df

    def _generateCells(self):

        """ This function serves to generate the cells and gives back a dic object.

        ====== PARAMETERS ======

        df: DataFrame containing the relevant positional information
        dCells: minimum distance between any two cells
        zRatio: ratio between the xy and z dimensions

        ====== RETURNS ======

        _Cells: dict object"""

        _Cells = {}

        df = self.NucFrame
        dCells = self.DCells
        zRatio = self.ZRatio

        df['label'] = df['label'].astype(int).astype(str)

        for label in df['label'].unique():

            dic = {}

            # All values are strings since json doesn't know ints

            dic['x'] = str(df.loc[df['label'] == label, 'x'].iloc[0])
            dic['y'] = str(df.loc[df['label'] == label, 'y'].iloc[0])
            dic['z'] = str(df.loc[df['label'] == label, 'z'].iloc[0])
            dic['neighbours'] = self._nearestNeighbour(df, label, dCells, zRatio)
            dic['state'] = 'Live'

            _Cells[str(int(label))] = dic

        return _Cells


    def _nearestNeighbour(self, df, label, dCells, zRatio):

        """Returns a list of float labels of the cells closest to the given label.
        This method is dependant only on a minimum distance given by the investigator."""

        x = df.loc[df['label'] == label, 'x'].iloc[0]
        y = df.loc[df['label'] == label, 'y'].iloc[0]
        z = df.loc[df['label'] == label, 'z'].iloc[0]

        lf = df.loc[df['label'] != label].copy()

        return lf.loc[np.sqrt((lf['x'] - x)**2 + (lf['y'] - y)**2 +
            (lf['z'] - z)**2/zRatio**2) < dCells, 'label'].values.tolist()

    def _makeG(self):

        G=nx.Graph()
        _Cells = self.Spheroid['cells']
        G.add_nodes_from(_Cells.keys())

        for key in _Cells.keys():

            neighbours = _Cells[key]['neighbours']

            for node in neighbours:

                G.add_edge(key, node)

        return G

    def _refineSph(self):

        G = self._makeG()

        A = nx.betweenness_centrality(G) # betweeness centrality
        B = nx.clustering(G)
        C = nx.degree(G)

        for v in G:

            self.Spheroid['cells'][v]['degree'] = C[v]
            self.Spheroid['cells'][v]['clustering'] = B[v]
            self.Spheroid['cells'][v]['centrality'] = A[v]

        self.Spheroid['N'] = len(self.Spheroid['cells'])
        self.Spheroid['assortativity'] = nx.degree_assortativity_coefficient(G)
        self.Spheroid['average degree'] = np.asarray([float(C[v]) for v in G]).mean()


    def _verifySegmentation(self):

        if not len(self.NucFrame):
            return print('Image doesnt exist')

        if not os.path.exists(self.Path + r'/filmstack/'):
            os.mkdir(self.Path + r'/filmstack/')

        zshape, _, _ = np.shape(self.NucImage)

        ImageAll = self.NucImage

        for n in range(zshape):

            Image = ImageAll[n,:,:]

            plt.figure(figsize=(12, 12))
            plt.subplot(111)
            plt.imshow(Image, vmin=0, vmax=1300, cmap=plt.cm.gray)
            plt.axis('off')

            scalebar = ScaleBar(0.0000003, location = 'lower right')
            # 1 pixel = 0.3 um
            plt.gca().add_artist(scalebar)

            r = self.RNoyau

            for cellLabel in self.Spheroid['cells'].keys():

                x = int(self.Spheroid['cells'][cellLabel]['x'])
                y = int(self.Spheroid['cells'][cellLabel]['y'])
                z = int(self.Spheroid['cells'][cellLabel]['z'])

                if (r**2 - (z -n)**2/self.ZRatio**2) > 0:

                    rloc = np.sqrt(r**2 - (z -n)**2/self.ZRatio**2)
                    s = np.linspace(0, 2*np.pi, 100)
                    x = rloc*np.sin(s) + x
                    y = rloc*np.cos(s) + y

                    if self.Spheroid['cells'][cellLabel]['state'] == 'Dead':

                        plt.plot(y, x, 'r-')

                    else: plt.plot(y, x, 'b-')

            plt.savefig(self.Path + r'/filmstack/im_' + str(n) +'.png')
            plt.close()
