import os
import numpy as np
import pandas
import glob
from tqdm import tqdm_notebook as tqdm
import json

# Maths
from scipy import ndimage, signal

# Image
from skimage import exposure
from skimage import io
from skimage import feature

# Plot
import matplotlib.pyplot as plt
from matplotlib_scalebar.scalebar import ScaleBar

# Graphs
import networkx as nx



class spheroid:

    """ Spheroid class containing the necessary info to build the spheroid.

    ====== NOTICE ======

     - All variables starting with a capital letter refer to class variables.
     - All variables starting with '_' refer to a dict


    ====== PARAMETERS ======

    path: string object, path to the folder storing images
    position: string object, well ID
    time: string object, time of experiment
    zRatio: ratio between x,y and z pixels
    rNoyau: assumed radius of the cell nuclei
    dCells: distance between two nuclei in order for them to be considered in
        contact"""

    def __init__(self, path, position, time, zRatio, rNoyau, dCells):

        self.Path = path
        self.Position = position
        self.Time = time
        self.ZRatio = zRatio
        self.RNoyau = rNoyau
        self.DCells = dCells
        self.NucImage = []
        self.DeadImage = []
        self.Thresh = 500
        self.ThreshCell = 200

    def _loadSpheroid(self, spheroidDict):

        self.spheroid = spheroidDict


    def _loadImageNuclei(self, channelNuc):

        """Function creating the 3D matrix image of the nuclei"""

        image_list = []
        for filename in sorted(glob.glob(path + r'/' + '*.tif')): #assuming tif

            im = io.imread(self.Path + '/' + filename)
            image_list.append(im[channelNuc])

        self.NucImage = np.reshape(image_list, (len(image_list),
            np.shape(image_list[0])[0], np.shape(image_list[0])[1]))

    def _loadImageDead(self, channelDead):

        """Function creating the 3D matrix image of the dead cells"""

        image_list = []
        for filename in sorted(glob.glob(path + r'/' + '*.tif')): #assuming tif

            im = io.imread(self.Path + '/' + filename)
            image_list.append(im[channelDead])

        self.DeadImage = np.reshape(image_list, (len(image_list),
            np.shape(image_list[0])[0], np.shape(image_list[0])[1]))


    def _getNuclei(self):

        if not len(self.NucImage):

            ### Check that image to study does exist

            print('Image doesnt exist')

        self._getMaximaFrame(self._getMaskImage(), self.ZRatio, self.RNoyau)
        self._duplicataClean()


    def _makeSpheroid(self):

        ### COPYPASTE ABOVE

        if not len(self.NucFrame):

            ### Check that image to study does exist

            print('Image doesnt exist')

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

        """Loops through the spheroid dict to modify the state of detected
         dead cells.

         1) We define a mask
         2)convolve through the self.DeadImage matrix
         3) consider any given cell dead if the value of the convolution
         is above a certain threshold."""

        # Creating the mask

        X = np.arange(0, 40)
        Y = np.arange(0, 40)
        Z = np.arange(0, 40)
        X, Y, Z = np.meshgrid(X, Y, Z)

        mask = np.sqrt((X-20)**2 + (Y-20)**2 + (Z-20)**2/self.ZRatio**2) < self.RNoyau
        mask = np.transpose(mask, (2,1,0)).astype(np.int)

        # Convolution

        deadConv = scipy.signal.fftconvolve(self.DeadImage, mask, mode='same')

        # Testing the cells

        for cellLabel in self.Spheroid['cells'].keys():

            x = self.Spheroid['cells'][cellLabel]['x']
            y = self.Spheroid['cells'][cellLabel]['y']
            z = self.Spheroid['cells'][cellLabel]['z']

            zlen, _, _ = np.nonzero(mask)

            if deadConv[z,x,y]/len(zlen) > self.Thresh:
                self.Spheroid['cells'][cellLabel]['state'] = 'Dead'

    def _verifySegmentation(self):

        """Function producing images of the spheroid for analysis."""

        if not len(self.NucFrame):
            return print('Image doesnt exist')

        if not os.path.exists(self.Path + r'/filmstack/'):
            os.mkdir(self.Path + r'/filmstack/')

        zshape, _, _ = np.shape(self.NucImage)

        ImageAll = self.NucImage

        for n in range(zshape):

            Image = ImageAll[n,:,:]

            plt.figure(figsize=(6, 6))
            plt.subplot(111)
            plt.imshow(Image, vmin=0, vmax=800, cmap=plt.cm.gray)
            plt.axis('off')

            scalebar = ScaleBar(0.0000006, location = 'lower right') # 1 pixel = 0.6 umeter
            plt.gca().add_artist(scalebar)

            r = self.RNoyau

            for item, row in self.NucFrame.iterrows():

                if (r**2 - (row['z'] -n)**2/self.ZRatio**2) > 0:

                    rloc = np.sqrt(r**2 - (row['z'] -n)**2/self.ZRatio**2)
                    s = np.linspace(0, 2*np.pi, 100)
                    x = rloc*np.sin(s) + row['x']
                    y = rloc*np.cos(s) + row['y']

                    plt.plot(y, x, 'r-')

            plt.savefig(self.Path + r'/filmstack/im_' + str(n) +'.png')
            plt.close()

    #### UTILITY FUNCTIONS ####

    # Question: how store utility functions in Python class?

    def _getMaskImage(self):

        blurred = ndimage.gaussian_filter(self.NucImage, sigma=2)
        mask = (blurred > np.percentile(blurred, 93)).astype(np.float)
        mask += 0.1

        binary_img = mask > 0.5
        binary_img = skimage.morphology.binary_closing(ndimage.binary_dilation(ndimage.binary_erosion(binary_img)).astype(np.int))

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

        a = mask_image_crop
        zDim, xDim, yDim = np.shape(mask_image_crop)

        X = np.arange(0, 20)
        Y = np.arange(0, 20)
        Z = np.arange(0, 20)
        X, Y, Z = np.meshgrid(X, Y, Z)

        mask = np.sqrt((X-10)**2 + (Y-10)**2 + (Z-10)**2/zRatio**2) < rNoyau
        mask = np.transpose(mask, (2,1,0))
        zlen, _, _ = np.nonzero(mask)

        conv = scipy.signal.fftconvolve(a, mask, mode='same')

        mask_conv = ndimage.gaussian_filter(np.multiply(conv, binary_img[min(z):max(z), min(x):max(x), min(y):max(y)]),
                                            sigma=2)/len(zlen)

        # Distance minimum entre deux noyaux repérés par l'algorithme de convolution
        # On prend exprès 10% de marge supplémentaire pour éviter les contacts inopportuns.
        # On fait aussi attention au liveau minimal d'intensite des pics.

        # We use an absolute threshold for minima detection based upon manual
        # verification of the images.

        coordinates = np.asarray(skimage.feature.peak_local_max(mask_conv,
            threshold_abs=self.ThreshCell, min_distance= self.RNoyau))

        coordinates[:, 0] += min(z)
        coordinates[:, 1] += min(x)
        coordinates[:, 2] += min(y)

        df = pandas.DataFrame(coordinates, columns = ['z', 'x', 'y'])

        for ind, row in df.iterrows():
            df.loc[ind, 'val'] = mask_conv[int(row['z']) - min(z), int(row['x']) - min(x), int(row['y']) - min(y)]
            df.loc[ind, 'label'] = int(ind)

        self.NucFrame = df


    def _duplicataClean(self):

        """Very bright cells tend to be visible much beyond their true position.
        Function cleans away these duplicatas by keeping just the brightest point."""

        df = self.NucFrame

        for ind, row in df.iterrows():

            lf = df.loc[(df['x'] - row['x'])**2 + (df['y'] - row['y'])**2 < 4]

            if len(lf) > 1:

                a = len(df)
                df = df.drop(lf.loc[lf['val'] < lf['val'].max()].index)

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

            dic['x'] = df.loc[df['label'] == label, 'x'].iloc[0]
            dic['y'] = df.loc[df['label'] == label, 'y'].iloc[0]
            dic['z'] = df.loc[df['label'] == label, 'z'].iloc[0]
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

        return lf.loc[np.sqrt((lf['x'] - x)**2 + (lf['y'] - y)**2 + (lf['z'] - z)**2/zRatio**2) < dCells, 'label'].values.tolist()
