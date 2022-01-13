#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 16 10:40:55 2020

@author: quenot
"""

from xml.dom import minidom
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from scipy.ndimage import convolve
from matplotlib import pyplot as plt
import time
from numba import jit
import xlrd
import pandas as pd
from scipy.signal import fftconvolve
from skimage.transform import resize as resize2

class Detector:
    def __init__(self, exp_dict):
        self.xmlDetectorFileName="xmlFiles/Detectors.xml"
        self.xmldocDetector = minidom.parse(self.xmlDetectorFileName)
        self.myName=""
        self.myDimensions=(0,0)
        self.myPixelSize=0. #en um
        self.myPSF=0. #en pixel
        self.myEfficiencyLimit=100. #en kev
        self.margins=exp_dict['margin']
        self.myEnergyLimit=200 #No longer useful?
        self.myBinsThersholds=[] #keX
        self.myScintillatorMaterial=None
        self.myScintillatorThickness=0. #um
        self.beta=[]
        
        
    def defineCorrectValuesDetector(self):
        """
        Define detector parameters from the xml file

        Raises:
            ValueError: detector not found in xml file.

        Returns:
            None.

        """
        for currentDetector in self.xmldocDetector.documentElement.getElementsByTagName("detector"):
            correctDetector = self.getText(currentDetector.getElementsByTagName("name")[0])
            if correctDetector == self.myName:
                self.myDimensions=self.getMyDimensions(currentDetector)
                self.myPixelSize=float(self.getText(currentDetector.getElementsByTagName("myPixelSize")[0]))
                self.myPSF=float(self.getText(currentDetector.getElementsByTagName("myPSF")[0]))
                
                for node in currentDetector.childNodes:
                    if node.localName=="myEnergyLimit":
                        self.myEnergyLimit=float(self.getText(currentDetector.getElementsByTagName("myEnergyLimit")[0]))
                    if node.localName=="myBinsThersholds":
                        myBinsThersholdsTmp=self.getText(currentDetector.getElementsByTagName("myBinsThersholds")[0])
                        myBinsThersholdsTmp=list(myBinsThersholdsTmp.split(","))
                        self.myBinsThersholds=[float(ele) for ele in myBinsThersholdsTmp]
                    if node.localName=="myScintillatorMaterial":
                        self.myScintillatorMaterial=(self.getText(currentDetector.getElementsByTagName("myScintillatorMaterial")[0]))
                        self.myScintillatorThickness=float(self.getText(currentDetector.getElementsByTagName("myScintillatorThickness")[0]))
                return
            
        raise ValueError("detector not found in xml file")
            
            
    def detection(self,incidentWave,effectiveSourceSize,samplingFactor):
        
        """
        Adds source and PSF blurrings, resamples to detector pixel size and add shot noise

        Args:
            incidentWave (2d numpy array): Intensity arriving to the detector.
            effectiveSourceSize (float): source projected FWHM.

        Returns:
            detectedImage (2d numpy array): detected image.

        """
        sampledpixeledsize=self.myPixelSize/samplingFactor
        # If NEAM effectiveSourceSize != effectiveSourceSize but factor to multiply source size (in um) by
        weights = self.source_from_file(20, effectiveSourceSize, sampledpixeledsize)

        
        if effectiveSourceSize:
            # sigmaSource=effectiveSourceSize/2.355 #from FWHM to std dev
            # incidentWave=gaussian_filter(incidentWave, sigmaSource,mode='wrap')
            # Neam's Changes
            padSize=np.max(weights.shape)
            incidentWave=np.pad(incidentWave,(padSize,padSize), mode='edge')
            incidentWave=fftconvolve(incidentWave, weights, mode='same')
            incidentWave=incidentWave[padSize:-padSize,padSize:-padSize]
            #incidentWave=convolve(incidentWave, weights, mode='nearest')
        intensityBeforeDetection=resize(incidentWave, self.myDimensions[0],self.myDimensions[1])
        seed       = int(np.floor(time.time()*100%(2**32-1)))
        rs         = np.random.RandomState(seed)
        if self.myPSF:
            detectedImage=gaussian_filter(intensityBeforeDetection, self.myPSF,mode='wrap')
        else:
            detectedImage=intensityBeforeDetection
        
        detectedImage = rs.poisson((detectedImage))

        detectedImage=detectedImage[self.margins:self.myDimensions[0]-self.margins,self.margins:self.myDimensions[1]-self.margins]
        
        return detectedImage
        
    
    def getText(self,node):
        return node.childNodes[0].nodeValue
    
    def getMyDimensions(self,node):
        dimX=int(self.getText(node.getElementsByTagName("dimX")[0]))+self.margins*2
        dimY=int(self.getText(node.getElementsByTagName("dimY")[0]))+self.margins*2
        return np.array([dimX ,dimY])
    
    
    def getBeta(self, sourceSpectrum):
        print("Materials :", self.myScintillatorMaterial)
        pathTablesDeltaBeta ='Samples/DeltaBeta/TablesDeltaBeta.xls'
        for sh in xlrd.open_workbook(pathTablesDeltaBeta).sheets():
            for col in range(sh.ncols):
                row=0
                myCell = sh.cell(row, col)
#                    print("\n\n Sample materials : ",self.myMaterials[imat])
                if myCell.value == self.myScintillatorMaterial:
                    row=row+3
                    for energy,_ in sourceSpectrum:
                        currentCellValue=sh.cell(row,col).value
                        if energy*1000<currentCellValue:
                            self.beta.append((energy,1))
                            print("No delta beta values under",currentCellValue, "eV")
                            continue
                        nextCellValue=sh.cell(row+1,col).value
                        while nextCellValue<energy*1e3: #find in which interval the energy is in the delta beta tables
                            row+=1
                            currentCellValue=sh.cell(row,col).value
                            nextCellValue=sh.cell(row+1,col).value
                        #Linear interpollation between those values
                        step=nextCellValue-currentCellValue
                        currentCellBeta=float(sh.cell(row,col+2).value)
                        nextCellBeta=float(sh.cell(row+1,col+2).value)
                        betaInterp=abs(nextCellValue-energy*1e3)/step*currentCellBeta+abs(currentCellValue-energy*1e3)/step*nextCellBeta
                        self.beta.append((energy,betaInterp))
                        
                    return
            raise ValueError("The scintillator material has not been found in delta beta tables")
    
<<<<<<< HEAD
    def source_geometry(self, distance, sigma,pixel_size):
    
            pixel_size=pixel_size*1e-6 
            x = 5*sigma/pixel_size
            y = distance/pixel_size+5*sigma/pixel_size
            y0, x0=y/2-distance/pixel_size/2, x/2
            y1, x1=y/2+distance/pixel_size/2, x/2
        
            x=np.arange(x)
            y=np.arange(y)
            sigma /= pixel_size
    
            gx = np.exp(-(x-x0)**2/(2*sigma**2))
            gx1 = np.exp(-(x-x1)**2/(2*sigma**2))

            gy = np.exp(-(y-y0)**2/(2*sigma**2))
            gy1 = np.exp(-(y-y1)**2/(2*sigma**2))

            g = np.outer(gx, gy)
            g1 = np.outer(gx1, gy1)
            g2=(g+g1)/2
            g=g2/np.sum(g2)
            # plt.figure()
            # plt.imshow(g2)
            # plt.colorbar()
            # plt.show()
            return g2

    def open_source_file(self, filepath = "sources/LFXT.txt"):
        headers = ['x', 'y', 'z', 'xs', 'ys', 'E']
        source = pd.read_csv(filepath)
        source.pop(" Charge: -1.5730242e-16 C")
        source = source[1:]
        source.columns = ['DATA']
        
        source[headers] = source.DATA.str.split(" ", expand=True)
        source.pop('DATA')

        data = {}
        for header in headers:
            data[header] = source[header].to_numpy().astype(float)
        return data

    def source_from_file(self, bins, scaling, sampled_pixel_size, filepath = "sources/LFXT.txt"):
        data = self.open_source_file(filepath)
        x = data['x']
        y = data['y']
        range_x=(np.max(x)-np.min(x))/1e-4
        range_y=(np.max(y)-np.min(y))/1e-4
        print(scaling, sampled_pixel_size)
        print(range_x,range_y)
        nbinx=round(range_x*scaling)
        nbiny=round(range_y*scaling)
        print(nbinx,nbiny)

        hist, pixel_x, pixel_y = np.histogram2d(x,y, (nbinx,nbiny))
        remove_columns = [all(k < 0.1*max(hist.flatten()) for k in j) for j in hist]
        hist = np.array([i for i, j in zip(hist, remove_columns) if j is False])
        remove_rows = [all(k < 0.1*max(hist.flatten()) for k in j) for j in hist.T]
        hist = [i for i, j in zip(hist.T, remove_rows) if j is False]
        hist /= np.sum(hist)
        extent = np.array([min(pixel_x), max(pixel_x), min(pixel_y), max(pixel_y)])

        plt.figure('Pre Stretch')
        plt.imshow(hist, extent=extent, aspect='auto')
        plt.colorbar()
        plt.xlabel('Cm')
        plt.show(block=True)
        
        # hist = np.array(hist)
        # x_size = abs(pixel_x[-1]-pixel_x[-2])
        # y_size = abs(pixel_y[-1]-pixel_y[-2])
        # x_size *= scaling
        # y_size *= scaling
        # scale_factor = np.round(sampled_pixel_size*1e-4/np.array(x_size, y_size))
        # hist = resize2(hist, np.array(hist.shape)*scale_factor,order=1)
        # print(scale_factor)
        # plt.figure('Post Stretch')
        # plt.imshow(hist, extent=extent*scaling*scale_factor, aspect='auto')
        # plt.colorbar()
        # plt.show(block=True)
        return hist


    def source_from_file2(self, bins, scaling, sampled_pixel_size, filepath = "sources/LFXT.txt"):
        data = self.open_source_file(filepath)
        x = data['x']
        y = data['y']
        hist, pixel_x, pixel_y = np.histogram2d(x,y, bins)
        remove_columns = [all(k < 0.1*max(hist.flatten()) for k in j) for j in hist]
        hist = np.array([i for i, j in zip(hist, remove_columns) if j is False])
        remove_rows = [all(k < 0.1*max(hist.flatten()) for k in j) for j in hist.T]
        hist = [i for i, j in zip(hist.T, remove_rows) if j is False]
        hist /= np.sum(hist)
        extent = np.array([min(pixel_x), max(pixel_x), min(pixel_y), max(pixel_y)])

        plt.figure('Pre Stretch')
        plt.imshow(hist, extent=extent, aspect='auto')
        plt.colorbar()
        plt.xlabel('Cm')
        plt.show(block=True)
        
        hist = np.array(hist)
        x_size = abs(pixel_x[-1]-pixel_x[-2])
        y_size = abs(pixel_y[-1]-pixel_y[-2])
        x_size *= scaling
        y_size *= scaling
        scale_factor = np.round(sampled_pixel_size*1e-4/np.array(x_size, y_size))
        hist = resize2(hist, np.array(hist.shape)*scale_factor,order=1)
        print(scale_factor)
        plt.figure('Post Stretch')
        plt.imshow(hist, extent=extent*scaling*scale_factor, aspect='auto')
        plt.colorbar()
        plt.show(block=True)
        return hist

      
=======
    def source_geometry(distance, sigma,pixel_size,x=2000,y=2000):
    
        y0, x0=y/2-distance/pixel_size, x/2
        y1, x1=y/2+distance/pixel_size, x/2
            
        x=np.arange(x)
        y=np.arange(y)
        sigma /= pixel_size
        
        gx = np.exp(-(x-x0)**2/(2*sigma**2))
        gx1 = np.exp(-(x-x1)**2/(2*sigma**2))

        gy = np.exp(-(y-y0)**2/(2*sigma**2))
        gy1 = np.exp(-(y-y1)**2/(2*sigma**2))

        g = np.outer(gx, gy)
        g1 = np.outer(gx1, gy1)

        
        return g+g1

>>>>>>> 9b4ea9c (idk)
    
@jit(nopython=True)
def resize(imageToResize,sizeX, sizeY):
    Nx, Ny=imageToResize.shape
    if Nx==sizeX and Ny==sizeY:
        return imageToResize
    
    resizedImage=np.ones((sizeX,sizeY))
    sampFactor=int(Nx/sizeX)
    
    for x0 in range(sizeX):
        for y0 in range(sizeY):
            resizedImage[x0,y0]=np.sum(imageToResize[int(np.floor(x0*sampFactor)):int(np.floor(x0*sampFactor+sampFactor)),int(np.floor(y0*sampFactor)):int(np.floor(y0*sampFactor+sampFactor))])
            
    return resizedImage


    



