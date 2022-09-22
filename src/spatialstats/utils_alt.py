
import numpy as np
import pandas as pd
import os
from math import ceil
import random

class dataset:

    def __init__(self, pathToData):
        head_tail = os.path.split(pathToData)
        file_string = head_tail[0].split('/')
        basePath = os.path.split(head_tail[0])

        self.roi = file_string[-2]
        self.sample = file_string[-3]
        self.indication = file_string[-4]
        self.baseFolderPath = basePath[0]
        #self.name = '_'.join([os.path.splitext(head_tail[1])[0],"_filename"])
        self.name = self.indication + "_" + self.sample + "_" + self.roi
        self.pathToCellData = pathToData
        #self.pathToWriteOutput = self.baseFolderPath + '/spatial_stats_outputs/'
        self.df = pd.read_csv(self.pathToCellData,delimiter='\t')
        self.points = np.asarray((self.df['x'], self.df['y'])).transpose()
        # We estimate domain size by taking most extreme values and rounding up to nearest 50 microns
        self.domainX = ceil(np.max(self.df['x'])/50.0)*50
        self.domainY = ceil(np.max(self.df['y'])/50.0)*50


class dataset_filterSampleID:
    # Alternative creation of dataset with the same format as the above, but generated
    # by passing a sample id. Needed for when people decide that we're going to
    # ignore the decisions we spent ages meeting over for the last year and just
    # merge all the data together into a giant heap.

    def __init__(self, pathToData, sample_id):
        df = pd.read_csv(pathToData,delimiter='\t')
        df = df[df['sample_id'] == sample_id]

        self.df = df

       
        self.sample = df.iloc[0].sample_id
       
        #self.baseFolderPath = basePath[0]
        #self.name = '_'.join([os.path.splitext(head_tail[1])[0],"_filename"])
        self.name = self.sample
        self.pathToCellData = pathToData
        #self.pathToWriteOutput = self.baseFolderPath + '/spatial_stats_outputs/'

        self.points = np.asarray((self.df['x'], self.df['y'])).transpose()
        # We estimate domain size by taking most extreme values and rounding up to nearest 50 microns
        self.domainX = ceil(np.max(self.df['x'])/50.0)*50
        self.domainY = ceil(np.max(self.df['y'])/50.0)*50

def getPCFContributionsWithinGrid(contributions, xmax, ymax, points):
    # Helper function to generate 95% CI around a PCF
    # Following Loh (2008 - https://web.njit.edu/~loh/Papers/ApJ.Loh.2008a.pdf)
    contributions = np.asarray(contributions)
    
    # Split domain into equal rectangles - fix at 100 microns square
    #TODO Make this length scale a parameter - with this fixed, method will become unsuitable 
    # for domains smaller than around 300um x 300um
    rectangleWidthX = 100
    rectangleWidthY = 100
    
    xRect = np.arange(0,xmax+1,rectangleWidthX)
    yRect = np.arange(0,ymax+1,rectangleWidthY)
    
    nRectanglesX = np.shape(xRect)[0]-1
    nRectanglesY = np.shape(yRect)[0]-1
    
    
    # Identify the rectangle that each point belongs to
    rectID = 0
    rectNs = np.zeros(nRectanglesX*nRectanglesY)
    rectContributions = np.zeros((nRectanglesX*nRectanglesY,np.shape(contributions)[1]))
    for i in range(nRectanglesX):
        for j in range(nRectanglesY):
            accept = (points[:,0] > xRect[i]) & (points[:,0] <= xRect[i+1])
            accept = accept & (points[:,1] > yRect[j]) & (points[:,1] <= yRect[j+1])
            if sum(accept) > 0:
                rectContributions[rectID,:] = np.sum(contributions[accept,:],axis=0)
                rectNs[rectID] = sum(accept)
            rectID = rectID + 1
    nRectangles = nRectanglesX*nRectanglesY
    return nRectangles, rectContributions, rectNs


def plotPCFWithBootstrappedConfidenceInterval(ax, radii, g, contributions, points, xmax, ymax, label=None, includeZero=True):
    # Helper function to plot a PCF with 95% CI
    # Bootstrapping - work out a confidence interval
    # Following Loh (2008 - https://web.njit.edu/~loh/Papers/ApJ.Loh.2008a.pdf)
    contributions = np.asarray(contributions)
    numContribs = np.shape(contributions)[0]

    # Split domain into equal rectangles - fix at 100 microns square
    #TODO Make this length scale a parameter - with this fixed, method will become unsuitable 
    # for domains smaller than around 300um x 300um
    rectangleWidthX = 100
    rectangleWidthY = 100

    xRect = np.arange(0,xmax+1,rectangleWidthX)
    yRect = np.arange(0,ymax+1,rectangleWidthY)
    
    nRectanglesX = np.shape(xRect)[0]-1
    nRectanglesY = np.shape(yRect)[0]-1

    # Identify the rectangle that each point belongs to
    # rectangleIDs = np.zeros(shape=(len(contributions)))
    rectID = 0
    rectNs = np.zeros(nRectanglesX*nRectanglesY)
    rectContributions = np.zeros((nRectanglesX*nRectanglesY,np.shape(contributions)[1]))
    for i in range(nRectanglesX):
        for j in range(nRectanglesY):
            accept = (points[:,0] > xRect[i]) & (points[:,0] <= xRect[i+1])
            accept = accept & (points[:,1] > yRect[j]) & (points[:,1] <= yRect[j+1])
            if sum(accept) > 0:
                rectContributions[rectID,:] = np.sum(contributions[accept,:],axis=0)
                rectNs[rectID] = sum(accept)
            # rectangleIDs[accept] = rectID
            rectID = rectID + 1


    # Each bootstrap sample, we select nRectanglesX*nRectanglesY rectangles and construct a PCF from them
    numBootstrapSims = 999
    samplePCFs = np.zeros(shape=(numBootstrapSims, np.shape(contributions)[1]))
    toSample = np.random.choice(nRectanglesX*nRectanglesY, size=(nRectanglesX*nRectanglesY,numBootstrapSims))


    # First sum down each of the 400 boxes (leaves 999 by 30)
    # Then sum down the 400 boxes for the Ns (leaving 999x1)
    # Then divide each entry along the remaining line
    sample = np.sum(rectContributions[toSample,:],axis=0)
    Ns = np.sum(rectNs[toSample],axis=0)

    samplePCFs = sample / Ns[:,np.newaxis]



    # Get 95% CI
    PCF_min = 2*np.mean(samplePCFs,axis=0) - np.percentile(samplePCFs, 97.5, axis=0)
    PCF_max = 2*np.mean(samplePCFs,axis=0) - np.percentile(samplePCFs, 2.5, axis=0)
#    PCF_min = np.insert(PCF_min, 0, 0)
#    PCF_max = np.insert(PCF_max, 0, 0)


    if not includeZero:
        radii = radii[1:]
        g = g[1:]
        PCF_min = PCF_min[1:]
        PCF_max = PCF_max[1:]

    ax.plot(radii, g, label=label)
    ax.fill_between(radii, PCF_min, PCF_max, alpha=0.4)
    ax.set_xlabel('r ($\\mu$m)')
    ax.set_ylabel('g(r)')
    ax.axhline(y=1, c=[0.5, 0.5, 0.5], linestyle=':')
    return 0



def CalculateBootstrapAroundCSRForPValues(N_A, N_B, domainX, domainY, dr_mum, maxR_mum, numBootstrapSamples):
    from scipy.spatial.distance import cdist
    #from multiprocessing import Pool
    #import matplotlib.pyplot as plt

    bootstrapSamples = np.zeros(shape=(numBootstrapSamples, len(np.arange(0, maxR_mum, dr_mum))))

    density_B = N_B / (domainX * domainY)
    PCF_radii = np.arange(0, maxR_mum, dr_mum)
    for boot in range(numBootstrapSamples):
        if boot % 10 == 0:
            print(str(100*(boot)/numBootstrapSamples) + '%')
        points_A = np.random.rand(N_A, 2) * np.asarray([domainX, domainY])
        points_B = np.random.rand(N_B, 2) * np.asarray([domainX, domainY])

        # First calculate areas
        areas_A = getAnnulusAreasAroundPoints(points_A, dr_mum, maxR_mum, domainX, domainY)
        areas_B = getAnnulusAreasAroundPoints(points_B, dr_mum, maxR_mum, domainX, domainY)

        # Shape (N_A, N_B)
        distances_AtoB = cdist(points_A, points_B, metric='euclidean')
        radii, g, contributions = crossPCF(distances_AtoB, areas_A, density_B, dr_mum, maxR_mum)
        bootstrapSamples[boot, :] = np.transpose(g)

    # # Plot halos
    # plt.figure()
    # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),5,axis=1),color=[1,0,0])
    # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),95,axis=1),color=[1,0,0])
    # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),1,axis=1),color=[1,0.6,0.6])
    # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),99,axis=1),color=[1,0.6,0.6])
    # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),0.5,axis=1),color=[1,0.8,0.8])
    # plt.plot(PCF_radii,np.percentile(bootstrapSamples.transpose(),99.5,axis=1),color=[1,0.8,0.8])
    # plt.gca().axhline(1,color=[0,0,0],linestyle=':')
    return PCF_radii, bootstrapSamples


def returnIntersectionPoints(x0, y0, r, domainX, domainY):
    # Calculate the points of intersection between a circle of radius r centred at (x0,y0)
    # and the box boundaries x = 0, y = 0, x = domainX and y = domainY
    # This also includes corners which are within the domain

    # Circle: (x - x0)**2 + (y - y0)**2 = r**2

    intersectionPoints = []
    # Move around anti-clockwise from upper left corner
    # Start with top left corner. Is it in the domain?
    # (0,domainY)
    if x0 ** 2 + (domainY - y0) ** 2 < r ** 2:
        intersectionPoints.append([0, domainY])

    # x = 0
    if x0 < r:
        # find intersection coords with LHS
        if 0 < y0 + np.sqrt(r ** 2 - x0 ** 2) < domainY:
            intersectionPoints.append([0, np.sqrt(r ** 2 - x0 ** 2) + y0])
        if 0 < y0 - np.sqrt(r ** 2 - x0 ** 2) < domainY:
            intersectionPoints.append([0, y0 - np.sqrt(r ** 2 - x0 ** 2)])

    # Bottom left corner, (0,0)
    if x0 ** 2 + y0 ** 2 < r ** 2:
        intersectionPoints.append([0, 0])

    # y = 0
    if y0 < r:
        # find intersection points with bottom
        if 0 < x0 - np.sqrt(r ** 2 - y0 ** 2) < domainX:
            intersectionPoints.append([x0 - np.sqrt(r ** 2 - y0 ** 2), 0])
        if 0 < x0 + np.sqrt(r ** 2 - y0 ** 2) < domainX:
            intersectionPoints.append([x0 + np.sqrt(r ** 2 - y0 ** 2), 0])

    # Bottom right corner, (domainX,0)
    if (domainX - x0) ** 2 + y0 ** 2 < r ** 2:
        intersectionPoints.append([domainX, 0])

    # x = domainX
    if domainX - x0 < r:
        # find intersection points with RHS
        if 0 < y0 - np.sqrt(r ** 2 - (domainX - x0) ** 2) < domainY:
            intersectionPoints.append([domainX, y0 - np.sqrt(r ** 2 - (domainX - x0) ** 2)])
        if 0 < y0 + np.sqrt(r ** 2 - (domainX - x0) ** 2) < domainY:
            intersectionPoints.append([domainX, y0 + np.sqrt(r ** 2 - (domainX - x0) ** 2)])

    # Top right corner, (domainX,domainY)
    if (domainX - x0) ** 2 + (domainY - y0) ** 2 < r ** 2:
        intersectionPoints.append([domainX, domainY])

    # y = domainY
    if domainY - y0 < r:
        # find intersection points with top
        if 0 < x0 + np.sqrt(r ** 2 - (domainY - y0) ** 2) < domainX:
            intersectionPoints.append([x0 + np.sqrt(r ** 2 - (domainY - y0) ** 2), domainY])
        if 0 < x0 - np.sqrt(r ** 2 - (domainY - y0) ** 2) < domainX:
            intersectionPoints.append([x0 - np.sqrt(r ** 2 - (domainY - y0) ** 2), domainY])
    return intersectionPoints



def returnAreaOfCircleInDomain(x0, y0, r, domainX, domainY):
    intersectionPoints = returnIntersectionPoints(x0, y0, r, domainX, domainY)

    if not intersectionPoints:
        area = np.pi * r ** 2
    else:
        # Need to calculate area from intersection Points
        intersectionPoints.append(intersectionPoints[0])
        area = 0
        for v in range(len(intersectionPoints) - 1):
            a = intersectionPoints[v]
            b = intersectionPoints[v + 1]
            # Find out if this is a segment or a triangle
            isTriangle = False

            # Check if point b is anticlockwise from point a on the same line
            if a[0] == b[0] or a[1] == b[1]:
                if a[0] == b[0]:  # On a vertical line
                    if a[0] == 0:
                        # LHS
                        if b[1] < a[1]:
                            isTriangle = True
                    else:
                        # RHS
                        if b[1] > a[1]:
                            isTriangle = True
                else:  # On a horizontal line
                    if a[1] == 0:
                        # bottom
                        if b[0] > a[0]:
                            isTriangle = True
                    else:
                        # top
                        if a[0] > b[0]:
                            isTriangle = True

            # If points are on the same line moving anticlockwise, then return the area of the triangle formed by a, b and the centre
            if isTriangle:
                # Points share a border: return area of triangle between them
                area = area + 0.5 * np.abs(a[0] * (b[1] - y0) + b[0] * (y0 - a[1]) + x0 * (a[1] - b[1]))
            else:
                # Else, return the area of the circle segment between them
                # We need to be careful to take the angle between v1 and v2 in an anticlockwise direction
                v1 = [x0 - a[0], y0 - a[1]]
                v2 = [x0 - b[0], y0 - b[1]]

                theta = np.arctan2(v2[1], v2[0]) - np.arctan2(v1[1], v1[0])
                # Normalise to 0, 2pi
                if theta < 0:
                    theta = theta + 2 * np.pi

                area = area + 0.5 * theta * r ** 2
    return area


def returnAreaOfCircleInDomainAroundPoint(index, points, r, domainX, domainY):
    point = points[index,:]
    area = returnAreaOfCircleInDomain(point[0], point[1], r, domainX, domainY)
    return area




def getAnnulusAreasAroundPoints(points_i, dr, maxR, domainX, domainY):
    # We want to populate a table the same size as distances, which contains the area of the annulus containing that contribution
    # i.e., "at distance D(i->j) from point i, what is area of containing annulus?"

    vfunc_returnAreaOfCircleInDomainAroundPoint = np.vectorize(returnAreaOfCircleInDomainAroundPoint,excluded=['points'])

    PCF_radii_lower = np.arange(0, maxR, dr)
    PCF_radii_upper = np.arange(dr, maxR + dr, dr)

    allAreas = np.zeros(shape=(len(points_i),len(PCF_radii_lower)))

    for annulus in range(len(PCF_radii_lower)):
        inner = PCF_radii_lower[annulus]
        outer = PCF_radii_upper[annulus]

        areas_in = vfunc_returnAreaOfCircleInDomainAroundPoint(index=np.arange(len(points_i)), points=points_i, r=inner, domainX=domainX, domainY=domainY)
        areas_out = vfunc_returnAreaOfCircleInDomainAroundPoint(index=np.arange(len(points_i)), points=points_i, r=outer, domainX=domainX, domainY=domainY)
        areas = areas_out - areas_in

        allAreas[:,annulus] = areas
    return allAreas



def crossPCF(distances_AtoB, areas_A, density_B, dr_mum, maxR_mum):
    N_A = np.shape(distances_AtoB)[0]

    PCF_radii_lower = np.arange(0, maxR_mum, dr_mum)
    PCF_radii_upper = np.arange(dr_mum, maxR_mum + dr_mum, dr_mum)

    crossPCF_AtoB = np.ones(shape=(len(PCF_radii_lower),1))
    contributions = np.zeros(shape=(N_A,len(PCF_radii_lower)))
    for annulus in range(len(PCF_radii_lower)):
        inner = PCF_radii_lower[annulus]
        outer = PCF_radii_upper[annulus]

        # Find pairwise distances within this radius
        distanceMask = np.logical_and((distances_AtoB > inner),(distances_AtoB <= outer))
        for i in range(N_A):
            # For each point in pA
            # Find pairwise distances to points in pB within this radius
            fillIndices = np.where(distanceMask[i,:])[0]
            contribution = len(fillIndices)/(density_B*areas_A[i,annulus])
            crossPCF_AtoB[annulus] = crossPCF_AtoB[annulus] + contribution
            contributions[i,annulus] = contributions[i,annulus] + contribution
        crossPCF_AtoB[annulus] = crossPCF_AtoB[annulus] / N_A
    return PCF_radii_lower, crossPCF_AtoB, contributions

def changeSomeElements(matrix):
    
    # Select elements of a submatrix (a b; c d) such that elements in the same row/column are from the same row/column in matrix
    n, m = np.shape(matrix)
    rows = random.sample(range(n), 2)
    cols = random.sample(range(m), 2)
    a = matrix[rows[0],cols[0]]
    b = matrix[rows[0],cols[1]]
    c = matrix[rows[1],cols[0]]
    d = matrix[rows[1],cols[1]]
    
    # Now find the smallest values on the diagonals
    minDiag1 = min(a,d)
    minDiag2 = min(b,c)
    if minDiag1 == 0 and minDiag2 == 0:
        return matrix,False
    else:
        # At least one diagonal doesn't include 0. We subtract from that diagonal. wlog pick diag1 if both are fine
        if minDiag1 > 0:
            if minDiag1 == 1:
                k = 1
            else:
                # Choose k between 1 and minDiag1
                k = np.random.randint(1,minDiag1+1)
            new_a = a - k
            new_d = d - k
            new_b = b + k
            new_c = c + k
        else:
            if minDiag2 == 1:
                k = 1
            else:
                # Choose k between 1 and minDiag2
                k = np.random.randint(1,minDiag2+1)
            new_a = a + k
            new_d = d + k
            new_b = b - k
            new_c = c - k
        #print(rows, cols)
        matrix[rows[0],cols[0]] = new_a
        matrix[rows[0],cols[1]] = new_b
        matrix[rows[1],cols[0]] = new_c
        matrix[rows[1],cols[1]] = new_d
        return matrix,True










