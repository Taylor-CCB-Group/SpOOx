import numpy as np
import pandas as pd
import os

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
        self.pathToWriteOutput = self.baseFolderPath + '/spatial_stats_outputs/'
        self.df = pd.read_csv(self.pathToCellData,delimiter='\t')
        self.points = np.asarray((self.df['x'], self.df['y'])).transpose()
        self.domainX = np.round(np.max(self.df['x']), -3)
        self.domainY = np.round(np.max(self.df['y']), -3)
        

def plotPCF(ax, radii, g, label=None, includeZero=True):
    # Helper function to plot a PCF
    if not includeZero:
        radii = radii[1:]
        g = g[1:]
    ax.plot(radii, g, label=label)
    ax.set_xlabel('r ($\mu$m)')
    ax.set_ylabel('g(r)')
    ax.axhline(y=1,c=[0.5,0.5,0.5],linestyle=':')

def plotPCFWithBootstrappedConfidenceInterval(ax, radii, g, contributions, points, xmax, ymax, label=None, includeZero=True):
    # Helper function to plot a PCF with 95% CI
    # Bootstrapping - work out a confidence interval
    # Following Loh (2008 - https://web.njit.edu/~loh/Papers/ApJ.Loh.2008a.pdf)
    contributions = np.asarray(contributions)
    numContribs = np.shape(contributions)[0]

    # Split domain into equal rectangles
    nRectanglesX = 20
    nRectanglesY = 20
    xRect = np.linspace(0,xmax,nRectanglesX+1)
    yRect = np.linspace(0,ymax,nRectanglesY+1)

    # Identify the rectangle that each point belongs to
    rectangleIDs = np.zeros(shape=(len(contributions)))
    rectID = 0
    for i in range(nRectanglesX):
        for j in range(nRectanglesX):
            accept = (points[:,0] > xRect[i]) & (points[:,0] <= xRect[i+1])
            accept = accept & (points[:,1] > yRect[j]) & (points[:,1] <= yRect[j+1])
            rectangleIDs[accept] = rectID
            rectID = rectID + 1

    # Each bootstrap sample, we select nRectanglesX*nRectanglesY rectangles and construct a PCF from them
    numBootstrapSims = 999
    samplePCFs = np.zeros(shape=(numBootstrapSims, np.shape(contributions)[1]))
    for i in range(numBootstrapSims):
        N = 0
        toSample = np.random.choice(nRectanglesX*nRectanglesY, nRectanglesX*nRectanglesY)
        for j in toSample:
            sample = contributions[rectangleIDs == j]
            samplePCFs[i, :] = samplePCFs[i, :] + sum(sample)
            N = N + np.shape(sample)[0]
        samplePCFs[i, :] = samplePCFs[i, :] / N

    # Get 95% CI
    PCF_min = 2*np.mean(samplePCFs,axis=0) - np.percentile(samplePCFs, 97.5, axis=0)
    PCF_max = 2*np.mean(samplePCFs,axis=0) - np.percentile(samplePCFs, 2.5, axis=0)
    PCF_min = np.insert(PCF_min, 0, 0)
    PCF_max = np.insert(PCF_max, 0, 0)

    if not includeZero:
        radii = radii[1:]
        g = g[1:]
        PCF_min = PCF_min[1:]
        PCF_max = PCF_max[1:]

    ax.plot(radii, g, label=label)
    ax.fill_between(radii, PCF_min, PCF_max, alpha=0.4)
    ax.set_xlabel('r ($\mu$m)')
    ax.set_ylabel('g(r)')
    ax.axhline(y=1, c=[0.5, 0.5, 0.5], linestyle=':')
    return 0


def CalculateBootstrapAroundCSRForPValues(N_A, N_B, domainX, domainY, dr, maxR, numBootstrapSamples):
    from scipy.spatial.distance import cdist
    from multiprocessing import Pool

    bootstrapSamples = np.zeros(shape=(numBootstrapSamples, len(np.arange(0, maxR, dr))))
    for boot in range(numBootstrapSamples):
        print(str(100*(boot+1)/numBootstrapSamples) + '%')
        pointsA = np.random.rand(N_A, 2) * np.asarray([domainX, domainY])
        pointsB = np.random.rand(N_B, 2) * np.asarray([domainX, domainY])

        # Shape (N_A, N_B)
        distances_AtoB = cdist(pointsA, pointsB, metric='euclidean')
        PCF_radii = np.arange(0, maxR, dr)

        pointDensity_A = N_A / (domainX * domainY)
        pointDensity_B = N_B / (domainX * domainY)

        PCF_contributions = []

        # Run at most nParallelPoints points simultaneously
        nParallelPoints = 100
        pointIndexRanges = np.arange(0, N_A, nParallelPoints)
        pointIndexRanges = np.append(pointIndexRanges, N_A)

        for ind_A in range(len(pointIndexRanges) - 1):
            pool = Pool(maxtasksperchild=200)
            minInd = pointIndexRanges[ind_A]
            maxInd = pointIndexRanges[ind_A + 1]
            results = [pool.apply_async(returnPCFContribution, args=(
                point_ind, pointDensity_B, PCF_radii, distances_AtoB[point_ind], pointsA[point_ind, 0], pointsA[point_ind, 1],
                domainX, domainY)) for point_ind in range(minInd, maxInd)]
            output = [p.get() for p in results]

            PCF_contributions.extend(output)
            PCF_values = sum(PCF_contributions)
            # print(str(round(100 * (ind_A + 1) / (len(pointIndexRanges)-1), 2)) + '%')
            pool.close()
            pool.join()

        PCF_values = PCF_values / N_A
        PCF_values = np.insert(PCF_values, 0, 1)
        bootstrapSamples[boot, :] = PCF_values
    return bootstrapSamples

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


def returnAreaOfCircleInDomain(x0, y0, r, domainX, domainY, index):
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


def returnPCFContribution(point_ind, pointDensity, PCF_radii, d, x0, y0, domainX, domainY):
    totalAreas = np.zeros(len(PCF_radii) - 1)
    totalPointsObserved = np.zeros(len(PCF_radii) - 1)

    # areas = [pool.apply(returnAreaOfCircleInDomain, args=(x0, y0, r, domainX, domainY, r_ind)) for r_ind, r in enumerate(PCF_radii)]

    areas = np.zeros(len(PCF_radii))
    for r_ind, r in enumerate(PCF_radii):
        area = returnAreaOfCircleInDomain(x0, y0, r, domainX, domainY, r_ind)
        areas[r_ind] = area

    for r_ind in range(len(PCF_radii) - 1):
        include_points = (d > PCF_radii[r_ind]) & (d <= PCF_radii[r_ind + 1])
        n = sum(include_points)
        area = areas[r_ind + 1] - areas[r_ind]
        totalAreas[r_ind] = totalAreas[r_ind] + area
        totalPointsObserved[r_ind] = totalPointsObserved[r_ind] + n
        # assert(area > 0)
        # PCF_temp[r_ind] = 1.0*n / (area*pointDensity)
    # PCF_values = PCF_values + PCF_temp
    # PCF_values = [v/N for v in PCF_values]
    assert (totalAreas * pointDensity).all() > 0, "Error for point " + str(point_ind) + " at coords (" + str(x0) + "," + str(y0) + ")"
    PCF_contribution = totalPointsObserved / (totalAreas * pointDensity)
    return PCF_contribution


def EfficientPCF(points, domainX, domainY, dr, maxR):
    from scipy.spatial.distance import cdist
    from multiprocessing import Pool

    # D = 500
    # points = ds.points[(ds.points[:, 0] < D) & (ds.points[:, 1] < D)]

    N = len(points)
    distances = cdist(points, points, metric='euclidean')
    PCF_radii = np.arange(0, maxR, dr)
    # PCF_values = np.zeros(len(PCF_radii)-1)
    # corners = [[domainX,0],[0,0],[0,domainY],[domainX,domainY]]

    pointDensity = N / (domainX * domainY)
    # print(pointDensity)

    # tic = time.perf_counter()
    # PCF_contributions = [pool.apply(returnPCFContribution, args=(point_ind, pointDensity, PCF_radii, distances[point_ind], points[point_ind, 0], points[point_ind, 1], domainX, domainY)) for point_ind in range(N)]
    #
    # toc = time.perf_counter()
    # print(toc - tic)

    PCF_contributions = []

    # count = 0

    def log_result(result):
        # This is called whenever returnPCFContribution returns a result.
        # result_list is modified only by the main process, not the pool workers.
        PCF_contributions.append(result)
        # count = count + 1
        # print(count)

    # Run at most nParallelPoints points simultaneously
    nParallelPoints = 100
    pointIndexRanges = np.arange(0, N, nParallelPoints)
    pointIndexRanges = np.append(pointIndexRanges, N)

    for ind in range(len(pointIndexRanges) - 1):
        pool = Pool(maxtasksperchild=200)
        minInd = pointIndexRanges[ind]
        maxInd = pointIndexRanges[ind + 1]
        results = [pool.apply_async(returnPCFContribution, args=(
        point_ind, pointDensity, PCF_radii, distances[point_ind], points[point_ind, 0], points[point_ind, 1], domainX,
        domainY)) for point_ind in range(minInd, maxInd)]
        output = [p.get() for p in results]
        # for point_ind in range(minInd, maxInd):
        #     pool.apply_async(returnPCFContribution, args=(
        #         point_ind, pointDensity, PCF_radii, distances[point_ind], points[point_ind, 0], points[point_ind, 1],
        #         domainX,
        #         domainY), callback=log_result)
        PCF_contributions.extend(output)
        PCF_values = sum(PCF_contributions)
        # print(str(round(100 * (ind + 1) / len(pointIndexRanges), 2)) + '%')
        pool.close()
        pool.join()

    # for point_ind in range(1, nParallelPoints):
    #     pool.apply_async(returnPCFContribution, args=(
    #     point_ind, pointDensity, PCF_radii, distances[point_ind], points[point_ind, 0], points[point_ind, 1], domainX,
    #     domainY), callback=log_result)

    # for point_ind in range(N):
    #     if point_ind % 10 == 0:
    #         print(point_ind, N)
    #
    #     d = distances[point_ind]
    #     x0 = points[point_ind, 0]
    #     y0 = points[point_ind, 1]
    #
    #     PCF_contribution = returnPCFContribution(point_ind, pointDensity, PCF_radii, d, x0, y0, domainX, domainY)
    #     PCF_values = PCF_values + PCF_contribution

    #    PCF_values = sum(PCF_contributions)
    PCF_values = PCF_values / N
    PCF_values = np.insert(PCF_values, 0, 0)
    return PCF_values, PCF_radii


def EfficientPCF_AtoB(pointsA, pointsB, domainX, domainY, dr, maxR):
    from scipy.spatial.distance import cdist
    from multiprocessing import Pool

    N_A = len(pointsA)
    N_B = len(pointsB)
    # Shape (N_A, N_B)
    distances_AtoB = cdist(pointsA, pointsB, metric='euclidean')
    PCF_radii = np.arange(0, maxR, dr)

    pointDensity_A = N_A / (domainX * domainY)
    pointDensity_B = N_B / (domainX * domainY)

    PCF_contributions = []

    # Run at most nParallelPoints points simultaneously
    nParallelPoints = 100
    pointIndexRanges = np.arange(0, N_A, nParallelPoints)
    pointIndexRanges = np.append(pointIndexRanges, N_A)

    for ind_A in range(len(pointIndexRanges) - 1):
        pool = Pool(maxtasksperchild=200)
        minInd = pointIndexRanges[ind_A]
        maxInd = pointIndexRanges[ind_A + 1]
        results = [pool.apply_async(returnPCFContribution, args=(
        point_ind, pointDensity_B, PCF_radii, distances_AtoB[point_ind], pointsA[point_ind, 0], pointsA[point_ind, 1],
        domainX, domainY)) for point_ind in range(minInd, maxInd)]
        output = [p.get() for p in results]

        PCF_contributions.extend(output)
        PCF_values = sum(PCF_contributions)
        # print(str(round(100 * (ind_A + 1) / (len(pointIndexRanges)-1), 2)) + '%')
        pool.close()
        pool.join()

    PCF_values = PCF_values / N_A
    PCF_values = np.insert(PCF_values, 0, 0)
    return PCF_values, PCF_radii, PCF_contributions
