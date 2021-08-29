#!/usr/bin/env python3
"""
Author: Junda Huang 
Student number: 910203370050
Implementation of the k-means clustering algorithm

Hints:
- write a function to obtain Euclidean distance between two points.
- write a function to initialize centroids by randomly selecting points 
    from the initial set of points. You can use the random.sample() method
- write a function to find the closest centroid to each point and assign 
    the points to the clusters.     
- write a function to calculate the centroids given the clusters
- write a function to implement k-means
- write a function to calculate WGSS given the clusters and the points
- write a function to calculate BGSS given the clusters and the points
"""

# import statements
from math import sqrt
from random import choices, sample
from copy import deepcopy
from statistics import mean, stdev

# Functions
def csv_parser(lines):
    """Return list of point coordinates as [[x1,y1,z1,...],[x2,y2,z2,...]]
    
    lines: open file or list of lines. Expected format:
        The file has a single header line. 
        Each line contains the coordinates for one data point, starting
        with a label. A data point can be specified in arbitrary dimensions.

    Output: List of lists with the coordinates of the points.
    This function does not capture the labels of the data points. In case
    they are needed later, this function should be adjusted. 
    """ 

    data_points = []
    for line in lines:
        items = line.strip().split(",")
        try: #will fail on header line in file
            data_points.append(list(map(float, items[1:]))) #skip label
        except ValueError: #must be the header
            continue
    return data_points

def average_point(points):
    """
    calculate the average point of given data points

    input:
        points: (list) of list as data point coords of floats
    output:
        average_p: (list) of average point coords of floats
    """
    average_p = [sum(m)/len(points) for m in zip(*points)]
    return average_p

def euclidean(p1, p2):
    """
    Calculate the euclidean distance between two data points

    input: 
        p1, p2: (list or tuple) of data points coordinates in float or int
            e.g. p1 = [5.3, 4.6]; p2 = (3.2, 7.1)
    output:
        dis: (float) distance of two input data points

    """
    if len(p1) != len(p2):
        raise ValueError
    else:
        dist = sqrt(sum([(p1[c]-p2[c]) ** 2 for c in range(len(p1))]))
    return dist

def centroids_init(data_points, k):
    """
    Randomly select k points in data as the centroids

    input:
        data_points: (list) of list as data point coords of floats
        k: (int) number of clusters
    output:
        centr_points: (list) of k lists of randomly selected datapoints
    """
    centr_points = list(choices(data_points, k = k))
    return centr_points

def clustering_init(data_points, centr_points):
    """
    Find the closest centroid to each point and assign to clusters

    input:
        data_points: (list) of list as data point coords of floats
        centr_points: (list) of k lists of randomly selected datapoints
    output:
        clusters: (list) of list as index of dataset
            e.g. [[1, 3, 5], [2, 4, 6], [7, 8, 9]] 
        centroids: (list) of clusters' centroids
    """
    k = len(centr_points)
    clusters = []
    cps = [] #initalise a list for in-between centr_points
    for i in range(k):
        if centr_points[i] in data_points:
            clusters.append([data_points.index(centr_points[i])])
        else:
            clusters.append([])
            # data points as orignal centroids are added in clusters 
        cps.append(centr_points[i]) 
    for j, dp in enumerate(data_points):
        if dp in centr_points: # data points as orignal centroids are ignored
            continue
        else:
            cluster = 0
            closest = 0 
            for l, cp in enumerate(cps):
                dis = euclidean(cp, dp)
                if l == 0:
                    closest = dis
                elif dis < closest:
                    closest = dis
                    cluster = l
                else:
                    continue
        clusters[cluster].append(data_points.index(dp))
        # append to cluster
        points = [data_points[clusters[cluster][n]] \
            for n in range(len(clusters[cluster]))]
        cps[cluster] = average_point(points)
        # calculate new average for edited cluster and changing centr_point
    centroids = deepcopy(cps)
    return clusters, centroids

def kmean_cluster(data_points, centroids):
    """
    iterate kmean cluster until centroids are stablised

    input:
        data_points: (list) of list as data point coords of floats
        centroids: (list) of clusters' centroids initially generated
    output:
        kmean_clusters: (list) of list as index of dataset 
            e.g. [[1, 3, 5], [2, 4, 6], [7, 8, 9]]
        new_centroid: (list) of clusters' final centroids
        iteration_n: (int) iteration number
    """
    new_centroids = [[], centroids, []]
    i = 1
    while new_centroids[-3] != new_centroids[-2]:
        # if the centroids stays the same then stop
        kmean_clusters, new_centroids[i + 1] = \
            clustering_init(data_points, new_centroids[i])
        i += 1
        new_centroids.append([])
    iteration_n = i
    new_centroid = new_centroids[iteration_n]
    return kmean_clusters, new_centroid, iteration_n

def wgss_c(cluster, centroid, data_points):
    """
    calculate WGSS for given cluster 

    input:
        cluster: (list) of int as index of dataset. 
            e.g. [1, 3, 5]
        centroid: (list) of given clusters' centroid
        data_points: (list) of list as data point coords of floats
    output:
        wgss_c: (int) value of WGSS of given cluster
    """
    total_dis = sum([euclidean(data_points[cluster[j]],\
             centroid) ** 2 for j in range(len(cluster))])
    wgss_c = total_dis/len(data_points)
    return wgss_c

def wgss(kmean_clusters, new_centroid, data_points):
    """
    calculate total WGSS for given clusters

    input:
        kmean_clusters: (list) of list as index of dataset
            e.g. [[1, 3, 5], [2, 4, 6], [7, 8, 9]] 
        new_centroid: (list) of clusters' final centroid
        data_points: (list) of list as data point coords of floats
    output:
        wgss: (int) total value of WGSS of given clusters
    """
    wgss = 0
    for i, cluster in enumerate(kmean_clusters):
        wgss += wgss_c(cluster, new_centroid[i], data_points)
    return wgss

def bgss(kmean_clusters, new_centroid, data_points):
    """
    calculate BGSS for given clusters

    input:
        kmean_clusters: (list) of list as index of dataset
            e.g. [[1, 3, 5], [2, 4, 6], [7, 8, 9]] 
        new_centroid: (list) of clusters' final centroid
        data_points: (list) of list as data point coords of floats
    output:
        bgss: (int) total value of BGSS of given clusters
    """
    bgss = 0
    for i, centroid in enumerate(new_centroid):
        bgss += len(kmean_clusters[i]) * \
            (euclidean(centroid, average_point(data_points)) ** 2)
    return bgss

def clusters_quality(clusters_list, centroid_list, data_points):
    """
    Evaluate clusters quality and return W value

    input:
        clusters_list: (list) of clusters set as
            list of lists of of indexes of dataset
            e.g. [[1, 3, 5], [2, 4, 6], [7, 8, 9]], 
                [[3, 4, 5], [1, 2, 6], [7, 8, 9]]  
        centroid_list: (list) of clusters sets' final centroids as lists
        data_points: (list) of list as data point coords of floats
    output:
        clusters_qual: (list) clusters as list of lists and 
            W value as int as tuple
             e.g. [([[1, 3, 5], [2, 4, 6], [7, 8, 9]], 3.4),
                 ([[3, 4, 5], [1, 2, 6], [7, 8, 9]], 5.7)]
    """
    clusters_qual = []
    for clusters, centroids in zip(clusters_list, centroid_list):
        w = wgss(clusters, centroids, data_points)/\
            bgss(clusters, centroids, data_points)
        clusters_qual.append((clusters, w))
    return clusters_qual

def result_output(filename, k, rep, resampling):
    """
    To obtain answer needed for answer questions by using defined functions

    input:
        filename: (string) of the csv filename
        k: (int) number of clusters
        rep: (int) number of reptitions of kmeans required
        resampling: (string) of resampling methods:
            jackknifing, boottrapping or subsampling
    output:
        average_iter: (float) average iterations needed for converge
        stdev_iter: (float) standerd deviation of iterations
        kmeans: (list) of list as index of dataset
            e.g. [[1, 3, 5], [2, 4, 6], [7, 8, 9]]
        average_w: (float) average W value of the given clusters set
        stdev_w: (float) standerd deviation of W values
        best_w: (tuple) of a cluster with the lowest W value and its W value 
    """
    if resampling == None: # check if resampling is activated
        with open(filename) as lines:
            data_points = csv_parser(lines)
    else: 
        with open(filename) as lines:
            data = csv_parser(lines)
        data_points = resample(data, resampling = resampling)
    iterations = [] # store numbers of iteration in a list of int
    kmeans = [] # store different sets of clusters in list if lists 
    centers = [] # store centroids of each set of clusters
    for i in range(rep): # run rep times kmean
        centr_points = centroids_init(data_points, k = k)
        clusters, centroids = clustering_init(data_points, centr_points)
        kmean_clusters, new_centroids, iteration\
             = kmean_cluster(data_points, centroids)
        iterations.append(iteration)
        kmeans.append(kmean_clusters)
        centers.append(new_centroids)
    average_iter = mean(iterations)
    stdev_iter = stdev(iterations)
    clusters_qual = clusters_quality(kmeans, centers, data_points)
    w_list = []
    for i in clusters_qual:
        w_list.append(i[1])
    average_w = mean(w_list)
    stdev_w = stdev(w_list)
    best_w = min(clusters_qual, key = lambda t: t[1])
    return average_iter, stdev_iter, kmeans, centers, \
        w_list, average_w, stdev_w, best_w  

def k_optimization(filename, k_list, rep, resampling):
    """
    Using defined functions to find out the best k suited for clusters

    input:
        filename: (string) of the csv filename
        k_list: (list) of int of numbers of clusters to be tested
        rep: (int) number of repetition of kmeans required
        resampling: (string) of resampling methods:
            jackknifing, boottrapping or subsampling
    output:
        k_dict: (dictionary)
        best_k: (int) best numbers for clustering
        best_clusters: (tuple) best clustered clusters set with its W value. 
    """
    k_dict = {}
    for k in k_list:
        k_dict[k] = []
        k_dict[k].extend((result_output\
            (filename2, k = k, rep = rep, resampling = resampling)))
    kmin_list = []
    test = []
    for key, value in k_dict.items():
        kmin_list.append((key, value[5]))
        test.append(value[-1])
    best_k = min(kmin_list, key = lambda t: t[1])[0]
    best_clusters_k = k_dict[best_k][-1]
    best_clusters_all = min(test, key = lambda t: t[1])
    if resampling == None:
        return k_dict, best_k, best_clusters_k, best_clusters_all
    else: # if resampling activated, this function only return the k value
        return best_k

def boottrapping(data_points):
    """
    Resample data set by boottrapping - sampling with replacement

    input:
        data_points: (list) of list as data point coords of floats
    output:
        boottrap_points: (list) of list as data point coords of floats
    """
    boottrap_points = choices(data_points, k = len(data_points))
    return boottrap_points

def subsampling(data_points):
    """
    Resample data set by subsampling - sampling without replacement

    input:
        data_points: (list) of list as data point coords of floats
    output:
        subsample_points: (list) of list as data point coords of floats
    """
    subsample_points = sample(data_points, k = len(data_points))
    return subsample_points

def jackknifing(data_points):
    """
    Resample data set by jackknifing - randomly remove one obeservation

    input:
        data_points: (list) of list as data point coords of floats
    output:
        jackknife_points: (list) of list as data point coords of floats
    """
    jackknife_points = sample(data_points, k = len(data_points) - 1)
    return jackknife_points

def resample(data, resampling):
    """
    Choose esample data set by jackknifing, boottrapping or subsampling

    input:
        data: (list) of list as data point coords of floats
        resampling: (string) of resampling methods:
            jackknifing, boottrapping or subsampling
    output:
        data_points: (list) of list as data point coords of floats
    """
    data_points = []
    if resampling == 'jackknifing':
        data_points = jackknifing(data)
    elif resampling == 'boottrapping':
        data_points = boottrapping(data)
    elif resampling == 'subsampling':
        data_points = subsampling(data)
    else:
        raise ValueError
    return data_points

if __name__ == "__main__":

    # the code below should produce the results necessary to answer
    # the questions. In other words, if we run your code, we should see 
    # the data that you used to answer the questions.
    
    # Question 1: 
    filename1 = '2dtest.csv'
    average_iter, stdev_iter, kmeans, centers, w_list, average_w, \
        stdev_w, best_w = result_output\
            (filename1, k = 3, rep = 20, resampling = None)
    # 1a:
    print('Answer to question 1:\n1a:\n\
    The average iterations needed for converge is {:.3f}.\n\
    The standerd deviation of iterations is {:.3f}.'.format\
                (average_iter, stdev_iter))
    # 1b:
    print('1b:\nThe clusters are as follow:')
    for index, clusters in enumerate(kmeans):
        print(index + 1, clusters)
    # 1c:
    print('1c:\nThe W values are {}.\n\
        The average W value is {:.3f}.\n\
        The standerd deviation of W value is {:.3f}.'\
            .format(w_list, average_w, stdev_w))
    # 1d:
    print('1d:\nBest clusters with lowest w value is {}'.format(best_w))

    # Question 2:
    filename2 = 'LargeSet_1.csv'
    k_large1 = (2, 3, 4, 5, 6)
    k_large1_dict, best_k, best_clusters_k, best_clusters_all = \
        k_optimization(filename2, k_large1, rep = 10, resampling = None)
    # 2a:
    print('Answer to question 2:\n2a:')
    for key, value in k_large1_dict.items():
        print('The average iterations needed for converge k = {} is {:.3f}.\n\
    The standerd deviation of iterations of k = {} is {:.3f}.'.format\
                (key, value[0], key, value[1]))
    # 2b:
    print('2b:\n')
    for key, value in k_large1_dict.items():
        print('The average W values of k = {} is {:.3f}.\n\
        The standerd deviation of W value of k = {} is {:.3f}.'\
            .format(key, value[-3], key, value[-2]))
    # 2c:
    print('2c:\nThe value k = {} leads to the best clustering'.format\
        (best_k))
    # 2d:
    print('2d(1):\nThe best set of clusters is:\n{}\nIts W value is {:.3f}.'\
        .format(best_clusters_k[0], best_clusters_k[1]))
    print('2d(2):\nThe best set of clusters is:\n{}\nIts W value is {:.3f}.'\
        .format(best_clusters_all[0], best_clusters_all[1]))

    # Question 8:
    # boottrap
    boottrap_k = k_optimization\
        (filename2, k_large1, rep = 10, resampling = 'boottrapping')
    average_iter_bt, stdev_iter_bt, kmeans_bt, centers_bt, w_list_bt, \
        average_w_bt, stdev_w_bt, best_w_bt = result_output\
            (filename2, k = boottrap_k, rep = 10, resampling = None)
    # jackknife
    jackknife_k = k_optimization\
        (filename2, k_large1, rep = 10, resampling = 'jackknifing')
    average_iter_jk, stdev_iter_jk, kmeans_jk, centers_jk, w_list_jk, \
        average_w_jk, stdev_w_jk, best_w_jk = result_output\
            (filename2, k = jackknife_k, rep = 10, resampling = None)
    # subsample
    subsample_k = k_optimization\
        (filename2, k_large1, rep = 10, resampling = 'subsampling')
    average_iter_ss, stdev_iter_ss, kmeans_ss, centers_ss, w_list_ss, \
        average_w_ss, stdev_w_ss, best_w_ss = result_output\
            (filename2, k = subsample_k, rep = 10, resampling = None)
    print('Answer to question 8:\n',\
    'The optimal k from boottrap is {}.\nThe clusters are:\n{}\n\
    The optimal k from jackknife is {}.\nThe clusters are:\n{}\n\
    The optimal k from subsample is {}.\nThe clusters are:\n{}\n'.format\
        (boottrap_k, best_w_bt, jackknife_k, best_w_jk, \
            subsample_k, best_w_ss))
