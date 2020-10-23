# DT-DBSCAN

## Background

DBSCAN is a density based clustering algorithm that works by creating clusters from initial points [1]. A point is considered a "core point" if it is locally dense, and "local density" is defined using the two input parameters *minPts* and *epsilon*. A set of closely connected core points are categorized into the same custer. This form of density based clustering is helpful for finding outliers and arbitrarily shaped clusters.

DT-DBSCAN (Delaunay Triangulated DBSCAN)is an iteration of DBSCAN that aims to improve the expected time complexity of the aogirthm for larger sets of data. DT-DBSCAN utilizes edges formed by the Delaunay triangulation of the initial seeds to efficiently search for points within a locality. This results in faster classification of core points and noise points, which improves the run-time of the algorithm. Within this algorithm, the Delaunay triangulation is computed by using a linear run time algorithm created by Jyrki Katajainen and Markku Koppinen [2].

  [1] https://en.wikipedia.org/wiki/DBSCAN
  [2] http://hjemmesider.diku.dk/~jyrki/Paper/KK88.pdf

## How to Use

**DTDBSCAN.java**
> Computes DBSCAN clusters in O(n) expected time.
> There are two ways to input points into the clustering algorithm. The first option is to generate random datasets using the built-in method dataGenerator() (Default; used during experiment). The second option is to input points into data.in (see below).
> The variables epsilon and minPts may be changed by altering the variables at the top of DTDBSCAN.main()
> Cluster results are output to data.out (see below).

**data.in**
> Input for DTDBSCAN.java. First line denotes the number of points within the seed set. After the first line, every subsequent line denotes the coordinates for the input data set. The *n*-th point's x-coordinate lies on line 2*n*, and the y-coordinate lies on the line 2*n*+1.

**data.out**
> Outputs the cluster results created by DTDBSCAN.java. Points within each cluster are grouped together, and points within different cluster are divided by an empty lines. If there are noise points, then the points between the first line and the first empty line denote the noise points. If there are no noise points, then the first line will be empty.

## About

* This algorithm is a part of a submission to PeerJ
* **Note that this current iteration only works for data sets with 2D points. In addition, the algorithm operates if there are no two points with the same x and y coordinates.**
# DT-DBSCAN-Implementation
