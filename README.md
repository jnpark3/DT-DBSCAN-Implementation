# DT-DBSCAN

## Background

DBSCAN is a density based clustering algorithm that works by creating clusters from initial seed points [1]. A seed point is considered a "core point" if it is locally dense, and "local density" is defined using the two input parameters *minPts* and *epsilon*. A set of closely connected core points are categorized into the same custer. This form of density based clustering is helpful for finding outliers and arbitrarily shaped clusters.

DT-DBSCAN (Delaunay Triangulated DBSCAN)is an iteration of DBSCAN that aims to improve the expected time complexity of the aogirthm for larger sets of data. DT-DBSCAN utilizes edges formed by the Delaunay triangulation of the initial seeds to efficiently search for points within a locality. This results in faster classification of core points and noise points, which improves the run-time of the algorithm. Within this algorithm, the Delaunay triangulation is computed by using a linear run time algorithm created by Jyrki Katajainen and Markku Koppinen [2].

  [1] https://en.wikipedia.org/wiki/DBSCAN
  [2] http://hjemmesider.diku.dk/~jyrki/Paper/KK88.pdf

## How to Use

**DTDBSCAN.java**
> Computes DBSCAN clusters by using Delaunay triangulation edges. Inital seeds come from *data.in*, and cluster results are produced in *cluster.out*

**data.in**
> Input for DTDBSCAN.java. First line denotes the number of points within the seed set. After the first line, every subsequent line denotes the coordinates for the input data set. The *n*-th point's x-coordinate lies on line 2*n*, and the y-coordinate lies on the line 2*n*+1.

**cluster.out**
> Outputs the cluster results created by DTDBSCAN.java. Points within each cluster are grouped together, and points within different cluster are divided by a line. If there are noise points, then the points between the first line and the first dividing line denote the noise points. If there are no noise points, then the first line will be a space.

## About

* This algorithm is a part of a submission to the Orange County Science and Engineering Fair
* This algorithm is a part of a submission to the California State Science Fair
* **Note that this current iteration only works for data sets with 2D points. Currently developing an expansion to higher dimensions.**
# DT-DBSCAN-Implementation
