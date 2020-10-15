import java.util.*;
import java.io.*;

public class DTDBSCAN {

    public static class Triangulation {

        int min_index;
        int max_index;
        int left;
        int right;
        int bottom;
        int top;
        int[] convex;

        public Triangulation() {
        }
    }

    public static class Point {

        double x;
        double y;

        public Point() {
        }
    }

    public static void main(String[] args) throws Exception {

        BufferedReader f = new BufferedReader(new FileReader("dataGenerator.in"));  // Read file data.in
        int N = Integer.parseInt(f.readLine());                                     // Set N as size of dataset
        Point[] points = new Point[N];

        for (int i = 0; i < N; i++) {                                               // Generate Point set using coordinate information from data.in
            Point P = new Point();
            P.x = Double.parseDouble(f.readLine());
            P.y = Double.parseDouble(f.readLine());
            points[i] = P;
        }

        runClustering(points, 1000, 10, 10);                                        // Run DTDBSCAN Clustering algorithm

        f.close();
    }

    /*
    * This method computes the DBSCAN clusters of the set inputPoints in O(n) expected time if inputPoints is semi-uniform in density.
    *   Parameter inputPoints - array of Points to be clustered
    *   Parameter N - integer for size of dataset (points 1 ... N from inputPoints will be clustered)
    *   Parameter epsilon - a constant used in DBSCAN to determine clusters
    *   Parameter minPts - a constant used in DBSCAN to determine clusters
    */

    public static int runClustering(Point[] inputPoints, int N, double epsilon, int minPts)
            throws Exception {                                                      

        int start, finish;
        start = (int) System.currentTimeMillis();                                   // Measure Time Start

        // Step 1: Preparation for Delaunay Triangulation:
        //      - Outlined in [Katajainen & Koppinen 1988]
        //      - Insert every point in inputPoints into an O(sqrt(n)) x O(sqrt(n)) square grid on R^2
        //      - To facilitate with later quad-tree merging process, inputPoint array is re-ordered 

        double xMin, xMax, yMin, yMax;

        if (inputPoints.length == 0) {
            return 0;
        }                                                                           

        xMin = xMax = inputPoints[0].x;
        yMin = yMax = inputPoints[0].y;
        for (int i = 1; i < N; i++) {
            if (inputPoints[i].x > xMax){
                xMax = inputPoints[i].x;
            }
            if (inputPoints[i].x < xMin){
                xMin = inputPoints[i].x;
            }
            if (inputPoints[i].y > yMax){
                yMax = inputPoints[i].y;
            }
            if (inputPoints[i].y < yMin){
                yMin = inputPoints[i].y;
            }
        }

        int gridSpaces = (int) Math.pow(4, (int) (Math.log(N) / Math.log(4)) + 1);
        int gridSize = (int) Math.sqrt(gridSpaces);                                 
        double width = Math.max((yMax - yMin), (xMax - xMin)) / gridSize;           

        ArrayList<ArrayList<ArrayList<Integer>>> grid = new ArrayList<ArrayList<ArrayList<Integer>>>(); 

        for (int i = 0; i < gridSize; i++){                                         
            grid.add(new ArrayList<ArrayList<Integer>>());
            for (int j = 0; j < gridSize; j++){
                grid.get(i).add(new ArrayList<Integer>());
            }
        }
        for (int i = 0; i < N; i++) {                                               
            int x = Math.min((int) ((inputPoints[i].x - xMin) / width), gridSize - 1);
            int y = Math.min((int) ((inputPoints[i].y - yMin) / width), gridSize - 1);
            grid.get(x).get(y).add(i);
        }

        int[] morton = morton_order(gridSpaces, gridSize);                          
        Point[] points = new Point[N];                                              
        int n = 0;                                                                  

        for (int i = 0; i < gridSpaces; i++) {                                      
            for (int j = 0; j < grid.get(morton[2 * i]).get(morton[(2 * i) + 1]).size(); j++) {
                points[n] = inputPoints[grid.get(morton[2 * i]).get(morton[(2 * i) + 1]).get(j)];
                grid.get(morton[2 * i]).get(morton[(2 * i) + 1]).set(j, n);
                n++;
            }
        }

        // Step 2: Linear Time Delaunay Triangulation:
        //      - Outlined in [Katajainen & Koppinen 1988]
        //      - Triangulates the points in each grid space
        //      - Merges the triangulations in the grid spaces in a quad-tree order to form completed Triangulation

        delaunayTriangulation process = new delaunayTriangulation();
        ArrayList<ArrayList<Integer>> Adjacency = process.triangulation(grid, points, N);

        // Step 3: Core Point Classification
        //      - Algorithm 1 in paper
        //      - Checks to see if there are minPts points within epsilon distance of each point by tree-searching through Delaunay Triangulation
        //      - Each point is labeled as 0 if it is a core point and 2 if it is a noise point
        //      - For noise points, the neighborhood is stored for later use

        int[] labels = new int[N];                                                 
        ArrayList<ArrayList<Integer>> neighborhoods = new ArrayList<ArrayList<Integer>>();  

        for (int i = 0; i < N; i++) {                                              
            ArrayList<Integer> locality = classifyCorePoints(points, Adjacency, i, epsilon, minPts);
            if (locality.size() > minPts) {
                labels[i] = 0;
                neighborhoods.add(new ArrayList<Integer>());
            } else {
                labels[i] = 2;
                neighborhoods.add(locality);
            }
        }

        // Step 4: Cluster Identification
        //      - Algorithm 2 in paper
        //      - Boundary points have their labels changed from 2 to 1
        //      - An edge is constructed for every pair of core points which contain 
        //      - Clusters are generated by floodfilling through core points in the graph

        for (int i = 0; i < N; i++) {
            if (labels[i] == 2) {
                ArrayList<Integer> coreNeighbors = new ArrayList<Integer>();
                for (int j : neighborhoods.get(i)) {
                    if (labels[j] == 0) {
                        labels[i] = 1;
                        coreNeighbors.add(j);
                    }
                }
                for (int j : coreNeighbors) {
                    for (int k : coreNeighbors) {
                        if (j > k) {
                            Adjacency.get(j).add(k);
                            Adjacency.get(k).add(j);
                        }
                    }
                }
            }
        }

        int[] clusters = new int[N];
        int cluster = 1;
        for (int i = 0; i < N; i++) {
            if (clusters[i] != 0 || labels[i] != 0)
                continue;
            Queue<Integer> queue = new PriorityQueue<Integer>();
            queue.add(i);
            clusters[i] = cluster;
            while (queue.size() > 0) {
                int v = queue.remove();
                for (int q : Adjacency.get(v)) {
                    if (dist(points[v], points[q]) < epsilon && clusters[q] == 0) {
                        clusters[q] = cluster;
                        if (labels[v] == 0);
                            queue.add(q);
                    }
                }

            }
            cluster++;
        }

        // End of Algorithm

        finish = (int) System.currentTimeMillis();                                  // Measure Time Finish

        //clusterOut(points, clusters, N, "clusterDTDBSCAN.out");                   // Outputs clusters to clusterDTDBSCAN.out

        return finish - start;
    }

    // Note: Annotations are incomplete beyond this point

    public static ArrayList<Integer> classifyCorePoints(Point[] Points, ArrayList<ArrayList<Integer>> Adjacency, int a,
            double epsilon, int minPts) {

        ArrayList<Integer> query = new ArrayList<Integer>();
        ArrayList<Integer> finished = new ArrayList<Integer>();

        finished.add(a);

        for (int n : Adjacency.get(a))
            query.add(n);

        while (query.size() > 0) {
            int check = query.get(query.size() - 1);
            query.remove(query.size() - 1);

            Boolean condition = false;
            for (int num : finished)
                if (num == check) {
                    condition = true;
                    break;
                }
            if (condition)
                continue;
            if (dist(Points[a], Points[check]) > epsilon)
                continue;
            finished.add(check);
            if (finished.size() > minPts)
                return finished;
            for (int num : Adjacency.get(check)) {
                query.add(num);
            }
        }
        return finished;
    }

    // Auxiliary functions

    public static Point[] dataDgenerator(int k) throws Exception {
        Point[] solution = new Point[k];
        double[] centers = new double[14];
        for(int i = 0; i < 14; i++){
            centers[i] = Math.random()*800 + 100;
        }
        for(int i = 0; i < k; i++){
            Point in = new Point();
            if(Math.random()<1){   // Set to 0 for r=0, 0.359 for r=2, 0.834 for r=10, and 0.965 for r=50
                int rand = (int)(7*Math.random());
                in.x = centers[2*rand] + Math.random()*200 - 100;
                in.y = centers[2*rand+1] + Math.random()*200 - 100;
            }
            else{
                in.x = Math.random()*1000;
                in.y = Math.random()*1000;
            }
            solution[i] = in;
        }
        return solution;
    }

    public static void clusterOut(Point[] Points, int[] clusters, int size, String fname) throws Exception { 
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fname)));

        int max = 0;
        for (int i = 0; i < clusters.length; i++)
            if (clusters[i] > max)
                max = clusters[i];

        for (int i = 0; i <= max; i++) {
            for (int j = 0; j < clusters.length; j++)
                if (clusters[j] == i)
                    out.println("(" + Points[j].x + "," + Points[j].y + ")");
            out.println();
            out.println();
            out.println();
            out.println();
            out.println();
            out.println();
        }

        out.close();
    }

    public static double dist(Point a, Point b) { // Return euclidian distance between Points a and b
        return Math.pow(Math.pow(a.x - b.x, 2) + Math.pow(a.y - b.y, 2), 0.5);
    }

    public static int[] morton_order(int K, int G) { // computes morton ordering of grid [*** No Change]
        int[] morton_ordering = new int[2 * K];
        int[] morton_number = new int[2 * G];

        int a = 0;
        int b = 0;

        for (int i = 0; i < K; i++) {
            morton_ordering[2 * i] = a;
            morton_ordering[(2 * i) + 1] = b;
            for (int k = 0; k < morton_number.length; k++) {
                morton_number[k]++;
                if (morton_number[k] > 1) {
                    morton_number[k] = 0;
                    if (k % 2 == 0)
                        a = a - (int) Math.pow(2, (k / 2));
                    else
                        b = b - (int) Math.pow(2, ((k - 1) / 2));
                } else {
                    if (k % 2 == 0)
                        a = a + (int) Math.pow(2, (k / 2));
                    else
                        b = b + (int) Math.pow(2, ((k - 1) / 2));
                    break;
                }
            }
        }
        return morton_ordering;
    }
}