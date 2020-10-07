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

        int j=0;
        runClustering(points, 190000, 30, 10);
            for(int i = 0; i < 4; i++){
                j += (runClustering(points, 190000, 30, 10));                       // Run DTDBSCAN Clustering algorithm
            }
            System.out.println(j/4.0);

        f.close();
    }

    /*
    * This method computes the DBSCAN clusters of the set inputPoints in O(n) expected time if inputPoints is semi-uniform in density.
    *   Parameter inputPoints - array of Points to be clustered
    *   Parameter N - integer for size of dataset (points 1 ... N from inputPoints will be clustered)
    *   Parameter epsilon - a constant used in DBSCAN to determine clusters
    *   Parameter minPts - a constant used in DBSCAN to determine clusters
    */

    public static int runClustering(Point[] inputPoints, int N, int epsilon, int minPts)
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

        ArrayList<ArrayList<Integer>> Adjacency = delaunayTriangulation(grid, points, N);

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

    // Delaunay Triangulation Algorithm in Linear Expected Time [Katajainen &
    // Kopinen 1988]

    public static ArrayList<ArrayList<Integer>> delaunayTriangulation(ArrayList<ArrayList<ArrayList<Integer>>> grid,
            Point[] Points, int N) { // Computes Katajainen Koppinen O(n) triangulation*

        if (Points.length == 0) { // Eliminate case where points set is empty
            ArrayList<ArrayList<Integer>> Adjacency = new ArrayList<ArrayList<Integer>>();
            return Adjacency;
        }

        int K = (int) Math.pow(4, (int) (Math.log(N) / Math.log(4)) + 1);
        int[] Morton = morton_order(K, (int) Math.sqrt(K));

        Triangulation[] Trigs = new Triangulation[K];

        ArrayList<ArrayList<Integer>> Adjacency = new ArrayList<ArrayList<Integer>>();

        for (int i = 0; i < N; i++) {
            Adjacency.add(new ArrayList<Integer>());
        }

        for (int i = 0; i < K; i++) {
            Trigs[i] = gridInsert(Adjacency, Points, grid.get(Morton[2 * i]).get(Morton[(2 * i) + 1]));
        }

        for (int i = (int) (K / 4); i > 0; i /= 4) { // Quad-tree merging of Triangulations
            Triangulation[] Trigs_temp = new Triangulation[i];
            for (int j = 0; j < Trigs_temp.length; j++) {
                Triangulation TBottom = mergeTrigs(Adjacency, Points,
                        hullCalc(Points, Trigs[4 * j], Trigs[(4 * j) + 1], false), Trigs[4 * j], Trigs[(4 * j) + 1],
                        false);
                int[] hull1 = hullCalc(Points, Trigs[(4 * j) + 2], Trigs[(4 * j) + 3], false);
                Triangulation TTop = mergeTrigs(Adjacency, Points, hull1, Trigs[(4 * j) + 2], Trigs[(4 * j) + 3],
                        false);
                int[] hull2 = hullCalc(Points, TTop, TBottom, true);
                Trigs_temp[j] = mergeTrigs(Adjacency, Points, hull2, TBottom, TTop, true);
            }
            Trigs = new Triangulation[i];
            for (int j = 0; j < Trigs_temp.length; j++)
                Trigs[j] = Trigs_temp[j];
        }

        return Adjacency;
    }

    public static Triangulation gridInsert(ArrayList<ArrayList<Integer>> Adjacency, Point[] orderedPoints,
            ArrayList<Integer> gridPoints) { // returns triangulation for grid [primitive]**

        Triangulation result = new Triangulation();

        if (gridPoints.size() == 0) { // resolve case that grid is empty
            result.min_index = result.max_index = result.right = result.left = result.top = result.bottom = -1;
            result.convex = new int[0];
            return result;
        }

        result.min_index = result.left = result.right = result.bottom = result.top = gridPoints.get(0); // calculate result.right/left/top/bottom
        result.max_index = gridPoints.get(gridPoints.size() - 1);

        for (int i : gridPoints) {
            if (orderedPoints[i].x > orderedPoints[result.right].x)
                result.right = i;
            if (orderedPoints[i].x < orderedPoints[result.left].x)
                result.left = i;
            if (orderedPoints[i].y > orderedPoints[result.top].y)
                result.top = i;
            if (orderedPoints[i].y < orderedPoints[result.bottom].y)
                result.bottom = i;
        }

        if (gridPoints.size() == 1) { // resolve case where grid has 1 point
            result.convex = new int[1];
            result.convex[0] = gridPoints.get(0);
            return result;
        }
        if (gridPoints.size() == 2) { // resolve case where grid has 2 points
            Adjacency.get(gridPoints.get(0)).add(gridPoints.get(1));
            Adjacency.get(gridPoints.get(1)).add(gridPoints.get(0));
            result.convex = new int[2];
            result.convex[0] = result.left;
            result.convex[1] = result.right;
            if (result.left == result.right) {
            }
            return result;
        }
        if (gridPoints.size() == 3) { // resolve case where grid has 3 points
            if (angle(orderedPoints[gridPoints.get(0)], orderedPoints[gridPoints.get(1)]) < angle(
                    orderedPoints[gridPoints.get(0)], orderedPoints[gridPoints.get(2)])) {
                Adjacency.get(gridPoints.get(0)).add(gridPoints.get(1));
                Adjacency.get(gridPoints.get(0)).add(gridPoints.get(2));
            } else {
                Adjacency.get(gridPoints.get(0)).add(gridPoints.get(2));
                Adjacency.get(gridPoints.get(0)).add(gridPoints.get(1));
            }
            if (angle(orderedPoints[gridPoints.get(1)], orderedPoints[gridPoints.get(0)]) < angle(
                    orderedPoints[gridPoints.get(1)], orderedPoints[gridPoints.get(2)])) {
                Adjacency.get(gridPoints.get(1)).add(gridPoints.get(0));
                Adjacency.get(gridPoints.get(1)).add(gridPoints.get(2));
            } else {
                Adjacency.get(gridPoints.get(1)).add(gridPoints.get(2));
                Adjacency.get(gridPoints.get(1)).add(gridPoints.get(0));
            }
            if (angle(orderedPoints[gridPoints.get(2)], orderedPoints[gridPoints.get(1)]) < angle(
                    orderedPoints[gridPoints.get(2)], orderedPoints[gridPoints.get(0)])) {
                Adjacency.get(gridPoints.get(2)).add(gridPoints.get(1));
                Adjacency.get(gridPoints.get(2)).add(gridPoints.get(0));
            } else {
                Adjacency.get(gridPoints.get(2)).add(gridPoints.get(0));
                Adjacency.get(gridPoints.get(2)).add(gridPoints.get(1));
            }
            result.convex = new int[3];
            result.convex[0] = result.left;
            result.convex[1] = Adjacency.get(result.left).get(0);
            result.convex[2] = Adjacency.get(result.left).get(1);
            return result;
        } else { // triangulate case where grid has >3 points (O(n^4) algorithm)
            for (int k = 0; k < gridPoints.size(); k++) {
                for (int l = k + 1; l < gridPoints.size(); l++) {
                    for (int m = l + 1; m < gridPoints.size(); m++) {
                        boolean condition = false;
                        for (int n = 0; n < gridPoints.size(); n++) {
                            if (inside(orderedPoints[gridPoints.get(n)], orderedPoints[gridPoints.get(k)],
                                    orderedPoints[gridPoints.get(l)], orderedPoints[gridPoints.get(m)])) {
                                condition = true;
                                break;
                            }
                        }
                        if (condition == false) {
                            Adjacency.get(gridPoints.get(k)).add(gridPoints.get(l));
                            Adjacency.get(gridPoints.get(l)).add(gridPoints.get(k));
                            Adjacency.get(gridPoints.get(m)).add(gridPoints.get(l));
                            Adjacency.get(gridPoints.get(l)).add(gridPoints.get(m));
                            Adjacency.get(gridPoints.get(k)).add(gridPoints.get(m));
                            Adjacency.get(gridPoints.get(m)).add(gridPoints.get(k));
                        }
                    }
                }
            }
        }
        for (int i = 0; i < gridPoints.size(); i++) { // re-order adjacency list for each point in clockwise order

            ArrayList<Integer> temp = new ArrayList<Integer>();
            ArrayList<Integer> query = new ArrayList<Integer>();

            for (int j = 0; j < Adjacency.get(gridPoints.get(i)).size(); j++)
                query.add(Adjacency.get(gridPoints.get(i)).get(j));

            while (query.size() > 0) {
                int min = query.get(0);
                int min_index = 0;
                for (int k = 0; k < query.size(); k++) {
                    if (angle(orderedPoints[gridPoints.get(i)], orderedPoints[min]) > angle(
                            orderedPoints[gridPoints.get(i)], orderedPoints[query.get(k)])) {
                        min = query.get(k);
                        min_index = k;
                    }
                }
                query.remove(min_index);
                if (temp.size() == 0 || temp.get(temp.size() - 1) != min)
                    temp.add(min);
            }
            Adjacency.get(gridPoints.get(i)).clear();
            for (int number : temp)
                Adjacency.get(gridPoints.get(i)).add(number);
        }

        ArrayList<Integer> convex_temp = new ArrayList<Integer>(); // generates convex for triangulation starting at result.left
        convex_temp.add(result.left);

        if (!Adjacency.get(result.left).isEmpty()) {
            int check_index = -1;
            int check = Adjacency.get(result.left).get(0);
            for (int i = 0; i < Adjacency.get(check).size(); i++)
                if (Adjacency.get(check).get(i) == result.left)
                    check_index = i;
            while (result.left != check) {
                convex_temp.add(check);
                check_index++;
                if (check_index == Adjacency.get(check).size())
                    check_index = 0;
                int pre_check = check;
                check = Adjacency.get(check).get(check_index);
                for (int i = 0; i < Adjacency.get(check).size(); i++)
                    if (pre_check == Adjacency.get(check).get(i))
                        check_index = i;
            }
        }
        result.convex = new int[convex_temp.size()];
        for (int i = 0; i < convex_temp.size(); i++)
            result.convex[i] = convex_temp.get(i);

        return result;
    }

    public static int[] hullCalc(Point[] Points, Triangulation T1, Triangulation T2, Boolean vertical) { // output the top and bottom convex hull edges **

        if (T1.right == -1 || T2.right == -1) { // resolve when Trig1 or Trig2 are empty
            int[] results = new int[4 + T2.convex.length + T1.convex.length];
            for (int i = 0; i < 4; i++)
                results[i] = -1;
            if (T1.right == -1 && T2.right == -1)
                return results;
            else if (T1.right == -1)
                for (int i = 4; i < 4 + T2.convex.length; i++)
                    results[i] = T2.convex[i - 4];
            else if (T1.left == -1)
                for (int i = 4; i < 4 + T1.convex.length; i++)
                    results[i] = T1.convex[i - 4];
            return results;
        }

        int[] results = new int[4];

        int point1; // computing result[1] and result[2]
        int p1_index = -1;
        int point2;
        int p2_index = -1;

        if (vertical) {
            point1 = T1.bottom;
            point2 = T2.top;
        } else {
            point1 = T1.right;
            point2 = T2.left;
        }

        for (int i = 0; i < T1.convex.length; i++)
            if (T1.convex[i] == point1)
                p1_index = i;
        for (int i = 0; i < T2.convex.length; i++)
            if (T2.convex[i] == point2)
                p2_index = i;

        boolean run = true;
        int c1_index;
        int c2_index;

        while (vertical && run) {
            c1_index = p1_index + 1;
            if (c1_index >= T1.convex.length)
                c1_index = 0;

            boolean correct = false;
            if ((Points[T1.convex[p1_index]].x == Points[T2.convex[p2_index]].x)) {
                if (Points[T1.convex[p1_index]].x > Points[T1.convex[c1_index]].x && T1.convex.length > 1)
                    correct = true;
            } else {
                double Y;
                double slope = slope(Points[T1.convex[p1_index]], Points[T2.convex[p2_index]], false); // Slope identically define twice
                Y = (slope * (Points[T1.convex[c1_index]].x - Points[T1.convex[p1_index]].x))
                        + Points[T1.convex[p1_index]].y;
                if (slope < 0 && Y > Points[T1.convex[c1_index]].y)
                    correct = true;
                if (slope > 0 && Y < Points[T1.convex[c1_index]].y)
                    correct = true;
            }
            if (correct && T1.convex.length > 1) {
                p1_index = c1_index;
            } else {
                c2_index = p2_index - 1;
                if (c2_index < 0)
                    c2_index = T2.convex.length - 1;
                if (Points[T1.convex[p1_index]].x == Points[T2.convex[p2_index]].x) {
                    if (Points[T1.convex[p1_index]].x > Points[T2.convex[c2_index]].x && T2.convex.length > 1)
                        correct = true;
                } else {
                    double Y;
                    double slope = slope(Points[T1.convex[p1_index]], Points[T2.convex[p2_index]], false); // Slope identically defined twice
                    Y = (slope * (Points[T2.convex[c2_index]].x - Points[T1.convex[p1_index]].x))
                            + Points[T1.convex[p1_index]].y;
                    if (slope < 0 && Y > Points[T2.convex[c2_index]].y)
                        correct = true;
                    if (slope > 0 && Y < Points[T2.convex[c2_index]].y)
                        correct = true;
                }
                if (correct && T2.convex.length > 1) {
                    p2_index = c2_index;
                } else {
                    results[0] = T1.convex[p1_index];
                    results[1] = T2.convex[p2_index];
                    run = false;
                }
            }
        }
        while (!vertical && run) {
            c1_index = p1_index + 1;
            if (c1_index >= T1.convex.length)
                c1_index = 0;

            double slope = slope(Points[T1.convex[p1_index]], Points[T2.convex[p2_index]], false);
            double Y = (slope * (Points[T1.convex[c1_index]].x - Points[T1.convex[p1_index]].x))
                    + Points[T1.convex[p1_index]].y;
            if (Y > Points[T1.convex[c1_index]].y && T1.convex.length > 1) {
                p1_index = c1_index;
            } else {
                c2_index = p2_index - 1;
                if (c2_index < 0)
                    c2_index = T2.convex.length - 1;
                Y = (slope * (Points[T2.convex[c2_index]].x - Points[T1.convex[p1_index]].x))
                        + Points[T1.convex[p1_index]].y;
                if (Y > Points[T2.convex[c2_index]].y && T2.convex.length > 1) {
                    p2_index = c2_index;
                } else {
                    results[0] = T1.convex[p1_index];
                    results[1] = T2.convex[p2_index];
                    run = false;
                }
            }
        }

        if (vertical) { // computing result[3] and result[4]
            point1 = T1.bottom;
            point2 = T2.top;
        } else {
            point1 = T1.right;
            point2 = T2.left;
        }
        p1_index = -1;
        p2_index = -1;

        for (int i = 0; i < T1.convex.length; i++)
            if (T1.convex[i] == point1)
                p1_index = i;

        for (int i = 0; i < T2.convex.length; i++)
            if (T2.convex[i] == point2)
                p2_index = i;

        run = true;
        while (vertical && run) {
            c1_index = p1_index - 1;
            if (c1_index < 0)
                c1_index = T1.convex.length - 1;

            boolean correct = false;
            if (Points[T1.convex[p1_index]].x == Points[T2.convex[p2_index]].x) {
                if (Points[T1.convex[p1_index]].x < Points[T1.convex[c1_index]].x && T1.convex.length > 1)
                    correct = true;
            } else {
                double slope = slope(Points[T1.convex[p1_index]], Points[T2.convex[p2_index]], false); 
                double Y = (slope * (Points[T1.convex[c1_index]].x - Points[T1.convex[p1_index]].x))
                        + Points[T1.convex[p1_index]].y;
                if (slope < 0 && Y < Points[T1.convex[c1_index]].y)
                    correct = true;
                if (slope > 0 && Y > Points[T1.convex[c1_index]].y)
                    correct = true;
            }
            if (correct && T1.convex.length > 1) {
                p1_index = c1_index;
            } else {
                c2_index = p2_index + 1;
                if (c2_index == T2.convex.length)
                    c2_index = 0;
                if (Points[T1.convex[p1_index]].x == Points[T2.convex[p2_index]].x) {
                    if (Points[T1.convex[p1_index]].x < Points[T2.convex[c2_index]].x && T2.convex.length > 1)
                        correct = true;
                } else {
                    double slope = slope(Points[T1.convex[p1_index]], Points[T2.convex[p2_index]], false); 
                    double Y = (slope * (Points[T2.convex[c2_index]].x - Points[T1.convex[p1_index]].x))
                            + Points[T1.convex[p1_index]].y;
                    if (slope < 0 && Y < Points[T2.convex[c2_index]].y)
                        correct = true;
                    if (slope > 0 && Y > Points[T2.convex[c2_index]].y)
                        correct = true;
                }
                if (correct && T2.convex.length > 1) {
                    p2_index = c2_index;
                } else {
                    results[2] = T1.convex[p1_index];
                    results[3] = T2.convex[p2_index];
                    run = false;
                }
            }
        }
        while (!vertical && run) {
            c1_index = p1_index - 1;
            if (c1_index < 0)
                c1_index = T1.convex.length - 1;

            double slope = slope(Points[T1.convex[p1_index]], Points[T2.convex[p2_index]], false);
            double Y = (slope * (Points[T1.convex[c1_index]].x - Points[T1.convex[p1_index]].x))
                    + Points[T1.convex[p1_index]].y;
            if (Y < Points[T1.convex[c1_index]].y && T1.convex.length > 1) {
                p1_index = c1_index;
            } else {
                c2_index = p2_index + 1;
                if (c2_index >= T2.convex.length)
                    c2_index = 0;
                Y = (slope * (Points[T2.convex[c2_index]].x - Points[T1.convex[p1_index]].x))
                        + Points[T1.convex[p1_index]].y;
                if (Y < Points[T2.convex[c2_index]].y && T2.convex.length > 1) {
                    p2_index = c2_index;
                } else {
                    results[2] = T1.convex[p1_index];
                    results[3] = T2.convex[p2_index];
                    run = false;
                }
            }
        }

        ArrayList<Integer> newConvex = new ArrayList<Integer>();

        if (T1.convex.length == 1 && T2.convex.length == 1) {
            if (Points[T1.convex[0]].x > Points[T2.convex[0]].x) {
                newConvex.add(T2.convex[0]);
                newConvex.add(T1.convex[0]);
            } else {
                newConvex.add(T1.convex[0]);
                newConvex.add(T2.convex[0]);
            }
        } else {

            int index1 = -1;
            int index2 = -1;
            int index3 = -1;

            if (Points[T1.left].x < Points[T2.left].x || !vertical) {

                for (int i = 0; i < T1.convex.length; i++) {
                    if (T1.convex[i] == T1.left)
                        index1 = i;
                    if (T1.convex[i] == results[0])
                        index3 = i;
                }
                for (int i = 0; i < T2.convex.length; i++) {
                    if (T2.convex[i] == results[3]) {
                        index2 = i;
                        break;
                    }
                }

                newConvex.add(T1.convex[index1]);
                while (T1.convex[index1] != results[2]) {
                    index1++;
                    if (index1 >= T1.convex.length)
                        index1 = 0;
                    newConvex.add(T1.convex[index1]);
                }
                newConvex.add(T2.convex[index2]);
                while (T2.convex[index2] != results[1]) {
                    index2++;
                    if (index2 >= T2.convex.length)
                        index2 = 0;
                    newConvex.add(T2.convex[index2]);
                }
                newConvex.add(T1.convex[index3]);
                while (T1.convex[index3] != T1.left) {
                    index3++;
                    if (index3 >= T1.convex.length)
                        index3 = 0;
                    newConvex.add(T1.convex[index3]);
                }
                newConvex.remove(newConvex.size() - 1);
            } else {
                for (int i = 0; i < T2.convex.length; i++) {
                    if (T2.convex[i] == T2.left)
                        index1 = i;
                    if (T2.convex[i] == results[3])
                        index3 = i;
                }
                for (int i = 0; i < T1.convex.length; i++) {
                    if (T1.convex[i] == results[0]) {
                        index2 = i;
                        break;
                    }
                }

                newConvex.add(T2.convex[index1]);
                while (T2.convex[index1] != results[1]) {
                    index1++;
                    if (index1 >= T2.convex.length)
                        index1 = 0;
                    newConvex.add(T2.convex[index1]);
                }
                newConvex.add(T1.convex[index2]);
                while (T1.convex[index2] != results[2]) {
                    index2++;
                    if (index2 >= T1.convex.length)
                        index2 = 0;
                    newConvex.add(T1.convex[index2]);
                }
                newConvex.add(T2.convex[index3]);
                while (T2.convex[index3] != T2.left) {
                    index3++;
                    if (index3 >= T2.convex.length)
                        index3 = 0;
                    newConvex.add(T2.convex[index3]);
                }
                newConvex.remove(newConvex.size() - 1);
            }
        }

        int[] return_Solution = new int[newConvex.size() + 4];

        for (int i = 0; i < 4; i++)
            return_Solution[i] = results[i];
        for (int i = 0; i < newConvex.size(); i++)
            return_Solution[i + 4] = newConvex.get(i);

        return return_Solution;
    }

    public static Triangulation mergeTrigs(ArrayList<ArrayList<Integer>> Adjacency, Point[] Points, int[] Hull,
            Triangulation T1, Triangulation T2, Boolean vertical) {

        if (T2.min_index == -1) // resolving cases where either triangulation is empty
            return T1;
        else if (T1.min_index == -1)
            return T2;

        Triangulation T3 = new Triangulation();

        T3.convex = new int[Hull.length - 4]; // setting solution.convex
        for (int i = 0; i < (Hull.length - 4); i++)
            T3.convex[i] = Hull[i + 4];

        T3.min_index = Math.min(T1.min_index, T2.min_index); // setting all other variables for solution
        T3.max_index = Math.max(T1.max_index, T2.max_index);

        if (Points[T1.left].x < Points[T2.left].x)
            T3.left = T1.left;
        else
            T3.left = T2.left;
        if (Points[T1.right].x > Points[T2.right].x)
            T3.right = T1.right;
        else
            T3.right = T2.right;
        if (Points[T1.top].y > Points[T2.top].y)
            T3.top = T1.top;
        else
            T3.top = T2.top;
        if (Points[T1.bottom].y < Points[T2.bottom].y)
            T3.bottom = T1.bottom;
        else
            T3.bottom = T2.bottom;

        if (T1.convex.length == 1 && T2.convex.length == 1) {
            T3.convex = new int[2];
            T3.convex[0] = T3.left;
            T3.convex[1] = T3.right;
            Adjacency.get(T1.convex[0]).add(T2.convex[0]);
            Adjacency.get(T2.convex[0]).add(T1.convex[0]);
            return T3;
        }

        int L = Hull[0]; // Initial Variables initialization
        int R = Hull[1];

        while (L != Hull[2] || R != Hull[3]) { // Edge connection loop

            boolean leftStay = false; // Initial variables
            boolean rightStay = false;
            int index_L = 0;
            int index_R = 0;

            if (Adjacency.get(L).size() == 0) // Insert L&R into Adjacency lists
                Adjacency.get(L).add(R);
            else {
                int P1 = Adjacency.get(L).get(0);
                int P2 = Adjacency.get(L).get(Adjacency.get(L).size() - 1);
                if (angle(Points[L], Points[R]) < angle(Points[L], Points[P1])) {
                    index_L = 0;
                    Adjacency.get(L).add(0, R);
                } else if (angle(Points[L], Points[R]) > angle(Points[L], Points[P2])) {
                    index_L = Adjacency.get(L).size();
                    Adjacency.get(L).add(R);
                } else {
                    for (int i = 1; i < Adjacency.get(L).size(); i++) {
                        if (angle(Points[L], Points[R]) <= angle(Points[L], Points[Adjacency.get(L).get(i)])) {
                            index_L = i;
                            Adjacency.get(L).add(i, R);
                            break;
                        }
                    }
                }
            }

            if (Adjacency.get(R).size() == 0)
                Adjacency.get(R).add(L);
            else {
                int P1 = Adjacency.get(R).get(0);
                int P2 = Adjacency.get(R).get(Adjacency.get(R).size() - 1);
                if (angle(Points[R], Points[L]) < angle(Points[R], Points[P1])) {
                    index_R = 0;
                    Adjacency.get(R).add(0, L);
                } else if (angle(Points[R], Points[L]) > angle(Points[R], Points[P2])) {
                    index_R = Adjacency.get(R).size();
                    Adjacency.get(R).add(L);
                } else {
                    for (int i = 1; i < Adjacency.get(R).size(); i++) {
                        if (angle(Points[R], Points[L]) <= angle(Points[R], Points[Adjacency.get(R).get(i)])) {
                            index_R = i;
                            Adjacency.get(R).add(i, L);
                            break;
                        }
                    }
                }
            }

            int index = index_R + 1; // Right Analysis Set up
            if (index >= Adjacency.get(R).size())
                index = 0;
            int R1 = Adjacency.get(R).get(index); // Right Analysis

            boolean condition1 = angle(Points[L], Points[R]) > angle(Points[L], Points[R1]);
            boolean condition2 = angle(Points[L], Points[R1]) > angle(Points[L], Points[R]) - 180;
            boolean condition3 = angle(Points[L], Points[R]) + 180 < angle(Points[L], Points[R1]);
            boolean conditionA;
            if (vertical)
                conditionA = T1.convex.length > 1;
            else
                conditionA = T2.convex.length > 1;

            if (((condition1 && condition2) || condition3) && conditionA) {

                int index2 = index + 1;
                if (index2 >= Adjacency.get(R).size())
                    index2 = 0;
                int R2 = Adjacency.get(R).get(index2);

                while (inside(Points[R2], Points[R1], Points[L], Points[R])) {

                    Adjacency.get(R1).remove(Adjacency.get(R1).indexOf(R));
                    Adjacency.get(R).remove(index);

                    if (index >= Adjacency.get(R).size())
                        index = 0;
                    if (index2 >= Adjacency.get(R).size())
                        index2 = 0;
                    R1 = R2;
                    R2 = Adjacency.get(R).get(index2);
                }
            } else
                leftStay = true;

            int index2 = index_L - 1; // Left Analysis Set-Up
            if (index2 < 0)
                index2 = Adjacency.get(L).size() - 1;
            int L1 = Adjacency.get(L).get(index2); // Left Analysis

            boolean condition4 = angle(Points[R], Points[L]) < angle(Points[R], Points[L1]);
            boolean condition5 = angle(Points[R], Points[L1]) - 180 < angle(Points[R], Points[L]);
            boolean condition6 = angle(Points[R], Points[L]) > 180 + angle(Points[R], Points[L1]);
            boolean conditionB;

            if (vertical)
                conditionB = T2.convex.length > 1;
            else
                conditionB = T1.convex.length > 1;

            if (((condition4 && condition5) || condition6) && conditionB) {

                int index3 = index2 - 1;
                if (index3 < 0)
                    index3 = Adjacency.get(L).size() - 1;
                int L2 = Adjacency.get(L).get(index3);

                while (inside(Points[L2], Points[L1], Points[R], Points[L])) {

                    Adjacency.get(L1).remove(Adjacency.get(L1).indexOf(L));
                    Adjacency.get(L).remove(index2);

                    index2--;
                    if (index2 == -1)
                        index2 = Adjacency.get(L).size() - 1;
                    index3--;
                    if (index3 == -1)
                        index3 = Adjacency.get(L).size() - 1;
                    L1 = L2;
                    L2 = Adjacency.get(L).get(index3);
                }
            } else
                rightStay = true;

            if (leftStay) // Update the values of R and L
                L = L1;
            else {
                if (rightStay)
                    R = R1;
                else {
                    if (inside(Points[L1], Points[L], Points[R], Points[R1]))
                        L = L1;
                    else
                        R = R1;
                }
            }
        }

        if (Adjacency.get(Hull[2]).size() == 0) // Insert tangent[2] into adjacency
            Adjacency.get(Hull[2]).add(Hull[3]);
        else {
            int P1 = Adjacency.get(Hull[2]).get(0);
            int P2 = Adjacency.get(Hull[2]).get(Adjacency.get(Hull[2]).size() - 1);
            if (angle(Points[Hull[2]], Points[Hull[3]]) < angle(Points[Hull[2]], Points[P1]))
                Adjacency.get(Hull[2]).add(0, Hull[3]);
            else if (angle(Points[Hull[2]], Points[Hull[3]]) > angle(Points[Hull[2]], Points[P2]))
                Adjacency.get(Hull[2]).add(Hull[3]);
            else {
                for (int i = 1; i < Adjacency.get(Hull[2]).size(); i++) {
                    if (angle(Points[Hull[2]], Points[Hull[3]]) <= angle(Points[Hull[2]],
                            Points[Adjacency.get(Hull[2]).get(i)])) {
                        Adjacency.get(Hull[2]).add(i, Hull[3]);
                        break;
                    }
                }
            }
        }

        if (Adjacency.get(Hull[3]).size() == 0) // Insert tangent[3] into adjacency
            Adjacency.get(Hull[3]).add(L);
        else {
            int P1 = Adjacency.get(Hull[3]).get(0);
            int P2 = Adjacency.get(Hull[3]).get(Adjacency.get(Hull[3]).size() - 1);
            if (angle(Points[Hull[3]], Points[Hull[2]]) < angle(Points[Hull[3]], Points[P1]))
                Adjacency.get(Hull[3]).add(0, Hull[2]);
            else if (angle(Points[Hull[3]], Points[Hull[2]]) > angle(Points[Hull[3]], Points[P2]))
                Adjacency.get(Hull[3]).add(Hull[2]);
            else {
                for (int i = 1; i < Adjacency.get(Hull[3]).size(); i++) {
                    if (angle(Points[Hull[3]], Points[Hull[2]]) <= angle(Points[Hull[3]],
                            Points[Adjacency.get(Hull[3]).get(i)])) {
                        Adjacency.get(Hull[3]).add(i, Hull[2]);
                        break;
                    }
                }
            }
        }

        return T3;
    }

    // Auxiliary functions

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

    public static double angle(Point a, Point b) { // calculates angle from a to b starting at positive y-axis and going
                                                   // clockwise[**]
        if (a.x == b.x) {
            if (a.y > b.y)
                return 180;
            else
                return 0;
        } else {
            double cosine = (180 / Math.PI)
                    * Math.acos((a.x - b.x) / (Math.sqrt(Math.pow((a.x - b.x), 2) + Math.pow((a.y - b.y), 2))));
            if (b.x >= a.x && b.y >= a.y)
                return cosine - 90;
            else if (b.x < a.x && b.y >= a.y)
                return cosine + 270;
            else
                return 270 - cosine;
        }
    }

    public static boolean inside(Point a, Point b, Point c, Point d) { // Check if point a is in circumcircle of triangle bcd[**]

        if (a.x == b.x && a.y == b.y)                                       // Resolve cases were two points are same
            return false; 
        if (a.x == c.x && a.y == c.y)
            return false;
        if (a.x == d.x && a.y == d.y)
            return false;
        if (b.x == c.x && b.y == c.y)
            return false;
        if (b.x == d.x && b.y == d.y)
            return false;
        if (c.x == d.x && c.y == d.y)
            return false;

        Point circumcenter = new Point();

        if (slope(b, c, true) == slope(c, d, true)) { // Resolve cases were bcd is a degenerate triangle
            return true;
        } else if (b.x == c.x) { // Solve for circumcenter if b.x=c.x to avoid divide by 0
            circumcenter.y = (b.y + c.y) / 2;
            circumcenter.x = -slope(d, c, false) * (circumcenter.y - ((d.y + c.y) / 2)) + ((d.x + c.x) / 2);
        } else if (c.x == d.x) { // Solve for circumcenter if c.x=d.x to avoid divide by 0
            circumcenter.y = (c.y + d.y) / 2;
            circumcenter.x = -slope(b, c, false) * (circumcenter.y - ((b.y + c.y) / 2)) + ((b.x + c.x) / 2);
        } else { // Solve for circumcenter for all other cases
            circumcenter.x = (-slope(d, c, true) * (c.x + d.x) + slope(d, b, true) * (b.x + d.x) + (c.y - b.y))
                    / (2 * (slope(d, b, true) - slope(d, c, true)));
            circumcenter.y = slope(d, c, true) * (circumcenter.x - ((d.x + c.x) / 2)) + ((d.y + c.y) / 2);
        }

        if (dist(b, circumcenter) > dist(a, circumcenter)){ // Check if a is within the circumcircle
            return true;
        }
        return false;
    }

    public static double dist(Point a, Point b) { // Return euclidian distance between Points a and b
        return Math.pow(Math.pow(a.x - b.x, 2) + Math.pow(a.y - b.y, 2), 0.5);
    }

    public static double slope(Point a, Point b, boolean inverse) { // Return slope line ab
        if ((a.x == b.x && !inverse) || (a.y == b.y && inverse)) {
            System.out.println(" NOTE: slope value approximated to double.max");
            return Double.MAX_VALUE;
        }
        if (inverse){
            return (-(a.x - b.x) / (a.y - b.y));
        }
        else{
            return ((a.y - b.y) / (a.x - b.x));
        }
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