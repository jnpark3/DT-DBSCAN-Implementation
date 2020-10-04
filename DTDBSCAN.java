import java.util.ArrayList;
import java.util.*;
import java.io.*;

public class DTDBSCAN {

    static double[] coordinates;
    static ArrayList<ArrayList<Integer>> Repeats;
    static ArrayList<ArrayList<Integer>> Adjacency;
    
    public static class Triangulation {

        int min_index;
        int max_index;
        int[] convex;
        int left;
        int right;
        int bottom;
        int top;

        public Triangulation() {
        }
    }

    public static void clusterOut( double[] coords, int[] labels, int size, String fname ) throws Exception {
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter( fname )));

        int max = 0;
        for(int i = 0; i < labels.length; i++){
            if(labels[i] > max)
                max = labels[i];
        }

        for(int i = 0; i <= max; i++){
            for(int j = 0; j < labels.length; j++){
                if(labels[j] == i){
                    out.println( "(" + coords[2 * j] + "," + coords[2 * j + 1] + ")" );
                }
            }
            out.println();
            out.println();
            out.println();
            out.println();
        }

        out.close();
    }

    public static void main(String[] args) throws Exception {

        long start, finish;  

        BufferedReader f = new BufferedReader(new FileReader( "data.in" ));
        int N = Integer.parseInt(f.readLine());
        double[] vertices = new double[2 * N];

        for (int i = 0; i < N; i++) {
            vertices[2 * i] = Double.parseDouble(f.readLine());
            vertices[2 * i + 1] = Double.parseDouble(f.readLine());
        }

        f.close();

        start = System.currentTimeMillis();                             //Measure Time Start

        int[] labels = runClustering( vertices, N, 30, 10, false );

        finish = System.currentTimeMillis();                            //Measure Time Finish
        System.out.println( Long.toString( finish - start )) ;

        clusterOut( coordinates, labels, N, "clusterDTDBSCAN.out" );        
    }

    public static int[] runClustering( double[] vertices, int N, int EPSILON, int MINPTS, Boolean verbose ) throws Exception {

        delaunayTriangulation( vertices, N );

        int[] identities = new int[N];
        ArrayList<ArrayList<Integer>> neighborhoodSet = new ArrayList<ArrayList<Integer>>();

        for(int i = 0; i < N; i++){
            ArrayList<Integer> pointSet = classifyCorePoints(i, EPSILON, MINPTS);
            neighborhoodSet.add(pointSet);
            if(pointSet.size() > MINPTS) {
                identities[i] = 0;
            }
            else {
                identities[i] = 2;
            }
        }

        for(int i = 0; i < N; i++){
            if(identities[i] == 2) {
                ArrayList<Integer> coreNeighbors = new ArrayList<Integer>();
                for(int j : neighborhoodSet.get(i)){
                    if(identities[j] == 0){
                        identities[i] = 1;
                        coreNeighbors.add(j);
                    }
                }
                for(int j : coreNeighbors){
                    for(int k : coreNeighbors){
                        if(j > k){
                            Adjacency.get(j).add(k);
                            Adjacency.get(k).add(j);
                        }
                    }   
                }
            }
        }

        int[] labels = new int[N];
        int cluster = 1;
        for(int i = 0; i < N; i++){
            if(labels[i] != 0 || identities[i] != 0)
                continue;
            Queue<Integer> queue = new PriorityQueue<Integer> ();
            queue.add(i);
            labels[i] = cluster;
            while(queue.size() > 0){
                int v = queue.remove();
                for(int q : Adjacency.get(v)){
                    if(Math.pow(coordinates[2*v]-coordinates[2*q],2) + Math.pow(coordinates[(2*v)+1]-coordinates[(2*q)+1],2) < Math.pow(EPSILON, 2) && labels[q] == 0){
                        labels[q] = cluster;
                        if(identities[v] == 0);
                            queue.add(q);
                    }
                }

            }
            cluster++;
        }

        return labels;
        
    }

    public static void delaunayTriangulation(double[] oldCoords, int N){
        double xMin; double xMax; double yMin; double yMax;
        if(oldCoords.length == 0){
            return;
        }
        else{
            xMin = oldCoords[0];
            xMax = oldCoords[0];
            yMin = oldCoords[1];
            yMax = oldCoords[1];
            for(int i = 1; i < N; i++){
                if(oldCoords[2 * i] > xMax)
                    xMax = oldCoords[2 * i];
                if(oldCoords[2 * i] < xMin)
                    xMin = oldCoords[2 * i];
                if(oldCoords[2 * i + 1] > yMax)
                    yMax = oldCoords[2 * i + 1];
                if(oldCoords[2 * i + 1] < yMin)
                    yMin = oldCoords[2 * i + 1];
            }
        }

        int K = (int) Math.pow(4, (int) (Math.log(N) / Math.log(4)) + 1);

        double gridWidth = Math.max((yMax - yMin), (xMax - xMin));

        int gridSize = (int) Math.sqrt(K);

        ArrayList<ArrayList<ArrayList<Integer>>> grid = new ArrayList<ArrayList<ArrayList<Integer>>>();

        for (int i = 0; i < gridSize; i++) {
            ArrayList<ArrayList<Integer>> subGrid = new ArrayList<ArrayList<Integer>>();
            grid.add(subGrid);
        }

        for (int i = 0; i < gridSize; i++) {
            for (int j = 0; j < gridSize; j++) {
                ArrayList<Integer> subGrid = new ArrayList<Integer>();
                grid.get(i).add(subGrid);
            }
        }

        double kTemp = gridWidth / gridSize;

        for (int i = 0; i < N; i++) {

            int x = Math.min((int) ((oldCoords[2 * i] - xMin) / kTemp), gridSize - 1);
            int y = Math.min((int) ((oldCoords[2 * i + 1] - yMin) / kTemp), gridSize - 1);

            grid.get(x).get(y).add(i);
        }

        int[] morton_ordering = morton_order(K, gridSize);

        coordinates = new double[2 * N];

        int new_index = 0;

        for (int i = 0; i < K; i++) {
            for (int j = 0; j < grid.get(morton_ordering[2 * i]).get(morton_ordering[(2 * i) + 1]).size(); j++) {
                coordinates[2 * new_index] = oldCoords[2 * grid.get(morton_ordering[2 * i]).get(morton_ordering[(2 * i) + 1]).get(j)];
                coordinates[(2 * new_index) + 1] = oldCoords[(2 * grid.get(morton_ordering[2 * i]).get(morton_ordering[(2 * i) + 1]).get(j)) + 1];

                grid.get(morton_ordering[2 * i]).get(morton_ordering[(2 * i) + 1]).set(j, new_index);
                new_index++;
            }
        }

        Triangulation[] triangulation_set = new Triangulation[K];

        Repeats = new ArrayList<ArrayList<Integer>>();
        Adjacency = new ArrayList<ArrayList<Integer>>();

        for (int i = 0; i < N; i++) {
            ArrayList<Integer> temporary = new ArrayList<Integer>();
            Adjacency.add(temporary);
        }
        for (int i = 0; i < N; i++) {
            ArrayList<Integer> temporary = new ArrayList<Integer>();
            Repeats.add(temporary);
        }

        for (int i = 0; i < K; i++) {
            triangulation_set[i] = gridInsert(coordinates,
                    grid.get(morton_ordering[2 * i]).get(morton_ordering[(2 * i) + 1]));
        }

        for (int i = (int) (K / 4); i > 0; i /= 4) {
            Triangulation[] temp_triangulation_set = new Triangulation[i];
            for (int j = 0; j < temp_triangulation_set.length; j++) {
                int[] hull1 = hullCalc(coordinates, triangulation_set[4 * j], triangulation_set[(4 * j) + 1], false);
                Triangulation trigBottom = mergeTrigs(coordinates, hull1, triangulation_set[4 * j], triangulation_set[(4 * j) + 1], false);
                int[] hull2 = hullCalc(coordinates, triangulation_set[(4 * j) + 2], triangulation_set[(4 * j) + 3], false);
                Triangulation trigTop = mergeTrigs(coordinates, hull2, triangulation_set[(4 * j) + 2], triangulation_set[(4 * j) + 3], false);
                int[] hull3 = hullCalc(coordinates, trigTop, trigBottom, true);
                temp_triangulation_set[j] = mergeTrigs(coordinates, hull3, trigBottom, trigTop, true);
            }

            triangulation_set = new Triangulation[i];
            for(int j = 0; j < temp_triangulation_set.length; j++){
                triangulation_set[j] = temp_triangulation_set[j];
            }
        }
    }

    public static Triangulation gridInsert(double[] coords, ArrayList<Integer> pre_points) {

        Triangulation result = new Triangulation();

        if (pre_points.size() == 0) {
            result.min_index = result.max_index = result.right = result.left = result.top = result.bottom = -1;
            return result;
        }

        result.min_index = result.left = result.right = result.bottom = result.top = pre_points.get(0);
        result.max_index = pre_points.get(pre_points.size() - 1);

        for (int i = 0; i < pre_points.size(); i++) {
            if (coords[2 * pre_points.get(i)] > coords[2 * result.right])
                result.right = pre_points.get(i);
            if (coords[2 * pre_points.get(i)] < coords[2 * result.left])
                result.left = pre_points.get(i);
            if ((coords[(2 * pre_points.get(i)) + 1] < coords[(2 * result.bottom) + 1]))
                result.bottom = pre_points.get(i);
            if (coords[(2 * pre_points.get(i)) + 1] > coords[(2 * result.top) + 1])
                result.top = pre_points.get(i);
        }

        for (int i = 0; i < pre_points.size(); i++) {
            for (int j = i + 1; j < pre_points.size(); j++) {
                if (coords[2 * pre_points.get(i)] == coords[2 * pre_points.get(j)]
                        && coords[(2 * pre_points.get(i)) + 1] == coords[(2 * pre_points.get(j)) + 1]) {
                    Repeats.get(pre_points.get(i)).add(pre_points.get(j));
                    pre_points.remove(j);
                    j--;
                }
            }
        }

        if (pre_points.size() == 1) {
            result.convex = new int[1];
            result.convex[0] = pre_points.get(0);
            return result;
        }
        if (pre_points.size() == 2) {
            Adjacency.get(pre_points.get(0)).add(pre_points.get(1));
            Adjacency.get(pre_points.get(1)).add(pre_points.get(0));
            result.convex = new int[2];
            result.convex[0] = result.left;
            result.convex[1] = result.right;
            if(result.left == result.right){
            }
            return result;
        }
        if (pre_points.size() == 3) {
            if (angle(coords[2 * pre_points.get(0)], coords[(2 * pre_points.get(0)) + 1], coords[2 * pre_points.get(1)],
                    coords[(2 * pre_points.get(1)) + 1]) < angle(coords[2 * pre_points.get(0)],
                            coords[(2 * pre_points.get(0)) + 1], coords[2 * pre_points.get(2)],
                            coords[(2 * pre_points.get(2)) + 1])) {
                Adjacency.get(pre_points.get(0)).add(pre_points.get(1));
                Adjacency.get(pre_points.get(0)).add(pre_points.get(2));
            } else {
                Adjacency.get(pre_points.get(0)).add(pre_points.get(2));
                Adjacency.get(pre_points.get(0)).add(pre_points.get(1));
            }

            if (angle(coords[2 * pre_points.get(1)], coords[(2 * pre_points.get(1)) + 1], coords[2 * pre_points.get(0)],
                    coords[(2 * pre_points.get(0)) + 1]) < angle(coords[2 * pre_points.get(1)],
                            coords[(2 * pre_points.get(1)) + 1], coords[2 * pre_points.get(2)],
                            coords[(2 * pre_points.get(2)) + 1])) {
                Adjacency.get(pre_points.get(1)).add(pre_points.get(0));
                Adjacency.get(pre_points.get(1)).add(pre_points.get(2));
            } else {
                Adjacency.get(pre_points.get(1)).add(pre_points.get(2));
                Adjacency.get(pre_points.get(1)).add(pre_points.get(0));
            }

            if (angle(coords[2 * pre_points.get(2)], coords[(2 * pre_points.get(2)) + 1], coords[2 * pre_points.get(1)],
                    coords[(2 * pre_points.get(1)) + 1]) < angle(coords[2 * pre_points.get(2)],
                            coords[(2 * pre_points.get(2)) + 1], coords[2 * pre_points.get(0)],
                            coords[(2 * pre_points.get(0)) + 1])) {
                Adjacency.get(pre_points.get(2)).add(pre_points.get(1));
                Adjacency.get(pre_points.get(2)).add(pre_points.get(0));
            } else {
                Adjacency.get(pre_points.get(2)).add(pre_points.get(0));
                Adjacency.get(pre_points.get(2)).add(pre_points.get(1));
            }

            result.convex = new int[3];

            result.convex[0] = result.left;
            result.convex[1] = Adjacency.get(result.left).get(0);
            result.convex[2] = Adjacency.get(result.left).get(1);
            return result;
        }

        else if (pre_points.size() > 3) {
            for (int k = 0; k < pre_points.size(); k++) {
                for (int l = k + 1; l < pre_points.size(); l++) {
                    for (int m = l + 1; m < pre_points.size(); m++) {
                        boolean condition = false;
                        for(int n = 0; n < pre_points.size(); n++){
                            if(inside(coords[2 * pre_points.get(n)], coords[2 * pre_points.get(n) + 1], coords[2 * pre_points.get(m)], coords[2 * pre_points.get(m) + 1], coords[2 * pre_points.get(l)], coords[2 * pre_points.get(l) + 1], coords[2 * pre_points.get(k)], coords[2 * pre_points.get(k) + 1])){
                                condition = true;
                                break;
                            }
                        }
                        if (condition == false) {
                            Adjacency.get(pre_points.get(k)).add(pre_points.get(l));
                            Adjacency.get(pre_points.get(l)).add(pre_points.get(k));
                            Adjacency.get(pre_points.get(m)).add(pre_points.get(l));
                            Adjacency.get(pre_points.get(l)).add(pre_points.get(m));
                            Adjacency.get(pre_points.get(k)).add(pre_points.get(m));
                            Adjacency.get(pre_points.get(m)).add(pre_points.get(k));
                        }
                    }
                }
            }
        }

        for (int i = 0; i < pre_points.size(); i++) {

            ArrayList<Integer> correctAdjacency = new ArrayList<Integer>();
            ArrayList<Integer> query = new ArrayList<Integer>();

            for (int j = 0; j < Adjacency.get(pre_points.get(i)).size(); j++) {
                query.add(Adjacency.get(pre_points.get(i)).get(j));
            }

            while (query.size() > 0) {
                int min = query.get(0);
                int min_index = 0;
                for (int k = 0; k < query.size(); k++) {
                    if (angle(coords[2 * pre_points.get(i)], coords[(2 * pre_points.get(i)) + 1], coords[2 * min],
                            coords[(2 * min) + 1]) > angle(coords[2 * pre_points.get(i)],
                                    coords[(2 * pre_points.get(i)) + 1], coords[2 * query.get(k)],
                                    coords[(2 * query.get(k)) + 1])) {
                        min = query.get(k);
                        min_index = k;
                    }
                }
                query.remove(min_index);
                if (correctAdjacency.size() == 0) {
                    correctAdjacency.add(min);
                } else if (correctAdjacency.get(correctAdjacency.size() - 1) != min) {
                    correctAdjacency.add(min);
                }
            }

            Adjacency.get(pre_points.get(i)).clear();

            for (int number : correctAdjacency) {
                Adjacency.get(pre_points.get(i)).add(number);
            }

        }

        ArrayList<Integer> convex_temp = new ArrayList<Integer>();
        convex_temp.add(result.left);

        int initial = result.left;

        if ( !Adjacency.get( initial ).isEmpty() ) {
            int check_index = -322;
            int check = Adjacency.get( initial ).get(0);
            for (int i = 0; i < Adjacency.get(check).size(); i++) {
                if (Adjacency.get(check).get(i) == initial) {
                    check_index = i;
                }
            }

            while (initial != check) {
                convex_temp.add(check);
                check_index++;
                if (check_index == Adjacency.get(check).size()) {
                    check_index = 0;
                }
                int pre_check = check;
                check = Adjacency.get(check).get(check_index);
                for (int i = 0; i < Adjacency.get(check).size(); i++) {
                    if (pre_check == Adjacency.get(check).get(i)) {
                        check_index = i;
                    }
                }
            }
        }

        result.convex = new int[convex_temp.size()];
        for (int i = 0; i < convex_temp.size(); i++) {
            result.convex[i] = convex_temp.get(i);
        }

        return result;
    }

    public static int[] hullCalc(double[] coords, Triangulation Trig1, Triangulation Trig2, Boolean vertical) {

        int[] results = new int[4];

        if (Trig1.right == -1 || Trig2.right == -1) {
            if (Trig1.right == -1 && Trig2.right == -1) {
                for (int i = 0; i < 4; i++) {
                    results[i] = -2;
                }
                return results;
            } else if (Trig1.right == -1) {
                results = new int[4 + Trig2.convex.length];
                for (int i = 0; i < 4; i++) {
                    results[i] = -1;
                }
                for (int i = 4; i < 4 + Trig2.convex.length; i++) {
                    results[i] = Trig2.convex[i - 4];
                }
                return results;
            } else {
                results = new int[4 + Trig1.convex.length];
                for (int i = 0; i < 4; i++) {
                    results[i] = -1;
                }
                for (int i = 4; i < 4 + Trig1.convex.length; i++) {
                    results[i] = Trig1.convex[i - 4];
                }
                return results;
            }
        }

        int point1;
        int point1_index = -1;
        int point2;
        int point2_index = -1;

        if(vertical){
            point1 = Trig1.bottom;
            point2 = Trig2.top;
        } else{
            point1 = Trig1.right;
            point2 = Trig2.left;
        }

        for (int i = 0; i < Trig1.convex.length; i++) {
            if (Trig1.convex[i] == point1) {
                point1_index = i;
            }
        }
        for (int i = 0; i < Trig2.convex.length; i++) {
            if (Trig2.convex[i] == point2) {
                point2_index = i;
            }
        }

        boolean run = true;

        int check1_index;
        int check2_index;

        if(vertical){
            while (run) {
                check1_index = point1_index + 1;
                if (check1_index >= Trig1.convex.length)
                    check1_index = 0;

                boolean correct = false;
                double expectedY;
                if ((coords[2 * Trig1.convex[point1_index]] == coords[2 * Trig2.convex[point2_index]])) {
                    expectedY = coords[2 * Trig1.convex[point1_index]];
                    if (expectedY > coords[(2 * Trig1.convex[check1_index])] && Trig1.convex.length > 1) {
                        correct = true;
                    }
                } else {
                    double slope = ((coords[(2 * Trig1.convex[point1_index]) + 1]
                            - coords[(2 * Trig2.convex[point2_index]) + 1])
                            / (coords[2 * Trig1.convex[point1_index]] - coords[2 * Trig2.convex[point2_index]]));
                    expectedY = (slope * (coords[2 * Trig1.convex[check1_index]] - coords[2 * Trig1.convex[point1_index]]))
                            + coords[(2 * Trig1.convex[point1_index]) + 1];
                    if (slope < 0 && expectedY > coords[(2 * Trig1.convex[check1_index]) + 1])
                        correct = true;
                    if (slope > 0 && expectedY < coords[(2 * Trig1.convex[check1_index]) + 1])
                        correct = true;
                }
                if (correct && Trig1.convex.length > 1) {
                    point1_index = check1_index;
                } else {
                    check2_index = point2_index - 1;
                    if (check2_index < 0)
                        check2_index = Trig2.convex.length - 1;

                    double expectedY2;
                    if ((coords[2 * Trig1.convex[point1_index]] == coords[2 * Trig2.convex[point2_index]])) {
                        expectedY2 = coords[2 * Trig1.convex[point1_index]];
                        if (expectedY2 > coords[(2 * Trig2.convex[check2_index])] && Trig2.convex.length > 1) {
                            correct = true;
                        }
                    } else {
                        double slope = ((coords[(2 * Trig1.convex[point1_index]) + 1]
                                - coords[(2 * Trig2.convex[point2_index]) + 1])
                                / (coords[2 * Trig1.convex[point1_index]] - coords[2 * Trig2.convex[point2_index]]));
                        expectedY2 = (slope
                                * (coords[2 * Trig2.convex[check2_index]] - coords[2 * Trig1.convex[point1_index]]))
                                + coords[(2 * Trig1.convex[point1_index]) + 1];
                        if (slope < 0 && expectedY2 > coords[(2 * Trig2.convex[check2_index]) + 1])
                            correct = true;
                        if (slope > 0 && expectedY2 < coords[(2 * Trig2.convex[check2_index]) + 1])
                            correct = true;
                    }
                    if (correct && Trig2.convex.length > 1) {
                        point2_index = check2_index;
                    } else {
                        results[0] = Trig1.convex[point1_index];
                        results[1] = Trig2.convex[point2_index];
                        run = false;
                    }
                }
            }
        }
        else{
            while (run) {
                check1_index = point1_index + 1;
                if (check1_index >= Trig1.convex.length)
                    check1_index = 0;
    
                double slope = ((coords[(2 * Trig1.convex[point1_index]) + 1]
                        - coords[(2 * Trig2.convex[point2_index]) + 1])
                        / (coords[2 * Trig1.convex[point1_index]] - coords[2 * Trig2.convex[point2_index]]));
                double expectedY = (slope
                        * (coords[2 * Trig1.convex[check1_index]] - coords[2 * Trig1.convex[point1_index]]))
                        + coords[(2 * Trig1.convex[point1_index]) + 1];
                if (expectedY > coords[(2 * Trig1.convex[check1_index]) + 1] && Trig1.convex.length > 1) {
                    point1_index = check1_index;
                } else {
                    check2_index = point2_index - 1;
                    if (check2_index < 0)
                        check2_index = Trig2.convex.length - 1;
                    double expectedY2 = (slope
                            * (coords[2 * Trig2.convex[check2_index]] - coords[2 * Trig1.convex[point1_index]]))
                            + coords[(2 * Trig1.convex[point1_index]) + 1];
                    if (expectedY2 > coords[(2 * Trig2.convex[check2_index]) + 1] && Trig2.convex.length > 1) {
                        point2_index = check2_index;
                    } else {
                        results[0] = Trig1.convex[point1_index];
                        results[1] = Trig2.convex[point2_index];
                        run = false;
                    }
                }
            }
        }

        if(vertical){
            point1 = Trig1.bottom;
            point2 = Trig2.top;
        } else{
            point1 = Trig1.right;
            point2 = Trig2.left;
        }
        point1_index = -1;
        point2_index = -1;

        for (int i = 0; i < Trig1.convex.length; i++) {
            if (Trig1.convex[i] == point1) {
                point1_index = i;
            }
        }
        for (int i = 0; i < Trig2.convex.length; i++) {
            if (Trig2.convex[i] == point2) {
                point2_index = i;
            }
        }

        run = true;
        if(vertical){
            while (run) {

                check1_index = point1_index - 1;
                if (check1_index < 0)
                    check1_index = Trig1.convex.length - 1;

                boolean correct = false;
                double expectedY;
                if ((coords[2 * Trig1.convex[point1_index]] == coords[2 * Trig2.convex[point2_index]])) {
                    expectedY = coords[2 * Trig1.convex[point1_index]];
                    if (expectedY < coords[(2 * Trig1.convex[check1_index])] && Trig1.convex.length > 1) {
                        correct = true;
                    }
                } else {
                    double slope = ((coords[(2 * Trig1.convex[point1_index]) + 1]
                            - coords[(2 * Trig2.convex[point2_index]) + 1])
                            / (coords[2 * Trig1.convex[point1_index]] - coords[2 * Trig2.convex[point2_index]]));

                    expectedY = (slope * (coords[2 * Trig1.convex[check1_index]] - coords[2 * Trig1.convex[point1_index]]))
                            + coords[(2 * Trig1.convex[point1_index]) + 1];
                    if (slope < 0 && expectedY < coords[(2 * Trig1.convex[check1_index]) + 1])
                        correct = true;
                    if (slope > 0 && expectedY > coords[(2 * Trig1.convex[check1_index]) + 1])
                        correct = true;
                }
                if (correct && Trig1.convex.length > 1) {
                    point1_index = check1_index;
                } else {
                    check2_index = point2_index + 1;
                    if (check2_index == Trig2.convex.length)
                        check2_index = 0;

                    double expectedY2;
                    if ((coords[2 * Trig1.convex[point1_index]] == coords[2 * Trig2.convex[point2_index]])) {
                        expectedY2 = coords[2 * Trig1.convex[point1_index]];
                        if (expectedY2 < coords[(2 * Trig2.convex[check2_index])] && Trig2.convex.length > 1) {
                            correct = true;
                        }
                    } else {
                        double slope = ((coords[(2 * Trig1.convex[point1_index]) + 1]
                                - coords[(2 * Trig2.convex[point2_index]) + 1])
                                / (coords[2 * Trig1.convex[point1_index]] - coords[2 * Trig2.convex[point2_index]]));
                        expectedY2 = (slope
                                * (coords[2 * Trig2.convex[check2_index]] - coords[2 * Trig1.convex[point1_index]]))
                                + coords[(2 * Trig1.convex[point1_index]) + 1];
                        if (slope < 0 && expectedY2 < coords[(2 * Trig2.convex[check2_index]) + 1])
                            correct = true;
                        if (slope > 0 && expectedY2 > coords[(2 * Trig2.convex[check2_index]) + 1])
                            correct = true;
                    }
                    if (correct && Trig2.convex.length > 1) {
                        point2_index = check2_index;
                    } else {
                        results[2] = Trig1.convex[point1_index];
                        results[3] = Trig2.convex[point2_index];
                        run = false;
                    }
                }
            }
        }else{
            while (run) {
                check1_index = point1_index - 1;
                if (check1_index < 0)
                    check1_index = Trig1.convex.length - 1;
                double slope = ((coords[(2 * Trig1.convex[point1_index]) + 1]
                        - coords[(2 * Trig2.convex[point2_index]) + 1])
                        / (coords[2 * Trig1.convex[point1_index]] - coords[2 * Trig2.convex[point2_index]]));
                double expectedY = (slope
                        * (coords[2 * Trig1.convex[check1_index]] - coords[2 * Trig1.convex[point1_index]]))
                        + coords[(2 * Trig1.convex[point1_index]) + 1];
                if (expectedY < coords[(2 * Trig1.convex[check1_index]) + 1] && Trig1.convex.length > 1) {
                    point1_index = check1_index;
                } else {
                    check2_index = point2_index + 1;
                    if (check2_index >= Trig2.convex.length) {
                        check2_index = 0;
                    }
                    double expectedY2 = (slope
                            * (coords[2 * Trig2.convex[check2_index]] - coords[2 * Trig1.convex[point1_index]]))
                            + coords[(2 * Trig1.convex[point1_index]) + 1];
                    if (expectedY2 < coords[(2 * Trig2.convex[check2_index]) + 1] && Trig2.convex.length > 1) {
                        point2_index = check2_index;
                    } else {
                        results[2] = Trig1.convex[point1_index];
                        results[3] = Trig2.convex[point2_index];
                        run = false;
                    }
                }
            }
        }

        ArrayList<Integer> newConvex = new ArrayList<Integer>();

        if (Trig1.convex.length == 1 && Trig2.convex.length == 1) {
            if(coords[2 * Trig1.convex[0]] > coords[2 * Trig2.convex[0]]){
                newConvex.add(Trig2.convex[0]);
                newConvex.add(Trig1.convex[0]);
            }
            else{
                newConvex.add(Trig1.convex[0]);
                newConvex.add(Trig2.convex[0]);
            }
        } else {

            int starting_condition = 1;
            if(coords[2 * Trig1.left] > coords[2 * Trig2.left]){
                starting_condition = 2;
            }

            if(starting_condition == 1 || !vertical){
                int checking_index = -23;
                for(int i = 0; i < Trig1.convex.length; i++){
                    if(Trig1.convex[i] == Trig1.left){
                        checking_index = i;
                        break;
                    }
                }
                newConvex.add(Trig1.convex[checking_index]);
                while(Trig1.convex[checking_index] != results[2]){
                    checking_index++;
                    if(checking_index >= Trig1.convex.length)
                        checking_index = 0;
                    newConvex.add(Trig1.convex[checking_index]);
                }

                checking_index = -23;
                for(int i = 0; i < Trig2.convex.length; i++){
                    if(Trig2.convex[i] == results[3]){
                        checking_index = i;
                        break;
                    }
                }
                newConvex.add(Trig2.convex[checking_index]);
                while(Trig2.convex[checking_index] != results[1]){
                    checking_index++;
                    if(checking_index >= Trig2.convex.length)
                        checking_index = 0;
                    newConvex.add(Trig2.convex[checking_index]);
                }

                checking_index = -23;
                for(int i = 0; i < Trig1.convex.length; i++){
                    if(Trig1.convex[i] == results[0]){
                        checking_index = i;
                        break;
                    }
                }
                newConvex.add(Trig1.convex[checking_index]);
                while(Trig1.convex[checking_index] != Trig1.left){
                    checking_index++;
                    if(checking_index >= Trig1.convex.length)
                        checking_index = 0;
                    newConvex.add(Trig1.convex[checking_index]);
                }
                newConvex.remove(newConvex.size()-1);
            }
            else{
                int checking_index = -23;
                for(int i = 0; i < Trig2.convex.length; i++){
                    if(Trig2.convex[i] == Trig2.left){
                        checking_index = i;
                        break;
                    }
                }
                newConvex.add(Trig2.convex[checking_index]);
                while(Trig2.convex[checking_index] != results[1]){
                    checking_index++;
                    if(checking_index >= Trig2.convex.length)
                        checking_index = 0;
                    newConvex.add(Trig2.convex[checking_index]);
                }

                checking_index = -23;
                for(int i = 0; i < Trig1.convex.length; i++){
                    if(Trig1.convex[i] == results[0]){
                        checking_index = i;
                        break;
                    }
                }
                newConvex.add(Trig1.convex[checking_index]);
                while(Trig1.convex[checking_index] != results[2]){
                    checking_index++;
                    if(checking_index >= Trig1.convex.length)
                        checking_index = 0;
                    newConvex.add(Trig1.convex[checking_index]);
                }

                checking_index = -23;
                for(int i = 0; i < Trig2.convex.length; i++){
                    if(Trig2.convex[i] == results[3]){
                        checking_index = i;
                        break;
                    }
                }
                newConvex.add(Trig2.convex[checking_index]);
                while(Trig2.convex[checking_index] != Trig2.left){
                    checking_index++;
                    if(checking_index >= Trig2.convex.length)
                        checking_index = 0;
                    newConvex.add(Trig2.convex[checking_index]);
                }
                newConvex.remove(newConvex.size()-1);
            }
        }

        int[] return_Solution = new int[newConvex.size() + 4];

        return_Solution[0] = results[0];
        return_Solution[1] = results[1];
        return_Solution[2] = results[2];
        return_Solution[3] = results[3];

        for (int i = 0; i < newConvex.size(); i++) {
            return_Solution[i + 4] = newConvex.get(i);
        }

        return return_Solution;
    }

    public static Triangulation mergeTrigs(double[] coords, int[] tangents, Triangulation trig1, Triangulation trig2, Boolean vertical){

        Triangulation solution = new Triangulation();

        //Convex Initialization

        if(tangents.length <= 4){
            solution.convex = new int[0];
        }
        else{
            solution.convex = new int[tangents.length - 4];
            for(int i = 0; i < (tangents.length - 4) ; i++){
                solution.convex[i] = tangents[i + 4];
            }
        }

        //Min & Max index Initialization (Sending off small cases)

        if(trig1.min_index == -1 || trig2.min_index == -1){
            if(trig1.min_index == -1 && trig2.min_index == -1){
                solution.min_index = -1;
                solution.max_index = -1;
                solution.right = -1;
                solution.left = -1;
                solution.top = -1;
                solution.bottom = -1;
                return solution;
            }
            else if(trig1.min_index == -1){
                solution.min_index = trig2.min_index;
                solution.max_index = trig2.max_index;
                solution.right = trig2.right;
                solution.left = trig2.left;
                solution.top = trig2.top;
                solution.bottom = trig2.bottom;
                return solution;
            }
            else{
                solution.min_index = trig1.min_index;
                solution.max_index = trig1.max_index;
                solution.right = trig1.right;
                solution.left = trig1.left;
                solution.top = trig1.top;
                solution.bottom = trig1.bottom;
                return solution;
            }
        }
        else{
            solution.min_index = Math.min(trig1.min_index, trig2.min_index);
            solution.max_index = Math.max(trig1.max_index, trig2.max_index);

            if(coords[2 * trig1.left] < coords[2 * trig2.left])
                solution.left = trig1.left;
            else   
                solution.left = trig2.left;

            if(coords[2 * trig1.right] > coords[2 * trig2.right])
                solution.right = trig1.right;
            else   
                solution.right = trig2.right;

            if(coords[(2 * trig1.top) + 1] > coords[(2 * trig2.top) + 1])
                solution.top = trig1.top;
            else   
                solution.top = trig2.top;

            if(coords[(2 * trig1.bottom) + 1] < coords[(2 * trig2.bottom) + 1])
                solution.bottom = trig1.bottom;
            else   
                solution.bottom = trig2.bottom;

            if(trig1.convex.length == 1 && trig2.convex.length == 1){
                solution.convex = new int[2];
                solution.convex[0] = solution.left;
                solution.convex[1] = solution.right;
                Adjacency.get(trig1.convex[0]).add(trig2.convex[0]);
                Adjacency.get(trig2.convex[0]).add(trig1.convex[0]);
                return solution;
            }
        }

        //Initial Variables initialization

        int L = tangents[0];
        int R = tangents[1];

        //Edge connection loop

        while(L != tangents[2] || R != tangents[3]){

            //Initial Variables for While Loop

            double angle_L = angle(coords[2 * L], coords[(2 * L) + 1], coords[2 * R], coords[(2 * R) + 1]);
            double angle_R = angle(coords[2 * R], coords[(2 * R) + 1], coords[2 * L], coords[(2 * L) + 1]);

            boolean a = false;
            boolean b = false;

            int insert_index_L = 0;
            int insert_index_R = 0;

            //Insert into Adjacency lists

            if(Adjacency.get(L).size() == 0){
                Adjacency.get(L).add(R);
            }
            else{
                int firstPoint = Adjacency.get(L).get(0);
                int lastPoint = Adjacency.get(L).get(Adjacency.get(L).size() - 1);
                if(angle_L < angle(coords[2 * L], coords[(2 * L) + 1], coords[2 * firstPoint], coords[(2 * firstPoint) + 1])){
                    insert_index_L = 0;
                    Adjacency.get(L).add(0, R);
                }
                else if(angle_L > angle(coords[2 * L], coords[(2 * L) + 1], coords[2 * lastPoint], coords[(2 * lastPoint) + 1])){
                    insert_index_L = Adjacency.get(L).size();
                    Adjacency.get(L).add(R);
                }
                else{
                    for(int i = 1; i < Adjacency.get(L).size(); i++){
                        int point = Adjacency.get(L).get(i);
                        if(angle_L <= angle(coords[2 * L], coords[(2 * L) + 1], coords[2 * point], coords[(2 * point) + 1])){
                            insert_index_L = i;
                            Adjacency.get(L).add(i, R);
                            break;
                        }
                    }
                }
            }
            if(Adjacency.get(R).size() == 0){
                Adjacency.get(R).add(L);
            }
            else{
                int firstPoint = Adjacency.get(R).get(0);
                int lastPoint = Adjacency.get(R).get(Adjacency.get(R).size() - 1);
                if(angle_R < angle(coords[2 * R], coords[(2 * R) + 1], coords[2 * firstPoint], coords[(2 * firstPoint) + 1])){
                    insert_index_R = 0;
                    Adjacency.get(R).add(0, L);
                }
                else if(angle_R > angle(coords[2 * R], coords[(2 * R) + 1], coords[2 * lastPoint], coords[(2 * lastPoint) + 1])){
                    insert_index_R = Adjacency.get(R).size();
                    Adjacency.get(R).add(L);
                }
                else{
                    for(int i = 1; i < Adjacency.get(R).size(); i++){
                        int point = Adjacency.get(R).get(i);
                        if(angle_R <= angle(coords[2 * R], coords[(2 * R) + 1], coords[2 * point], coords[(2 * point) + 1])){
                            insert_index_R = i;
                            Adjacency.get(R).add(i, L);
                            break;
                        }
                    }
                }
            }

            //Right Analysis Set up

            int temp_index = insert_index_R + 1;
            if (temp_index >= Adjacency.get(R).size()) {
                temp_index = 0;
            }
            int R1 = Adjacency.get(R).get(temp_index);

            //Right Analysis

            boolean condition1 = angle_L > angle(coords[2 * L], coords[(2 * L) + 1], coords[2 * R1], coords[(2 * R1) + 1]);
            boolean condition2 = angle(coords[2 * L], coords[(2 * L) + 1], coords[2 * R1],coords[(2 * R1) + 1]) > angle_L - 180;
            boolean condition3 = angle_L + 180 < (angle(coords[2 * L], coords[(2 * L) + 1], coords[2 * R1], coords[(2 * R1) + 1]));
            boolean conditionA;
            if(vertical)
                conditionA = trig1.convex.length > 1;
            else
                conditionA = trig2.convex.length > 1;

            if (((condition1 && condition2) || condition3) && conditionA) {

                int temp_index_2 = temp_index + 1;
                if (temp_index_2 >= Adjacency.get(R).size()) {
                    temp_index_2 = 0;
                }
                int R2 = Adjacency.get(R).get(temp_index_2);

                while (inside(coords[2 * R2], coords[(2 * R2) + 1], coords[(2 * R1)], coords[(2 * R1) + 1], coords[2 * L], coords[(2 * L) + 1], coords[2 * R], coords[(2 * R) + 1])) {

                    Adjacency.get(R1).remove(Adjacency.get(R1).indexOf(R));
                    Adjacency.get(R).remove(temp_index);

                    if (temp_index >= Adjacency.get(R).size()) {
                        temp_index = 0;
                    }

                    if (temp_index_2 >= Adjacency.get(R).size()) {
                        temp_index_2 = 0;
                    }

                    R1 = R2;
                    R2 = Adjacency.get(R).get(temp_index_2);
                }
            } else {
                a = true;
            }

            //Left Analysis Set-Up

            int temp_index2 = insert_index_L - 1;
            if (temp_index2 < 0) {
                temp_index2 = Adjacency.get(L).size() - 1;
            }
            int L1 = Adjacency.get(L).get(temp_index2);

            //Left Analysis

            boolean condition4 = angle_R < angle(coords[2 * R], coords[(2 * R) + 1], coords[2 * L1], coords[(2 * L1) + 1]);
            boolean condition5 = angle(coords[2 * R], coords[(2 * R) + 1], coords[2 * L1], coords[(2 * L1) + 1]) - 180 < angle_R;
            boolean condition6 = angle_R > (180 + angle(coords[2 * R], coords[(2 * R) + 1], coords[2 * L1], coords[(2 * L1) + 1]));
            boolean conditionB;
            if(vertical)
                conditionB = trig2.convex.length > 1;
            else
                conditionB = trig1.convex.length > 1;

            if (((condition4 && condition5) || condition6) && conditionB) {

                int temp_index2_2 = temp_index2 - 1;
                if (temp_index2_2 < 0) {
                    temp_index2_2 = Adjacency.get(L).size() - 1;
                }
                int L2 = Adjacency.get(L).get(temp_index2_2);

                while (inside(coords[2 * L2], coords[(2 * L2) + 1], coords[(2 * L1)], coords[(2 * L1) + 1], coords[2 * R], coords[(2 * R) + 1], coords[2 * L], coords[(2 * L) + 1])) {

                    Adjacency.get(L1).remove(Adjacency.get(L1).indexOf(L)); 
                    Adjacency.get(L).remove(temp_index2);

                    temp_index2 = temp_index2 - 1;
                    if (temp_index2 == -1) {
                        temp_index2 = Adjacency.get(L).size() - 1;
                    }

                    temp_index2_2 = temp_index2_2 - 1;
                    if (temp_index2_2 == -1) {
                        temp_index2_2 = Adjacency.get(L).size() - 1;
                    }

                    L1 = L2;
                    L2 = Adjacency.get(L).get(temp_index2_2);
                }
            } else {
                b = true;
            }

            //Update the values of R and L

            if (a) {
                L = L1;
            } else {
                if (b) {
                    R = R1;
                } else {
                    if (inside(coords[2 * L1], coords[(2 * L1) + 1], coords[2 * L], coords[(2 * L) + 1], coords[2 * R],
                            coords[(2 * R) + 1], coords[2 * R1], coords[(2 * R1) + 1]) == false) {
                        R = R1;
                    } else {
                        L = L1;
                    }
                }
            }
        }

        //Insert Tangent

        double angle_tangent1 = angle(coords[2 * tangents[2]], coords[(2 * tangents[2]) + 1], coords[2 * tangents[3]],
                coords[(2 * tangents[3]) + 1]);
        double angle_tangent2 = angle(coords[2 * tangents[3]], coords[(2 * tangents[3]) + 1], coords[2 * tangents[2]],
                coords[(2 * tangents[2]) + 1]);

        if(Adjacency.get(tangents[2]).size() == 0){
            Adjacency.get(tangents[2]).add(tangents[3]);
        }
        else{
            int firstPoint = Adjacency.get(tangents[2]).get(0);
            int lastPoint = Adjacency.get(tangents[2]).get(Adjacency.get(tangents[2]).size() - 1);
            if(angle_tangent1 < angle(coords[2 * tangents[2]], coords[(2 * tangents[2]) + 1], coords[2 * firstPoint], coords[(2 * firstPoint) + 1])){
                Adjacency.get(tangents[2]).add(0, tangents[3]);
            }
            else if(angle_tangent1 > angle(coords[2 * tangents[2]], coords[(2 * tangents[2]) + 1], coords[2 * lastPoint], coords[(2 * lastPoint) + 1])){
                Adjacency.get(tangents[2]).add(tangents[3]);
            }
            else{
                for(int i = 1; i < Adjacency.get(tangents[2]).size(); i++){
                    int point = Adjacency.get(tangents[2]).get(i);
                    if(angle_tangent1 <= angle(coords[2 * tangents[2]], coords[(2 * tangents[2]) + 1], coords[2 * point], coords[(2 * point) + 1])){
                        Adjacency.get(tangents[2]).add(i, tangents[3]);
                        break;
                    }
                }
            }
        }
        if(Adjacency.get(tangents[3]).size() == 0){
            Adjacency.get(tangents[3]).add(L);
        }
        else{
            int firstPoint = Adjacency.get(tangents[3]).get(0);
            int lastPoint = Adjacency.get(tangents[3]).get(Adjacency.get(tangents[3]).size() - 1);
            if(angle_tangent2 < angle(coords[2 * tangents[3]], coords[(2 * tangents[3]) + 1], coords[2 * firstPoint], coords[(2 * firstPoint) + 1])){
                Adjacency.get(tangents[3]).add(0, tangents[2]);
            }
            else if(angle_tangent2 > angle(coords[2 * tangents[3]], coords[(2 * tangents[3]) + 1], coords[2 * lastPoint], coords[(2 * lastPoint) + 1])){
                Adjacency.get(tangents[3]).add(tangents[2]);
            }
            else{
                for(int i = 1; i < Adjacency.get(tangents[3]).size(); i++){
                    int point = Adjacency.get(tangents[3]).get(i);
                    if(angle_tangent2 <= angle(coords[2 * tangents[3]], coords[(2 * tangents[3]) + 1], coords[2 * point], coords[(2 * point) + 1])){
                        Adjacency.get(tangents[3]).add(i, tangents[2]);
                        break;
                    }
                }
            }
        }

        return solution;
    }

    public static ArrayList<Integer> classifyCorePoints(int a, double EPSILON, int MINPTS){

        ArrayList<Integer> query = new ArrayList<Integer>();
        ArrayList<Integer> finished = new ArrayList<Integer>();

        finished.add(a);

        for(int n : Adjacency.get(a))
            query.add(n);
        
        while( query.size() > 0 ){
            int check = query.get(query.size() - 1);
            query.remove(query.size()-1);

            Boolean condition = false;
            for(int num : finished)
                if(num == check) {
                    condition = true;
                    break;
                }
            
            if(condition)
                continue;

            double distance = Math.sqrt(Math.pow(coordinates[2 * a] - coordinates[2 * check], 2) + Math.pow(coordinates[(2 * a) + 1] - coordinates[(2 * check) + 1], 2));

            if(distance > EPSILON)
                continue;

            finished.add(check); 

            if(finished.size() > MINPTS){
                return finished;
            }

            for(int num : Adjacency.get(check)){
                query.add(num);
            }  
        }
        return finished;
    }
    
    public static double angle(double a, double b, double c, double d) {
        if (a == c && b == d) {
            return 0;
        } else if (a == c) {
            if (b > d)
                return 180;
            else
                return 0;
        } else {

            double cosine = (180 / Math.PI)
                    * Math.acos((a - c) / (Math.sqrt(Math.pow((a - c), 2) + Math.pow((b - d), 2))));
            if (c >= a && d >= b)
                return cosine - 90;
            else if (c < a && d >= b)
                return cosine + 270;
            else
                return 270 - cosine;

        }
    }

    public static boolean inside(double a, double b, double c, double d, double e, double f, double g, double h) {

        double intersectionX = 0;
        double intersectionY = 0;

        if(a == c && b == d){
            return false;
        }
        if(a == e && b == f){
            return false;
        }
        if(a == g && b == h){
            return false;
        }
        if(c == e && d == f){
            return false;
        }
        if(c == g && d == h){
            return false;
        }
        if(e == g && f == h){
            return false;
        }

        if (c == e && e == g) {
            return true; // add colinear test (return false if all four points are colinear)
        } else if (c == e || e == g) {
            if (c == e) {
                intersectionY = (h + f)/2;
                double slope = ((g-e)/(h-f));
                intersectionX = ((intersectionY - ((d + f)/2))/slope) + ((c + e)/2);
            } else if (e == g) {
                intersectionY = (d + f)/2;
                double slope =  ((c - e) / (d - f));
                intersectionX = ((intersectionY - ((h + f)/2))/slope) + ((g + e)/2);
            }
        } else {
            double slope1 = -((c - e) / (d - f));
            double slope2 = -((e - g) / (f - h));
            if (slope1 == slope2) {
                return true; // add colinear test (return false if all four points are colinear)
            } else {
                intersectionX = ((-(slope2) * ((e + g) / 2)) + ((f + h) / 2) + (slope1 * ((e + c) / 2)) - ((d + f) / 2))
                        / (slope1 - slope2);

                intersectionY = (slope1 * (intersectionX - ((c + e) / 2))) + ((d + f) / 2);
            }
        }

        double circumDistance = Math.sqrt(Math.pow(c - intersectionX, 2) + Math.pow(d - intersectionY, 2));
        double aDistance = Math.sqrt(Math.pow(a - intersectionX, 2) + Math.pow(b - intersectionY, 2));

        if (circumDistance > aDistance) {
            return true;
        } else {
            return false;
        }
    }

    public static int[] morton_order(int K, int G) {
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
                    if (k % 2 == 0) {
                        a = a - (int) Math.pow(2, (k / 2));
                    } else {
                        b = b - (int) Math.pow(2, ((k - 1) / 2));
                    }
                } else {
                    if (k % 2 == 0) {
                        a = a + (int) Math.pow(2, (k / 2));
                    } else {
                        b = b + (int) Math.pow(2, ((k - 1) / 2));
                    }
                    break;
                }

            }

        }
        return morton_ordering;
    }
}