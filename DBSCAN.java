import java.util.*;
import java.io.*;

// O(n^2) DBSCAN Implementation for result comparison

public class DBSCAN {

    public static void main(String[] args) throws Exception {
        BufferedReader f = new BufferedReader(new FileReader( "dataGenerator.in" ));
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter( "clusterDBSCAN.out" )));
        int N = Integer.parseInt(f.readLine());
        double[] v = new double[2 * N];

        for (int i = 0; i < N; i++) {
            v[2 * i] = Double.parseDouble(f.readLine());
            v[2 * i + 1] = Double.parseDouble(f.readLine());
        }

        for(N = 0; N < 50000; N += 1000){
            int p = 0;
            int start, finish;
            start = (int) System.currentTimeMillis();     
            
            for(int m = 0; m < 10; m++){
                int[] label = new int[N];
                int C = 0;                                                /* Cluster counter */
                int minPts = 10;
                int eps = 30;
                    for (int i = 0; i < N; i++){
                        if(label[i] != 0) continue;
                        LinkedList<Integer> Neigh = RangeQuery(v, N, i, eps);
                        if (Neigh.size() < minPts) {     
                            label[i] = -1;
                            continue;
                        }
                        C++;
                        label[i] = C;
                        while(Neigh.size() > 0){
                            int a = Neigh.remove();
                            if (a == i) continue;
                            if (label[a] == -1) label[a] = C;
                            if (label[a] != 0) continue;
                            label[a] = C;
                            LinkedList<Integer> Neigh2 = RangeQuery(v, N, a, eps);
                            for(int b : Neigh2){
                                Neigh.add(b);
                            }
                        }
                    }   
                
            }
            finish = (int) System.currentTimeMillis(); // Measure Time Finish
            System.out.println((finish-start)/10.0);
        }    
            // int max = 0;
            // for(int i = 0; i < label.length; i++){
            //     if(label[i] > max)
            //         max = label[i];
            // }

            // for( int i = 0; i <= max; i++){
            //     for(int j = 0; j < label.length; j++){
            //         if(label[j] == i){
            //             out.println( "(" + v[2 * j] + "," + v[2 * j + 1] + ")" );    
            //         }
            //     }
            //     out.println();
            //     out.println();
            //     out.println();
            // }

            // out.println();
            // out.println();
            // out.println();

            // for(int j = 0; j < label.length; j++){
            //     if(label[j] == -1){
            //         out.println( "(" + v[2 * j] + "," + v[2 * j + 1] + ")" );    
            //     }
            // }

            /*
            for( int i = 0; i < label.length; i++){
                out.print( v[2 * i] + "," + v[2 * i + 1] + "," );
                out.println( label[i] );
            }
            */

        out.close();

        out.close();
        f.close();
    }

    public static LinkedList<Integer> RangeQuery(double[] v, int N, int i, int eps){
        LinkedList<Integer> solution = new LinkedList<Integer>();
        for(int j = 0; j < N; j++){
            if(Math.pow(v[2*i]-v[2*j], 2) + Math.pow(v[2*i+1]-v[2*j+1], 2) < Math.pow(eps, 2))
                solution.add(j);
        }
        return solution;
    }
    
}
