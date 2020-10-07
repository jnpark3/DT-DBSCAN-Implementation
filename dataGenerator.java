import java.io.*;

public class dataGenerator {
    public static void main(String[] args) throws Exception {
        PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter( "dataGenerator.out" )));

        for(int i = 0; i < 200000; i++){
            out.println(Math.random()*1000);
            out.println(Math.random()*1000);
        }

        out.close();
    }
}
