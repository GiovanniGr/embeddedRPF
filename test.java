import Jama.Matrix;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;

public class test {
    public static double[][] readMatrixFromPython(String filename, int r, int c){
        double [][] mat = new double[r][c];
        try  
        {  
        //the file to be opened for reading  
        FileInputStream fis=new FileInputStream(filename);       
        Scanner sc=new Scanner(fis);    //file to be scanned  
        //returns true if there is another line to read  
        int l=0;
        while(sc.hasNextLine())  
        {
            String line = sc.nextLine();
            String[] num_list = line.split(",");
            mat[l++] = Arrays.stream(num_list).mapToDouble(Double::parseDouble).toArray();
            //System.out.println(sc.nextLine());      //returns the line that was skipped
        }  
        sc.close();     //closes the scanner  
        }  
        catch(IOException e){e.printStackTrace();}
        return mat;

    }

    public static ArrayList<Matrix> readCovLinearFromPython(String filename, int dim){
        ArrayList<Matrix> cov = new ArrayList<Matrix>();

        try  
        {  
        //the file to be opened for reading  
        FileInputStream fis=new FileInputStream(filename);       
        Scanner sc=new Scanner(fis);    //file to be scanned  
        //returns true if there is another line to read  
        int l=0;
        Matrix aux = new Matrix(dim, dim);
        int dim2 = dim*dim;
        while(sc.hasNextLine())  
        {
            String line = sc.nextLine();
            if(l%dim2==0 && l!=0){
                cov.add(aux.copy());
            }
            aux.set(l%dim2/dim, l%dim2%dim, Double.parseDouble(line));
            //firstPotCov[l/4][(l%4)/2][(l%4)%2] = Double.parseDouble(line);
            l++;
            //System.out.println(sc.nextLine());      //returns the line that was skipped
        }
        cov.add(aux.copy());
        sc.close();     //closes the scanner  
        }  
        catch(IOException e){e.printStackTrace();}
        return cov;

    }


    public static void writeMatrixToPython(double[][] mat, String filename){
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(filename));
      
            for (int i = 0; i < mat.length; i++) {
              for (int j = 0; j < mat[i].length; j++) {
                bw.write(mat[i][j] + ((j == mat[i].length-1) ? "" : ","));
              }
                bw.newLine();
            }
            bw.flush();
            bw.close();
          } catch (IOException e) {}
    }

    public static void testCov(double[][] mat){
        double[][] filtered_mat = SignalProcAux.filter_bandpass(mat, 4, 30, new int[]{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32});
        double[][] mat_for_covariance = new double[33][1000];
        for(int r=0;r<33;r++){
        for(int c=0;c<1000;c++){
            mat_for_covariance[r][c] = mat[r][c];
        }
        }

        //double[][] cov_matrix_test = MatrixAux.computeCovarianceMatrix(new double[][]{{2.3,5.1,8.9,9.2,1.9,10.34,9.45,7.67,3.87,2.12,12,0.11,12.4},{1.3,7.54,8.12,3.98,5.26,8.34,19.2,0.012,8.23,78.65,11.23,12.67,87.45},{4.5,7.8,12.45,93.135,0.001,0.932,9.45,9.1,2,4.65,14.765,62.098,0.156}}).getArray();
        Matrix cov = MatrixAux.computeCovarianceMatrix(mat_for_covariance);
        Matrix cov_normTrace = MatrixAux.normalizeCovMatrix(cov,"trace");
        Matrix cov_normCov = MatrixAux.normalizeCovMatrix(cov,"cov");

        writeMatrixToPython(mat,"test/testBaseJava.txt");
        writeMatrixToPython(filtered_mat,"test/testFilteredJava.txt");
        writeMatrixToPython(cov.getArray(),"test/testCovJava.txt");
        writeMatrixToPython(cov_normTrace.getArray(),"test/testCovNormTraceJava.txt");
        writeMatrixToPython(cov_normCov.getArray(),"test/testCovNormCovJava.txt");
    }

    public static void testDistance_Geod_Mean(){
        double[][] mat = new double[][]{{15.40964107,   8.5875238 ,   9.1229583 ,  12.67034506,
            6.07197335,   7.09535784}, {8.5875238 ,  70.23148676, -12.75052337,  31.62069215,
              67.88629983,  70.08624275}, {9.1229583 , -12.75052337,  17.39204523,   5.17331657,
                -15.34425627, -15.50582496}, {12.67034506,  31.62069215,   5.17331657,  39.94475374,
                  40.49988185,  43.30090519}, {6.07197335,  67.88629983, -15.34425627,  40.49988185,
                    77.47320063,  79.56611252}, {7.09535784,  70.08624275, -15.50582496,  43.30090519,
                      79.56611252,  84.76104473}};
          
          double[][] mat2 = new double[][]{{5.66821762e+02, 4.55513333e+02, 4.45601375e+02, 3.32766550e+02,
            2.16859391e+02, 7.51203320e+01},{4.55513333e+02, 4.36102150e+02, 2.78248376e+02, 2.19497727e+02,
              1.67506595e+02, 1.00176972e+02},{4.45601375e+02, 2.78248376e+02, 4.88420961e+02, 3.45040016e+02,
                1.94516275e+02, 4.34459884e-01},{3.32766550e+02, 2.19497727e+02, 3.45040016e+02, 2.48745558e+02,
                  1.42971860e+02, 9.80899376e+00},{2.16859391e+02, 1.67506595e+02, 1.94516275e+02, 1.42971860e+02,
                    9.38859854e+01, 2.28484556e+01},{7.51203320e+01, 1.00176972e+02, 4.34459884e-01, 9.80899376e+00,
                      2.28484556e+01, 4.26424943e+01}};
      
          double[][] mat1 = new double[][]{{548.44055506, 469.91721369, 226.61698296, 152.8056371 ,
            101.86262846,  74.82005621},{469.91721369, 565.24965173,  47.31755651,  92.94547321,
              153.60768909, 206.62431451},{226.61698296,  47.31755651, 281.80538437, 152.07256519,
                16.66114356, -85.35062017},{152.8056371 ,  92.94547321, 152.07256519, 121.50963612,
                  59.35127044,  14.68795807},{101.86262846, 153.60768909,  16.66114356,  59.35127044,
                    83.23303058,  96.31357746},{74.82005621, 206.62431451, -85.35062017,  14.68795807,
                      96.31357746, 159.6450431}};
      
          double[][] mat3 = new double[][]{{26.4562668 , 21.12679316,  4.56359879,  9.60911787,  0.92571084,
            6.22640711},{21.12679316, 48.91315282, -2.95143365,  7.96438229,  6.92532842,
              14.68125834},{4.56359879, -2.95143365,  7.6917846 ,  5.19027487, -0.23674644,
                -1.07176858},{9.60911787,  7.96438229,  5.19027487,  8.92838589,  1.00468498,
                  2.44237471},{0.92571084,  6.92532842, -0.23674644,  1.00468498,  4.29701006,
                    3.68536888},{6.22640711, 14.68125834, -1.07176858,  2.44237471,  3.68536888,
                      8.20089656}};
      
          double[][] mat4 = new double[][]{{46.09186853, 24.80902188, 25.42182706, 18.55376353, 13.72774273,
            13.50110126},{24.80902188, 33.6152518 , 11.69779626, 17.95922033, 16.42389905,
              16.99451084},{25.42182706, 11.69779626, 21.91674414, 10.66258398,  9.50560219,
                6.99389858},{18.55376353, 17.95922033, 10.66258398, 16.76534964,  9.32654077,
                  9.49392435},{13.72774273, 16.42389905,  9.50560219,  9.32654077, 12.73004386,
                    10.36549815},{13.50110126, 16.99451084,  6.99389858,  9.49392435, 10.36549815,
                      12.06153423}};
      
          Matrix P = new Matrix(mat,6,6);
          Matrix Q = new Matrix(mat1,6,6);
      
          ArrayList<Matrix> set = new ArrayList<Matrix>();
          set.add(new Matrix(mat,6,6));
          set.add(new Matrix(mat1,6,6));
          set.add(new Matrix(mat2,6,6));
          set.add(new Matrix(mat3,6,6));
          set.add(new Matrix(mat4,6,6));
          Matrix eucl_mean = MatrixAux.euclideanMean(set);
          Matrix riem_mean = MatrixAux.riemannMean(set);
          test.writeMatrixToPython(eucl_mean.getArray(),"test/euclMeanTest.txt");
          test.writeMatrixToPython(riem_mean.getArray(),"test/riemMeanTest.txt");
          System.out.println(MatrixAux.RiemannianDistance(P, Q));
          Matrix geod = MatrixAux.RiemannianGeodesic(P,Q,0.1);
          test.writeMatrixToPython(geod.getArray(),"test/geodTest.txt");
    }

    public static void testPartialFit(){
        double[][] mat = new double[][]{{15.40964107,   8.5875238 ,   9.1229583 ,  12.67034506,
            6.07197335,   7.09535784}, {8.5875238 ,  70.23148676, -12.75052337,  31.62069215,
              67.88629983,  70.08624275}, {9.1229583 , -12.75052337,  17.39204523,   5.17331657,
                -15.34425627, -15.50582496}, {12.67034506,  31.62069215,   5.17331657,  39.94475374,
                  40.49988185,  43.30090519}, {6.07197335,  67.88629983, -15.34425627,  40.49988185,
                    77.47320063,  79.56611252}, {7.09535784,  70.08624275, -15.50582496,  43.30090519,
                      79.56611252,  84.76104473}};
          
          double[][] mat2 = new double[][]{{5.66821762e+02, 4.55513333e+02, 4.45601375e+02, 3.32766550e+02,
            2.16859391e+02, 7.51203320e+01},{4.55513333e+02, 4.36102150e+02, 2.78248376e+02, 2.19497727e+02,
              1.67506595e+02, 1.00176972e+02},{4.45601375e+02, 2.78248376e+02, 4.88420961e+02, 3.45040016e+02,
                1.94516275e+02, 4.34459884e-01},{3.32766550e+02, 2.19497727e+02, 3.45040016e+02, 2.48745558e+02,
                  1.42971860e+02, 9.80899376e+00},{2.16859391e+02, 1.67506595e+02, 1.94516275e+02, 1.42971860e+02,
                    9.38859854e+01, 2.28484556e+01},{7.51203320e+01, 1.00176972e+02, 4.34459884e-01, 9.80899376e+00,
                      2.28484556e+01, 4.26424943e+01}};
      
          double[][] mat1 = new double[][]{{548.44055506, 469.91721369, 226.61698296, 152.8056371 ,
            101.86262846,  74.82005621},{469.91721369, 565.24965173,  47.31755651,  92.94547321,
              153.60768909, 206.62431451},{226.61698296,  47.31755651, 281.80538437, 152.07256519,
                16.66114356, -85.35062017},{152.8056371 ,  92.94547321, 152.07256519, 121.50963612,
                  59.35127044,  14.68795807},{101.86262846, 153.60768909,  16.66114356,  59.35127044,
                    83.23303058,  96.31357746},{74.82005621, 206.62431451, -85.35062017,  14.68795807,
                      96.31357746, 159.6450431}};
      
          double[][] mat3 = new double[][]{{26.4562668 , 21.12679316,  4.56359879,  9.60911787,  0.92571084,
            6.22640711},{21.12679316, 48.91315282, -2.95143365,  7.96438229,  6.92532842,
              14.68125834},{4.56359879, -2.95143365,  7.6917846 ,  5.19027487, -0.23674644,
                -1.07176858},{9.60911787,  7.96438229,  5.19027487,  8.92838589,  1.00468498,
                  2.44237471},{0.92571084,  6.92532842, -0.23674644,  1.00468498,  4.29701006,
                    3.68536888},{6.22640711, 14.68125834, -1.07176858,  2.44237471,  3.68536888,
                      8.20089656}};
      
          double[][] mat4 = new double[][]{{46.09186853, 24.80902188, 25.42182706, 18.55376353, 13.72774273,
            13.50110126},{24.80902188, 33.6152518 , 11.69779626, 17.95922033, 16.42389905,
              16.99451084},{25.42182706, 11.69779626, 21.91674414, 10.66258398,  9.50560219,
                6.99389858},{18.55376353, 17.95922033, 10.66258398, 16.76534964,  9.32654077,
                  9.49392435},{13.72774273, 16.42389905,  9.50560219,  9.32654077, 12.73004386,
                    10.36549815},{13.50110126, 16.99451084,  6.99389858,  9.49392435, 10.36549815,
                      12.06153423}};
          
          Potato pot = new Potato(2, 0.1, 1, 50, new int[]{0,1});
          ArrayList<Matrix> calibration = new ArrayList<Matrix>();
          calibration.add(new Matrix(mat,6,6));
          calibration.add(new Matrix(mat1,6,6));
          calibration.add(new Matrix(mat2,6,6));
          calibration.add(new Matrix(mat3,6,6));
          pot.mean_cov = MatrixAux.riemannMean(calibration);
          pot.last_cov = new Matrix(mat4,6,6);
          double [] dist = new double[calibration.size()];
          for(int j=0;j<dist.length;j++) dist[j] = Math.log(MatrixAux.RiemannianDistance(pot.mean_cov, calibration.get(j)));
          double[] res = MatrixAux.meanAndStd(dist);
          pot.mean = res[0];
          pot.std = res[1];
          System.out.print("BEFORE PartialFit mean: ");
          System.out.println(pot.mean);
          
          System.out.print("BEFORE PartialFit std: ");
          System.out.println(pot.std);
      
      
          pot.partialFit();
      
          System.out.print("PartialFit mean: ");
          System.out.println(pot.mean);
          
          System.out.print("PartialFit std: ");
          System.out.println(pot.std);
          
          test.writeMatrixToPython(pot.mean_cov.getArray(), "test/testPartialFit.txt");
    }

    public static void testPredictProb(){
        ArrayList<Matrix> cov0 = test.readCovLinearFromPython("../TrialsMNE/cov0_10.txt",2);
    ArrayList<Matrix> cov1 = test.readCovLinearFromPython("../TrialsMNE/cov1_10.txt",2);
    ArrayList<Matrix> cov2 = test.readCovLinearFromPython("../TrialsMNE/cov2_10.txt",6);
    ArrayList<Matrix> calib0 = new ArrayList<Matrix>();
    ArrayList<Matrix> calib1 = new ArrayList<Matrix>();
    ArrayList<Matrix> calib2 = new ArrayList<Matrix>();

    for(int i=0;i<4;i++){
      calib0.add(cov0.remove(0));
      calib1.add(cov1.remove(0));
      calib2.add(cov2.remove(0));
    }
    ArrayList<ArrayList<Matrix>> calibration = new ArrayList<ArrayList<Matrix>>();
    calibration.add(calib0);
    calibration.add(calib1);
    calibration.add(calib2);

    ArrayList<ArrayList<Matrix>> covs = new ArrayList<ArrayList<Matrix>>();
    covs.add(cov0);
    covs.add(cov1);
    covs.add(cov2);

    ArrayList<Potato> potatoes = new ArrayList<Potato>();

    potatoes.add(new Potato(2, 0.1, 1, 50, new int[]{0,1}));
    potatoes.add(new Potato(2, 0.1, 1, 50, new int[]{0,1}));
    potatoes.add(new Potato(2, 0.1, 1, 50, new int[]{0,1}));

    for(int i=0;i<potatoes.size();i++){
      potatoes.get(i).mean_cov = MatrixAux.riemannMean(calibration.get(i));
      //potatoes.get(i).mean_cov.print(2,2);
      potatoes.get(i).last_cov = calibration.get(i).get(calibration.get(i).size()-1);
      double [] dist = new double[calibration.get(i).size()];
      for(int j=0;j<dist.length;j++) dist[j] = Math.log(MatrixAux.RiemannianDistance(potatoes.get(i).mean_cov, calibration.get(i).get(j)));
      //for(int j=0;j<dist.length;j++) System.out.println(dist[j]);
      double[] res = MatrixAux.meanAndStd(dist);
      potatoes.get(i).mean = res[0];
      potatoes.get(i).std = res[1];
      //System.out.print(potatoes.get(i).mean);
      //System.out.print(" ");
      //System.out.println(potatoes.get(i).std);
    }

    for(int j=0;j<covs.get(0).size();j++){
      double probabilities = 0;
      for(int i=0;i<potatoes.size();i++){
        double prob = potatoes.get(i).predictProbability(covs.get(i).get(j));
        //System.out.print(prob);
        //System.out.print("; ");
        if(prob<1e-10) probabilities += Math.log(1e-10);
        else probabilities += Math.log(prob);
      }
      //System.out.println(" ");
      //System.out.println(probabilities);
      System.out.println(MatrixAux.getChiProb(-2*probabilities,3));
    }
}
}
