import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.*;

import Jama.Matrix;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Scanner;  

import smile.stat.distribution.*;

public class Main {
  public static void main(String[] args) {

    //double[][] mat = test.readMatrixFromPython("../TrialsMNE/baseTest.txt",33,15360);
    //test.testCov(mat);

    //test.testDistance_Geod_Mean();

    //test.testPartialFit();

    //test.testPredictProb();

    //double[][] mat = test.readMatrixFromPython("../TrialsMNE/testCov.txt",10,500);
    //System.out.println(mat[0][0]);
    //Matrix reg = MatrixAux.computeCovarianceMatrixRegularized(mat);
    //reg.times(1e12).print(10,10);
    //mat1.times(1e12).print(10,10);


    //double[][] aux = new double[][]{{1,2,3,4,5},{1,2,3,4,5},{3,3,3,3,3}};
    //MatrixAux.computeCovarianceMatrix(aux).print(3,3);




    

    //RealMatrix cov = new Covariance(MatrixUtils.createRealMatrix(mat).transpose().getData()).getCovarianceMatrix();

    //System.out.println(cov.getColumnDimension());
    //System.out.println(cov.getRowDimension());
    
    // int[][][] oneSecChunks = new int[100][27][256];
    // int[][][] oneMinCalibration = new int[60][27][256];

    int[] min_fs = {1, 25, 1};
    int[] max_fs = {20, 45, 3};
    int[][] channels_idx = {
       {18, 7}, //FP1, FP2
       {26, 31}, //O1, O2
       {}
    };
    PotatoField rpf = new PotatoField(17.56, 0.00001, 0.01, min_fs, max_fs, channels_idx);
    double[][][] train = new double[][][]{
      {
        {1,2,-1,3},
        {5,-2,3,4},
        {1,-0.5,2,-0.3},
      },
      {
        {0.2,0.5,-0.2,-1.7},
        {2,3,-0.6,0.3},
        {0.1,-0.4,0.5,4},
      }
    };
    rpf.fitField(train);

    double[][][] test = new double[][][]{
      {
        {1.5,2.3,-1.2,2},
        {3,-1,0.3,0.4},
        {1.2,-1.5,0.2,0.3},
      },
      {
        {0.4,1.5,-1.2,-0.7},
        {2.5,1.3,-1.6,1.3},
        {0.3,-1.4,1.5,2.4},
      },
      {
        {10.4,-100.5,91.2,-100.7},
        {1.5,1.2,-2.6,-1.3},
        {22.3,-20.4,50.5,-200.4},
      }
    };
    boolean[] res=new boolean[test.length];
    for(int i=0;i<test.length;i++){
      res[i]=rpf.predictField(test[i]);
      //if(res[i]==true){
      //  rpf.partialFitField();
      //}
    }
    System.out.println(res[0]);
    System.out.println(res[1]);
    System.out.println(res[2]);
    // rpf.fit_potatoField(oneMinCalibration);
    // ArrayList<Boolean> predictions = new ArrayList<Boolean>();
    // for(int[][] oneSecChunk: oneSecChunks){
    //   boolean predic = rpf.predictPotatoField(oneSecChunk);
    //   if(predic) rpf.partialFit_potatoField();
    //   predictions.add(predic);
    // }
    // System.out.println("Hello World");

  }
}
