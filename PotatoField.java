import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

import Jama.Matrix;

public class PotatoField {
    double z_threshold;
    double p_threshold;
    double alpha;
    ArrayList<Potato> potatoes = new ArrayList<Potato>();
    int n_potatoes;
    static int n_iter_max = 10;
    static int pos_label = 1;
    static int neg_label = 0;
    //Map<String, Integer> elMap = new HashMap<String, Integer>();
    static int[] eeg_ch = {2,4,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21,22,23,24,25,26,27,28,29,30,31}; 
    
    

    public PotatoField(double z, double p, double alpha, int[] min_fs, int[] max_fs, int[][] channels)
    {/* 
      elMap.put("GROUND_ELECTRODE",0);
      elMap.put("HEART_R1",1);
      elMap.put("PZ_ELECTRODE",2);
      elMap.put("HEART_TIP",3);
      elMap.put("FZ_ELECTRODE",4);
      elMap.put("HEART_R2",5);
      elMap.put("CZ_ELECTRODE",6);
      elMap.put("FP2_ELECTRODE",7);
      elMap.put("T9_ELECTRODE",8);
      elMap.put("C4_ELECTRODE",15);
      elMap.put("C4_ELECTRODE",14);
      elMap.put("C4_ELECTRODE",11);
      elMap.put("C4_ELECTRODE",12);
      elMap.put("C4_ELECTRODE",13);
      elMap.put("C4_ELECTRODE",10);
      elMap.put("C4_ELECTRODE",9);
      elMap.put("C4_ELECTRODE",20);
      elMap.put("C4_ELECTRODE",17);
      elMap.put("C4_ELECTRODE",18);
      elMap.put("C4_ELECTRODE",23);
      elMap.put("C4_ELECTRODE",16);
      elMap.put("C4_ELECTRODE",21);
      elMap.put("C4_ELECTRODE",22);
      elMap.put("C4_ELECTRODE",19);
      elMap.put("C4_ELECTRODE",28);
      elMap.put("C4_ELECTRODE",25);
      elMap.put("C4_ELECTRODE",26);
      elMap.put("C4_ELECTRODE",31);
      elMap.put("C4_ELECTRODE",24);
      elMap.put("C4_ELECTRODE",29);
      elMap.put("C4_ELECTRODE",30);
      elMap.put("C4_ELECTRODE",27);
*/
      z_threshold=z;
      p_threshold=p;
      for(int n=0;n<min_fs.length;n++)
      {
        if(channels[n].length != 0){
          potatoes.add(new Potato(z, alpha, min_fs[n], max_fs[n], channels[n]));
        }
        else{
          potatoes.add(new Potato(z, alpha, min_fs[n], max_fs[n], eeg_ch));
        }
        
      }
      n_potatoes = min_fs.length;
    }


    public void fitField(double[][][] X){
      for (Potato p : potatoes) {
        ArrayList<Matrix> fittable = new ArrayList<Matrix>();
        for (double[][] m : X) {
          fittable.add(MatrixAux.computeCovarianceMatrix(m));
        }
        p.fit(fittable);
      }
    }

    public void partialFitField(){
      for (Potato p : potatoes) {
        p.partialFit();
      }
    }

    public double predictProbabilityField(double[][] M){
      double probabilities = 0;
      for (Potato p : potatoes) {
        Matrix m = MatrixAux.computeCovarianceMatrix(M);
        double prob = p.predictProbability(m);
        if(prob<1e-10) probabilities -= 23.0; //log(1e-10)
        else probabilities += Math.log(prob);
      }
      return MatrixAux.getChiProb(-2*probabilities,n_potatoes);
    }

    public boolean predictField(double[][] M){
      if(predictProbabilityField(M)>p_threshold) return true; //clean
      return false;
    }


}
