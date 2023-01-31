import Jama.Matrix;
import java.util.ArrayList;
import java.util.Arrays;

import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.Covariance;

//class potato
public class Potato{
    double threshold;
    double mean;
    double std;

    Matrix mean_cov;
    Matrix last_cov;

    double min_freq;
    double max_freq;
    double alpha;
    String[] channels;
    int[] ch_indexes;
    static int n_iter_max=100;
    static int pos_label=1;
    static int neg_label=0;
  
  
    public Potato(double th, double al, double minF, double maxF, int[] chan)
    {
      threshold=th;
      min_freq=minF;
      max_freq=maxF;
      ch_indexes=chan;
      alpha=al;

    }
    

    public void fit(ArrayList<Matrix> X){
      int[] y_old = new int[X.size()];
      Arrays.fill(y_old, 1);
      for(int i=0;i<n_iter_max;i++){
  
        ArrayList<Matrix> X_n = new ArrayList<Matrix>();
  
        for(int j=0;j<y_old.length;j++){
          if(y_old[j]==1) X_n.add(X.get(j));
        }
        mean_cov = MatrixAux.estimate_centroid(X_n);
        double [] dist = new double[X_n.size()];
        for(int j=0;j<dist.length;j++) dist[j] = Math.log(MatrixAux.RiemannianDistance(mean_cov, X_n.get(j)));
  
        double[] res = MatrixAux.meanAndStd(dist);
        mean = res[0];
        std = res[1];
        int[] y_n = new int[y_old.length];
        for(int j=0;j<dist.length;j++){
          if(MatrixAux.getZScore(dist[j],mean, std) < threshold) y_n[j]=1;
        }
  
        if(Arrays.equals(y_old,y_n)) break;
        else y_old = y_n.clone();
  
      }
    }
  
    public void partialFit(){
      mean_cov = MatrixAux.RiemannianGeodesic(mean_cov, last_cov, alpha); //rememeber to check alpha or 1/alpha
      double d = Math.log(MatrixAux.transform(mean_cov, last_cov));
      mean = (1 - alpha) * mean + alpha * d;
      std = Math.sqrt((1 - alpha) * std*std + alpha * (d - mean)*(d - mean));
    }
  
    public boolean predict(Matrix M){ //true if chunk is clean
      double z = MatrixAux.transform(mean_cov, M);
      if(z<threshold) return true;
      return false;
    }
    
    public double predictProbability(Matrix M){ //It is considered as abnormal/artifacted for low value of proba.
      double d = Math.log(MatrixAux.transform(mean_cov, M));
      double z = MatrixAux.getZScore(d, mean, std);
      return MatrixAux.getNormProb(z); 
    }
  
  }
  