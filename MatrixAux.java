import Jama.Matrix;
import Jama.EigenvalueDecomposition;
import java.util.ArrayList;
import smile.stat.distribution.*;
import org.apache.commons.math3.linear.MatrixUtils;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.*;
import org.apache.commons.math3.stat.descriptive.*;

//import dev.nm.*;
//import dev.nm.algebra.linear.matrix.doubles.matrixtype.dense.DenseMatrix;
//import dev.nm.stat.covariance.nlshrink.LedoitWolf2016;


public class MatrixAux {
    
    public static double sqr(double x){return x*x;}	

    // Non-integer power of a matrix
    public static Matrix power(Matrix M, double p)
    {
        EigenvalueDecomposition evd=M.eig();
        Matrix D=evd.getD();
        // power on the positive eigen values
        for(int i=0; i<D.getColumnDimension();i++)
            D.set(i,i,Math.pow(D.get(i,i),p));

        Matrix V=evd.getV();

        return V.times(D.times(V.transpose()));
    }
    // Exponential of a SPD matrix	
    public static Matrix exp(Matrix M)
    {
    EigenvalueDecomposition evd=M.eig();
    Matrix D=evd.getD();

    // exponential on the positive eigen values
    for(int i=0; i<D.getColumnDimension();i++)
        D.set(i,i,Math.exp(D.get(i,i)));

    Matrix V=evd.getV();

    return V.times(D.times(V.transpose()));
    }

    // Logarithm of a SPD matrix	
    public static Matrix log(Matrix M)
    {
    EigenvalueDecomposition evd=M.eig();
    Matrix D=evd.getD();

    // Logarithm on the positive eigen values
    for(int i=0; i<D.getColumnDimension();i++)
        D.set(i,i,Math.log(D.get(i,i)));

    Matrix V=evd.getV();

    return V.times(D.times(V.transpose()));
    }

    public static Matrix LogEuclideanMean(ArrayList<Matrix> set)
    {
        int i;
        int n=set.size(); 
        int d=set.get(0).getColumnDimension();
        Matrix R=new Matrix(d,d);

        for(i=0;i<n;i++){
            R=R.plus(log(set.get(i)).times(1.0/n))	;
        }

        return exp(R);
    }

    public static Matrix euclideanMean(ArrayList<Matrix> set)
    {
        int i;
        int n=set.size(); 
        int d=set.get(0).getColumnDimension();
        Matrix R=new Matrix(d,d);

        for(i=0;i<n;i++){
            R=R.plus(set.get(i).times(1.0/n));
        }

        return R;
    }

    public static Matrix riemannMean(ArrayList<Matrix> set)
    {
        int max_iter = 50;
        double tol = 10e-9;
        double nu = 1.0;
        double tau = Double.MAX_VALUE;
        double crit = Double.MAX_VALUE;

        Matrix C = euclideanMean(set);
        int d=set.get(0).getColumnDimension();

        for(int i=0;i<max_iter;i++)
        {
            Matrix C12 = power(C, 0.5);
            Matrix Cm12 = power(C.inverse(),0.5);

            Matrix J= new Matrix(d, d);
            for(Matrix M : set){
                //Cm12.times(M).times(Cm12).print(6,6);
                J = J.plus(log(Cm12.times(M).times(Cm12))); 
            }
            C = C12.times(exp(J.times(nu))).times(C12);
            //C.print(6,6);

            crit = J.normF();
            double h = nu*crit;
            if(h<tau){
                nu = 0.95*nu;
                tau=h;
            }
            else {
                nu=0.5*nu;
            }

            if(crit<=tol || nu<=tol){
                break;
            }

        }
        return C;
    }

    public static double RiemannianDistance(Matrix P, Matrix Q)
    {
    double result=0.0d;
    Matrix M=P.inverse().times(Q);
    EigenvalueDecomposition evd=M.eig();
    Matrix D=evd.getD();

    for(int i=0; i<D.getColumnDimension();i++)
        result+=sqr(Math.log(D.get(i,i)));

    return Math.sqrt(result);	
    }

    public static Matrix RiemannianGeodesic(Matrix P, Matrix Q, double lambda)
    {
    Matrix result;
    Matrix Phalf=power(P,0.5);
    Matrix Phalfinv=power(P,-0.5);

    result=Phalfinv.times(Q).times(Phalfinv);
    result=power(result,lambda);

    return (Phalf.times(result)).times(Phalf);	
    }

    public static double transform(Matrix P, Matrix M){
    return RiemannianDistance(P, M);
    }

    public static double[] meanAndStd(double[] sd) {
    
    double sum = 0;
    double newSum = 0; 

    for (int i = 0; i < sd.length; i++) {
        sum = sum + sd[i];
    }
    double mean = (sum) / (sd.length); 

    for (int j = 0; j < sd.length; j++) {
        // put the calculation right in there
        newSum = newSum + ((sd[j] - mean) * (sd[j] - mean)); 
    }
    double squaredDiffMean = (newSum) / (sd.length); 
    double standardDev = (Math.sqrt(squaredDiffMean)); 
    return new double[]{mean, standardDev};

    //System.out.println(mean);
    //System.out.println(standardDev);

        //DescriptiveStatistics stat = new DescriptiveStatistics(sd);
        //return new double[]{stat.getMean(),stat.getStandardDeviation()};
    }

    public static double getZScore(double dist, double mean, double std){
    return (dist-mean)/std;
    }

    // return pdf(x) = standard Gaussian pdf
    public static double pdf(double x) {
    return Math.exp(-x*x / 2) / Math.sqrt(2 * Math.PI);
    }

    // return cdf(z) = standard Gaussian cdf using Taylor approximation
    public static double normCDF(double z) {
    if (z < -8.0) return 0.0;
    if (z >  8.0) return 1.0;
    double sum = 0.0, term = z;
    for (int i = 3; sum + term != sum; i += 2) {
        sum  = sum + term;
        term = term * z * z / i;
    }
    return 0.5 + sum * pdf(z);
    }

    public static double getNormProb(double z){
        
    //return 1-normCDF(z);
      GaussianDistribution g = new GaussianDistribution(0,1);
      double dist = 1- g.cdf(z);
      return dist;
    }

    public static double getChiProb(double z, int n_potatoes){
        ChiSquareDistribution c = new ChiSquareDistribution(2*n_potatoes);
        return 1 - c.cdf(z);
    }

    public static Matrix estimate_centroid(ArrayList<Matrix> X){
        return LogEuclideanMean(X);
    }
    public static double returnMeanOfArray(double[] a){
        double sum = 0;
        for(int i=0;i<a.length;i++){
            sum+=a[i];
        }
        return sum/a.length;
        
    }

    public static Matrix regularizeCovarianceMatrix(double[][] sig){
        return new Matrix(3,3);
    }
    public static Matrix computeCovarianceMatrixRegularized(double[][] sig){ //TO BE IMPLEMENTED
        //convert channels into indexes: TO BE IMPLEMENTED
        int n_samples = sig[0].length;
        int n_features = sig.length;
        
        for(int i=0;i<n_features;i++){
            double av = returnMeanOfArray(sig[i]);
            for(int j=0;j<n_samples;j++){
                sig[i][j]-=av;
            }
        }
        Matrix aux = new Matrix(sig, n_features,n_samples);
        Matrix emp_cov = aux.times(aux.transpose()).times(1.0/n_samples);
        double mu = emp_cov.trace()/1.0/n_features;
        double alpha = returnMeanOfArray(emp_cov.arrayTimes(emp_cov).getRowPackedCopy());
        double num = alpha + mu*mu;
        double den = (n_samples + 1.0) * (alpha - (mu*mu) / n_features);
        double shrinkage = 1.0;
        if(den!=0)
            shrinkage = Math.min(num / den, 1.0);
        //System.out.println("mu: "+mu+"\nalpha: "+alpha+"\nnum: "+num+"\nden: "+den+"\nshrinkage: "+shrinkage);
        //emp_cov.times(1e12).print(10,10);
        Matrix shrunk_cov = emp_cov.times(1.0 - shrinkage);
        //shrunk_cov.times(1e12).print(10,10);
        //System.out.println(shrinkage * mu);
        for(int i=0;i<n_features;i++){
            shrunk_cov.set(i, i, shrunk_cov.get(i,i)+(shrinkage * mu));
        }
        //RealMatrix cov = new Covariance(MatrixUtils.createRealMatrix(sig).transpose().getData()).getCovarianceMatrix();
        
        //Matrix covar = new Matrix(cov.getData());
        return shrunk_cov;
    }

    public static Matrix computeCovarianceMatrix(double[][] sig){ //TO BE IMPLEMENTED
        //convert channels into indexes: TO BE IMPLEMENTED
        
        RealMatrix cov = new Covariance(MatrixUtils.createRealMatrix(sig).transpose().getData()).getCovarianceMatrix();
        Matrix covar = new Matrix(cov.getData());
        return covar;
    }

    public static Matrix normalizeCovMatrix(Matrix M, String method){
        int d=M.getColumnDimension();
        double[] stddev = new double[d];
        Matrix normM = new Matrix(d, d);
        double aux=0;
        if(method=="cov"){
            for(int i=0;i<d;i++){
                stddev[i]= Math.sqrt( Math.abs(M.get(i, i) ));
            }
            for(int i=0;i<d;i++){
                for(int j=0;j<d;j++){
                    aux = M.get(i,j)/(stddev[i]*stddev[j]);
                    if(Math.abs(aux)>1) normM.set(i,j,Math.signum(aux));
                    else normM.set(i,j,aux);
                }
            }
        }
        else{
            double trace = M.trace();
            for(int i=0;i<d;i++){
                for(int j=0;j<d;j++){
                    normM.set(i,j,M.get(i,j)/trace);
                }
            }
        }
        return normM;
    }


}

