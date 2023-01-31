import Jama.Matrix;
import com.github.psambit9791.jdsp.filter.*;




public class SignalProcAux {
    public static double[][] filter_bandpass(double[][] signal, int l_f, int h_f, int[] channels){
        int Fs=500;
        int order = 8;
        Butterworth flt = new Butterworth(Fs);
        int m = channels.length;
        int n = signal[0].length;
        double[][] procSignal = new double[m][n];
        for(int ch=0;ch<m;ch++){
            procSignal[ch] = flt.bandPassFilter(signal[channels[ch]], order, l_f, h_f);
        }
        //return new Matrix(procSignal, m, n);*/
        return procSignal;
    }

    public static double[][] notch_filter(double[][] signal){
        return signal;  //TODO: implement the notch filter
    }

    

}
