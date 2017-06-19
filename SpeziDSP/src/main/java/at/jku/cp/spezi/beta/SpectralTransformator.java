package at.jku.cp.spezi.beta;

/**
 * Created by Stefan on 31.05.2017.
 */
public class SpectralTransformator {


    static int numBands = 50;
    static float sampleRate = 44100;
    static int fftSize = 2048;
    static double l = 1.0;

    public static double[] toMel(double[] fftData) {
        double[] bandData = new double[numBands];
        int MinFreq = 0;
        int MaxFreq = (fftData.length - 1) * Math.round(sampleRate / fftSize);
        int[] centerFreqs = new int[numBands + 2];
        double melMaxFreq = 2959 * Math.log(1 + MaxFreq / 700.0);
        double melMinFreq = 2959 * Math.log(1 + MinFreq / 700.0);

        centerFreqs[0] = MinFreq;
        centerFreqs[numBands + 1] = MaxFreq;

        for (int i = 0; i < numBands; i++) {
            centerFreqs[i + 1] = Math.round((float) ((Math.exp((melMaxFreq - melMinFreq) / numBands * i / 2959) - 1) * 700)*fftSize/sampleRate);
        }

        for (int bandIdx = 1; bandIdx < numBands; bandIdx++) {
            int startFreqIdx = centerFreqs[bandIdx - 1];
            int centerFreqIdx = centerFreqs[bandIdx];
            int stopFreqIdx = centerFreqs[bandIdx + 1];
            int magnitudeScale;
            double height = 2.0 / (stopFreqIdx - startFreqIdx); //Normalize area under Triangles
            for (int i = startFreqIdx; i < centerFreqIdx; i++) {
                magnitudeScale = centerFreqIdx - startFreqIdx;
                bandData[bandIdx] += fftData[i] * height * (i - startFreqIdx) / magnitudeScale;
            }

            for (int i = centerFreqIdx; i <= stopFreqIdx; i++) {
                magnitudeScale = centerFreqIdx - stopFreqIdx;
                bandData[bandIdx] += fftData[i] * height * (i - stopFreqIdx) / magnitudeScale;
            }
        }
        //double l = 1.0;
        for (int i = 0; i < bandData.length; i++) {
            bandData[i] = Math.log(1 + l * bandData[i]);
        }

        return bandData;
    }
}
