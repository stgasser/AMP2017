package at.jku.cp.spezi.beta;

import at.jku.cp.spezi.dsp.AudioFile;
import at.jku.cp.spezi.dsp.Frame;
import at.jku.cp.spezi.dsp.Processor;
import org.math.plot.Plot2DPanel;

import javax.swing.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * @author andreas arzt
 * @author rainer kelz
 * @author franz strasser
 * @author stefan gasser
 */
public class Beta implements Processor {

    private static final boolean PLOT = true;

    /*******************************************************************************
     * 								Magic Numbers                                  *
     *******************************************************************************/
    //TODO Random Search for Parameters
    private static final int HOPSIZE = 512;
    private static final int FFTSIZE = 2 * HOPSIZE;           // For Reasons

    public static int MAXFILTER_TIME_WINDOWSIZE = 3;    //3
    public static int MAXFILTER_FREQ_WINDOWSIZE = 1;   //2

    public static double lambda = 1.0;   //49
    public static double PPTS = 1.0;    //0.87

    public static int w1, w2, w3, w4, w5;// 1 5 14 7 16

    public static int PEAK_PICKING_MEDIAN = 3;
    public static int PEAK_PICKING_LOCAL_MAX = 3;
    public static int PEAK_PICKING_AVG_MAX = 30;


    /*public static double PEAK_PICKING_THRESHOULD = 0.3;  //0.3
    public static double PEAK_PICKING_THRESHOULD_NEG_DECAY_LIMIT = 0.5; // 0.5
    public static int PEAK_PICKING_MEAN_WINDOWSIZE = 7;  //7
    public static int PEAK_PICKING_MAX_WINDOWSIZE = 5;   //4
    public static double PEAK_PICKING_MAX_WINDOWSIZE_TS = 0.0;
    public static int PEAK_PICKING_LOCAL_MAX_WINDOWSIZE = 6;   //4
    public static int PEAK_PICKING_POST_ONSET_IGNORE = 7;     //7

    private static final double PEAK_PICKING_THRESHOULD_REL = 1.3; //0.3
    private static final double PEAK_PICKING_THRESHOULD_NEG_DECAY_LIMIT_REL = -0.7; // 0.5
    private static final int PEAK_PICKING_POST_ONSET_IGNORE_REL = 11;//7
       */
    private static final int PP2_MEDIAN_WINDOWSIZE = 11;
    private static final double PP2_MEDIAN_SCALING_FACTOR = 1.12;
    private static final int PP2_LOCAL_WINDOWSIZE = 4;
    private static final int PP2_POST_ONSET_IGNORE = 4;

    private String filename;

    private AudioFile audioFile;

    /**
     * this list contains the results of the onset detection step
     * <p>
     * (time is in seconds)
     */
    private List<Double> onsets;

    /**
     * this list contain the results of the beat detection step
     * <p>
     * (beat times in seconds)
     */
    private List<Double> beats;

    public Beta() {
    }

    public void process(String filename) {
        //System.out.println("Initializing Processor '" + Beta.class.getName() + "'...");
        onsets = new ArrayList<Double>();
        beats = new ArrayList<Double>();
        this.filename = filename;
        // an AudioFile object is created with the following parameters
        // AudioFile(WAVFILENAME, fftSize, hopSize, OPTIONAL: window function)
        // sizes are in samples

        // the WAV files you were provided with are all sampled at 44100 Hz

        // if you would like to work with multiple DFT resolutions, you would
        // simply create multiple AudioFile objects with different parameters
        //System.out.println("Computing STFT ...");
        this.audioFile = new AudioFile(filename, 2048, 1024);

        //System.out.println("Running Analysis...");
        //onsetDetection();
        onsetDetection1();
        onsetDetection2();
        System.out.println(onsets.size());
        beatDetection();
        tempoEstimation();
    }

    public List<Double> getOnsets() {
        return onsets;
    }

    public List<Double> getBeats() {
        return beats;
    }

    public List<Double> getTempo() {
        List<Double> tempo = new ArrayList<>();
        tempo.add(150d);
        return tempo;
    }

    /**
     * this is a 'Signal Envelope' implementation
     * <p>
     * TODO: you have to implement *at least* 2 more different onset detection
     * functions have a look at the class 'Frame' - it contains the magnitudes,
     * the phases, and more which you can use to implement your detection
     * function
     */
    private void onsetDetection() {
        onsetDetection1();
        System.out.println("Starting Onset Detection ...");

        // this is the time difference between two consecutive samples
        double sampleTimeInSeconds = 1d / audioFile.getSampleRate();

        // this list stores the audio samples from the WAV file
        List<Double> samples = audioFile.getSamples();

        // if we wanted the list of STFT frames, we'd call
        // List<Frame> frames = audioFile.getFrames();

        // average samples in a window extending 400 samples into the past and
        // the future
        // 'x' is the current sample
        // |..............x..............|
        int w = 400;

        // where did this threshold come from? for the purpose of this example,
        // we pulled it out of our hat. for the competition you should
        // definitely
        // try to tune any such 'magic numbers' ...
        double threshold = 0.35;

        // - run over all samples from the signal
        // - compute the average over the energy in a window
        // - report everything larger than that average plus a threshold
        for (int i = w; i < samples.size() - w; i++) {
            double mean = 0d;
            for (int j = -w; j < w; j++) {
                mean = mean + Math.abs(samples.get(i + j));
            }
            mean = mean / (2 * w + 1);

            // if the current sample-value is greater than
            //
            // threshold + mean(window)
            //
            // we report an onset ...

            if (threshold + mean < samples.get(i)) {
                onsets.add(i * sampleTimeInSeconds);
            }
        }
    }

    private void onsetDetection1() {

        AudioFile audioFile = new AudioFile(filename, FFTSIZE, HOPSIZE);

        double sampleTime = 1d / audioFile.getSampleRate();

        List<Frame> frames = audioFile.getFrames();
        SpectralTransformator.sampleRate = audioFile.getSampleRate();
        SpectralTransformator.fftSize = FFTSIZE;
        SpectralTransformator.l = lambda;
        List<double[]> magdiff = new ArrayList<>(frames.size());
        List<double[]> mel = frames.stream()
                .map(frame -> frame.magnitudes)
                .map(SpectralTransformator::toMel)
                .collect(Collectors.toList());
        double[] first = new double[mel.get(0).length];
        for (int i = 0; i < first.length; i++) {
            first[i] = mel.get(0)[i] / 4.0;
        }
        magdiff.add(first);
        for (int i = 0; i < mel.size() - 1; i++) {
            double[] tmp = new double[mel.get(i).length];
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = mel.get(i + 1)[j] - (mel.get(i)[j] - magdiff.get(i)[j]);
                //tmp[j] = mel.get(i)[j];// - (mel.get(i)[j] - magdiff.get(i)[j]);
            }
            magdiff.add(tmp);
        }

        // Max Filtering
        double[][] spec = new double[magdiff.size()][magdiff.get(0).length];
        for (int i = 0; i < magdiff.size(); i++) {                         //i = Current Time Slot
            for (int j = 0; j < magdiff.get(i).length; j++) {              //j = Current Freq Slot
                double max = Double.NEGATIVE_INFINITY;
                for (int k = 0; k < MAXFILTER_TIME_WINDOWSIZE; k++) {
                    for (int l = 0; l < MAXFILTER_FREQ_WINDOWSIZE; l++) {
                        int currTimeIdx = i + k - (MAXFILTER_TIME_WINDOWSIZE) / 2;
                        int currFreqIdx = j + l - (MAXFILTER_FREQ_WINDOWSIZE) / 2;
                        if (currTimeIdx >= 0 && currTimeIdx < magdiff.size() && currFreqIdx >= 0 && currFreqIdx < magdiff.get(currTimeIdx).length) {
                            double currVal = magdiff.get(currTimeIdx)[currFreqIdx];
                            if (currVal > max) max = currVal;
                        }
                    }
                }
                spec[i][j] = max;
            }
        }

        List<double[]> magdiff2 = new ArrayList<>(spec.length);
        magdiff2.add(spec[0]);
        for (int i = 0; i < spec.length - 1; i++) {
            double[] tmp = spec[i];
            for (int j = 0; j < tmp.length; j++) {
                //tmp[j] = Math.sqrt(Math.abs(tmp[j] * tmp[j] - spec[i + 1][j] * spec[i + 1][j]));
                //tmp[j] = H(spec[i + 1][j] - tmp[j]);
                if(tmp[j]>0.01)
                    tmp[j] = spec[i+1][j]/tmp[j];
                else
                    tmp[j] = spec[i+1][j];
            }
            magdiff2.add(tmp);
        }


        //magdiff2 = Arrays.asList(spec);
        List<Double> myList = magdiff2.stream()
                .map((doubles -> {
                    double sum = 0;
                    for (double d :
                            doubles) {
                        sum += d;
                    }
                    return Math.abs(sum);
                })).sorted(Double::compareTo)
                .collect(Collectors.toList());
        double max = myList.get(myList.size() - 1);
        myList = magdiff2.stream()
                .map((doubles -> {
                    double sum = 0;
                    for (double d :
                            doubles) {
                        sum += d;
                    }
                    return Math.abs(sum);
                })).collect(Collectors.toList());
        double[][] data = new double[2][magdiff2.size()];
        for (int i = 0; i < magdiff2.size(); i++) {
            data[0][i] = i * sampleTime * HOPSIZE;
            data[1][i] = (myList.get(i));/// max;
        }
        /*
        int lastOnset = 0;
        for (int i = 0; i < data[0].length; i++) {
            double mean = 0, windowmax = Double.NEGATIVE_INFINITY;
            int datacnt = 0;
            for (int j = -PEAK_PICKING_MEAN_WINDOWSIZE / 2; j < PEAK_PICKING_MEAN_WINDOWSIZE / 2; j++) {
                int currindex = i + j;
                if (currindex < 0) continue;
                if (currindex >= data[0].length) continue;
                mean += data[1][currindex];
                datacnt++;
                if (windowmax < data[1][currindex] && j >= -PEAK_PICKING_MAX_WINDOWSIZE / 2 && j < PEAK_PICKING_MAX_WINDOWSIZE / 2) {
                    windowmax = data[1][currindex];
                }
            }
            mean = mean / datacnt;
            if (data[1][i] >= windowmax*PEAK_PICKING_MAX_WINDOWSIZE_TS) {
                if (/*data[1][i] / mean > PEAK_PICKING_THRESHOULD_REL * Math.max((lastOnset + PEAK_PICKING_POST_ONSET_IGNORE_REL - i) / PEAK_PICKING_POST_ONSET_IGNORE_REL, -PEAK_PICKING_THRESHOULD_NEG_DECAY_LIMIT_REL) ||*/
                        /*data[1][i] > mean + PEAK_PICKING_THRESHOULD * Math.max((lastOnset + PEAK_PICKING_POST_ONSET_IGNORE - i) / PEAK_PICKING_POST_ONSET_IGNORE, -PEAK_PICKING_THRESHOULD_NEG_DECAY_LIMIT)) {
                    boolean isMaxima = true;
                    for (int j = 0; j < PEAK_PICKING_LOCAL_MAX_WINDOWSIZE; j++) {
                        int idx = j + i - PEAK_PICKING_LOCAL_MAX_WINDOWSIZE / 2;
                        if (idx >= 0 && idx < data[1].length)
                            if (data[1][idx] > data[1][i]) isMaxima = false;
                    }
                    if (isMaxima){
                        onsets.add(data[0][i]);
                        lastOnset = i;
                    }
                    //i += PEAK_PICKING_POST_ONSET_IGNORE;
                }
            }
        }*/
        // Peak picking
        /*double[] newvals = new double[data[1].length];
        for (int i = 0; i < data[1].length; i++) {
            double[] arr = Arrays.copyOfRange(data[1], Math.max(i - PEAK_PICKING_MEDIAN, 0), Math.min(i + PEAK_PICKING_MEDIAN, data[1].length - 1));
            Arrays.sort(arr);
            newvals[i] = data[1][i] - arr[(Math.min(i + 3, data[1].length - 1) - Math.max(i - 3, 0)) / 2];
            if (newvals[i] < 0) newvals[i] = 0;
        }
        data[1] = newvals;
        newvals = new double[data[1].length];
        for (int i = 0; i < data[1].length; i++) {
            double[] arr = Arrays.copyOfRange(data[1], Math.max(i - PEAK_PICKING_AVG_MAX, 0), Math.min(i + PEAK_PICKING_AVG_MAX, data[1].length - 1));
            Arrays.sort(arr);
            if (arr[arr.length - 1] == 0) newvals[i] = 0;
            else newvals[i] = data[1][i] / arr[arr.length - 1];
            //System.out.println(newvals[i]);
        }

        for (int i = 0; i < newvals.length; i++) {
            boolean ismax = true;
            for (int j = 1; i - j >= 0 && j < PEAK_PICKING_LOCAL_MAX && ismax; j++) {
                if (newvals[i] < newvals[i - j]) ismax = false;
            }
            for (int j = 1; i + j < newvals.length && j < PEAK_PICKING_LOCAL_MAX && ismax; j++) {
                if (newvals[i] < newvals[i + j]) ismax = false;
            }
            if (ismax && newvals[i] > PPTS) onsets.add(data[0][i]);
        }*/
        //Peak Picking v3
        int lastonset = -w5;
        for (int i = 0; i < data[1].length; i++) {
            boolean ismax = true;
            for (int j = 1; j < w1 && i - j >= 0; j++) {
                if (data[1][i - j] > data[1][i]) ismax = false;
            }
            for (int j = 1; j < w2 && i + j < data[1].length; j++) {
                if (data[1][i + j] > data[1][i]) ismax = false;
            }
            int datacnt = 0;
            double sum = 0;
            for (int j = 0; j < w3 && i - j >= 0; j++) {
                datacnt++;
                sum += data[1][i - j];
            }
            for (int j = 1; j < w4 && i + j < data[1].length; j++) {
                datacnt++;
                sum += data[1][i + j];
            }
            double mean = sum / datacnt;
            if (ismax && i - lastonset >= w5 && data[1][i] >= mean + PPTS) onsets.add(data[0][i]);
        }
        // Maxfiltering
        /*for (int i = 0; i <= data[1].length - 5; i++) {
            double[] tmp = new double[5];
            for (int j = 0; j < 5; j++) {
                tmp[j] = data[1][i + j];
            }
            data[1][i] = Math.max(Math.max(tmp[0], tmp[1]), tmp[2]);
        }*/
        if (PLOT) {
            Plot2DPanel panel = new Plot2DPanel();
            panel.addLinePlot("data", data[0], data[1]);
            //panel.addLinePlot("data",data[0],newvals);
            JFrame frame = new JFrame(filename);
            frame.setContentPane(panel);
            frame.setSize(1000, 500);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        }
    }

    /**
     * HFC - High Frequency Content
     */
    private void onsetDetection2() {
        List<Frame> frames = audioFile.getFrames();

        /**
         * calculate hfc values
         */
        double hfcMax = 0;
        double[] hfc = new double[audioFile.getNrOfFrames()];
        for (int i = 0; i < audioFile.getNrOfFrames(); i++) {
            Frame frame = frames.get(i);
            double d = 0;
            for (int j = 0; j < frame.size; j++) {
                d += (j + 1) * Math.pow(frame.magnitudes[j], 2);
            }
            hfc[i] = d / frame.size;

            if (hfc[i] > hfcMax) {
                hfcMax = hfc[i];
            }
        }

        /**
         * calculate wpd values
         */
        double wpdMax = 0;
        double[] wpd = new double[audioFile.getNrOfFrames()];
        for (int i = 0; i < audioFile.getNrOfFrames(); i++) {
            Frame frame0 = frames.get(Math.min(i, audioFile.getNrOfFrames()));
            Frame frame1 = frames.get(Math.max(i, 0));
            Frame frame2 = frames.get(Math.max(i - 1, 0));
            double d = 0;
            for (int j = 0; j < frame1.size; j++) {
                d += Math.abs(frame1.magnitudes[j] * (frame0.magnitudes[j] - 2 * frame1.magnitudes[j] + frame2.magnitudes[j]));
            }
            wpd[i] = d / frame1.size;

            if (wpd[i] > wpdMax) {
                wpdMax = wpd[i];
            }
        }

        performPeakPickingFronz(hfc);
        performPeakPickingFronz(wpd);

        /**
         * eliminate double onsets
         */
        onsets = onsets.stream().sorted().collect(Collectors.toList());
        List<Double> tmpOnsets = new ArrayList<>();
        double sampleTime = 1d / audioFile.getSampleRate();
        double ignoreOnsetTime = sampleTime * PEAK_PICKING_POST_ONSET_IGNORE;
        for (int i = 0; i < onsets.size(); i++) {
            double first = onsets.get(i);

            while (i < onsets.size() && Math.abs(onsets.get(i) - first) <= 2 * ignoreOnsetTime) {
                i++;
            }
            i--;

            double second = onsets.get(i);
            double mean = (first + second) / 2;

            tmpOnsets.add(mean);
        }
        onsets = tmpOnsets;
    }

    private void performPeakPickingFronz(double[] detection) {
        double max = Arrays.stream(detection).max().getAsDouble();

        /**
         * map values to time values
         */
        double sampleTime = 1d / audioFile.getSampleRate();
        double[][] data = new double[2][detection.length];
        for (int i = 0; i < detection.length; i++) {
            data[0][i] = i * sampleTime * 1024;
            data[1][i] = detection[i] / max;
        }

        /**
         * perform peak picking
         */
        int lastOnset = -PP2_POST_ONSET_IGNORE;
        for (int i = 0; i < data[0].length; i++) {
            if (true) {
                if ((i - lastOnset) < PP2_POST_ONSET_IGNORE) {
                    i = lastOnset + PP2_POST_ONSET_IGNORE - 1;
                    continue;
                }

                double[] values = new double[PP2_MEDIAN_WINDOWSIZE];
                for (int j = 0; j < PP2_MEDIAN_WINDOWSIZE; j++) {
                    int index = i + j - (PP2_MEDIAN_WINDOWSIZE - 1) / 2;
                    if (index < 0 || index >= data[0].length) {
                        break;
                    }

                    values[j] = data[1][index];
                }

                values = Arrays.stream(values).sorted().toArray();
                double median;
                if (PP2_MEDIAN_WINDOWSIZE % 2 == 0) {
                    median = (values[PP2_MEDIAN_WINDOWSIZE / 2 - 1] + values[PP2_MEDIAN_WINDOWSIZE / 2]) / 2;
                } else {
                    median = values[(PP2_MEDIAN_WINDOWSIZE - 1) / 2];
                }
                median *= PP2_MEDIAN_SCALING_FACTOR;

                if (data[1][i] > median) {
                    boolean isMaxima = true;
                    for (int j = 0; j < PP2_LOCAL_WINDOWSIZE; j++) {
                        int idx = j + i - PP2_LOCAL_WINDOWSIZE / 2;
                        if (idx >= 0 && idx < data[1].length) {
                            if (data[1][idx] > data[1][i]) {
                                isMaxima = false;
                            }
                        }
                    }
                    if (isMaxima) {
                        onsets.add(data[0][i]);
                        lastOnset = i;
                    }
                }
            } else {
                double mean = 0;
                double windowmax = Double.NEGATIVE_INFINITY;
                int datacnt = 0;
                for (int j = -PEAK_PICKING_MEAN_WINDOWSIZE / 2; j < PEAK_PICKING_MEAN_WINDOWSIZE / 2; j++) {
                    int currindex = i + j;
                    if (currindex < 0) {
                        j = -i - 1;
                        continue;
                    }
                    if (currindex >= data[0].length) {
                        break;
                    }

                    mean += data[1][currindex];
                    datacnt++;
                    if (windowmax < data[1][currindex]) {
                        windowmax = data[1][currindex];
                    }
                }
                mean = mean / datacnt;
                if (data[1][i] >= windowmax * PEAK_PICKING_MAX_WINDOWSIZE_TS) {
                    if (data[1][i] > (mean + PEAK_PICKING_THRESHOULD * Math.max((lastOnset + PEAK_PICKING_POST_ONSET_IGNORE - i) / PEAK_PICKING_POST_ONSET_IGNORE, -PEAK_PICKING_THRESHOULD_NEG_DECAY_LIMIT))) {
                        boolean isMaxima = true;
                        for (int j = 0; j < PEAK_PICKING_LOCAL_MAX_WINDOWSIZE; j++) {
                            int idx = j + i - PEAK_PICKING_LOCAL_MAX_WINDOWSIZE / 2;
                            if (idx >= 0 && idx < data[1].length)
                                if (data[1][idx] > data[1][i]) isMaxima = false;
                        }
                        if (isMaxima) {
                            onsets.add(data[0][i]);
                            lastOnset = i;
                        }
                    }
                }
            }
        }

        /**
         * plot values
         */
        if (PLOT) {
            Plot2DPanel panel = new Plot2DPanel();
            panel.addLinePlot("data", data[0], data[1]);

            JFrame frame = new JFrame(filename);
            frame.setContentPane(panel);
            frame.setSize(1000, 500);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        }
    }

    /**
     * TODO: we do not provide any beat detection example implementation. you
     * ned to implement *at least* two different beat detection functions.
     */
    private void beatDetection() {
        //System.out.println("Starting Beat Detection (NOT IMPLEMENTED!) ...");
    }

    /**
     * TODO: we do not provide any beat detection example implementation. you
     * ned to implement *at least* two different tempo estimation functions.
     */
    private void tempoEstimation() {
        //System.out.println("Starting Tempo Estimation (NOT IMPLEMENTED!) ...");
    }

    private double H(double d) {
        return (d + Math.abs(d)) / 2.0;
    }
}
