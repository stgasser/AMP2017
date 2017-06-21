package at.jku.cp.spezi.beta;

import at.jku.cp.spezi.dsp.AudioFile;
import at.jku.cp.spezi.dsp.Frame;
import at.jku.cp.spezi.dsp.Processor;
import org.math.plot.Plot2DPanel;

import javax.swing.*;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author andreas arzt
 * @author rainer kelz
 * @author franz strasser
 * @author stefan gasser
 */
public class Beta implements Processor {

    private static final boolean PLOT = false;
    private static final boolean MIX = false;

    /*******************************************************************************
     * 								Magic Numbers                                  *
     *******************************************************************************/
    //TODO Random Search for Parameters
    private static final int HOPSIZE = 512;
    private static final int FFTSIZE = 2 * HOPSIZE;           // For Reasons

    public static int MAXFILTER_TIME_WINDOWSIZE = 4;        //4
    public static int MAXFILTER_FREQ_WINDOWSIZE = 1;        //1

    public static double lambda = 52.610014112036495;                      //52.610014112036495
    public static double PPTS = 1.32846349236593;                        //1.32846349236593

    public static int w1 = 1, w2 = 2, w3 = 17, w4 = 7, w5 = 18;                   // 1 2 17 7 18

    public static int PEAK_PICKING_MEDIAN = 3;
    public static int PEAK_PICKING_LOCAL_MAX = 3;
    public static int PEAK_PICKING_AVG_MAX = 30;


    public static double PEAK_PICKING_THRESHOULD = 0.3;  //0.3
    public static double PEAK_PICKING_THRESHOULD_NEG_DECAY_LIMIT = 0.5; // 0.5
    public static int PEAK_PICKING_MEAN_WINDOWSIZE = 7;  //7
    //public static int PEAK_PICKING_MAX_WINDOWSIZE = 5;   //4
    public static double PEAK_PICKING_MAX_WINDOWSIZE_TS = 0.0;
    public static int PEAK_PICKING_LOCAL_MAX_WINDOWSIZE = 6;   //4
    public static int PEAK_PICKING_POST_ONSET_IGNORE = 7;     //7

    //private static final double PEAK_PICKING_THRESHOULD_REL = 1.3; //0.3
    //private static final double PEAK_PICKING_THRESHOULD_NEG_DECAY_LIMIT_REL = -0.7; // 0.5
    //private static final int PEAK_PICKING_POST_ONSET_IGNORE_REL = 11;//7

    // Onset Detection 2
    private static final double PP2_HFC_WEIGHT = 0.2795;
    private static final int PP2_MEDIAN_WINDOWSIZE = 11; // fft != hop: 9
    private static final double PP2_MEDIAN_SCALING_FACTOR = 1.5350923217009074; // fft != hop: 1.38070923, MIX: 1.0769
    private static final int PP2_LOCAL_WINDOWSIZE = 0; // fft != hop: 2
    private static final int PP2_POST_ONSET_IGNORE = 5; // fft != hop: 4
    private static final int PP2_FFTSIZE = 1024;
    private static final int PP2_HOPSIZE = 1024;

    // Beat detection
    private static final int BT_TEMPO_RANGE_LOW = 60;
    private static final int BT_TEMPO_RANGE_HIGH = 200;
    private static final int BT_LOCAL_WINDOWSIZE = 10;
    private static final int BT_BEAT_WINDOW_FACTOR = 8;    // hoptime * factor = distance to beat in both directions
    private static final double BT_OCTAVE_TOLERANCE = 0.40729824993564334;
    private static final double BT_FIRST_BEAT_WINDOW = 1.7; //in seconds
    private static final int BT_TEMPO_PULS_NR_FRAMES = 3;
    // Beat detection search
//    public static int BT_LOCAL_WINDOWSIZE;
//    public static int BT_BEAT_WINDOW_FACTOR;
//    public static double BT_OCTAVE_TOLERANCE;
//    public static double BT_FIRST_BEAT_WINDOW;
//    public static int BT_TEMPO_PULS_NR_FRAMES;

    public static double f3;
    public static int f1, f2, f4, f5;

    private String filename;
    private double sampleTime;

    /**
     * this list contains the results of the onset detection step
     * <p>
     * (time is in seconds)
     */
    private List<Double> onsets1;
    private List<Double> onsets2;
    private double[][] onsetDetectionFunction;

    /**
     * this list contain the results of the beat detection step
     * <p>
     * (beat times in seconds)
     */
    private List<Double> beats;

    private double tempo = 0.0;

    public Beta() {
    }

    public void process(String filename) {
        if (PLOT) {
            System.out.println("Initializing Processor '" + Beta.class.getName() + "'...");
        }

        onsets1 = new ArrayList<Double>();
        onsets2 = new ArrayList<Double>();
        beats = new ArrayList<Double>();
        tempo = new ArrayList<Double>();
        this.filename = filename;

        if (PLOT) {
            System.out.println("Running Analysis...");
        }

        onsetDetection1();
        if (PLOT) {
            System.out.println("Onset detection 1: " + onsets1.size());
        }

        onsetDetection2();
        if (PLOT) {
            System.out.println("Onset detection 2: " + onsets2.size());
        }

        beatDetection();
        tempoEstimation();
    }

    public List<Double> getOnsets() {
        return onsets1;
    }

    public List<Double> getBeats() {
        return beats;
    }

    public List<Double> getTempo() {
        List<Double> tempo = new ArrayList<>();
        //System.out.println(this.tempo);
        tempo.add(this.tempo);
        return tempo;
    }

    /**
     * SF => fmeasure = 89.59 %
     */
    private void onsetDetection1() {
        AudioFile audioFile = new AudioFile(filename, FFTSIZE, HOPSIZE);

        double sampleTime = 1d / audioFile.getSampleRate();
        this.sampleTime = sampleTime;

        List<Frame> frames = audioFile.getFrames();
        SpectralTransformator.sampleRate = audioFile.getSampleRate();
        SpectralTransformator.fftSize = FFTSIZE;
        SpectralTransformator.l = lambda;
        List<double[]> magdiff = new ArrayList<>(frames.size());
        List<double[]> mel = frames.stream()
                .map(frame -> {
                    double[] vals = new double[frame.magnitudes.length];
                    for (int i = 0; i < vals.length; i++) {
                        vals[i] = frame.magnitudes[i];//* Math.cos(frame.phases[i]);
                        //vals[1][i] = frame.magnitudes[i] * Math.sin(frame.phases[i]);
                    }
                    return vals;
                })
                .map(SpectralTransformator::toMel)
                .collect(Collectors.toList());
        double[] first = new double[mel.get(0).length];
        for (int i = 0; i < first.length; i++) {
            first[i] = mel.get(0)[i] / 4.0;
            //first[1][i] = mel.get(0)[1][i] / 4.0;
        }
        magdiff.add(first);
        for (int i = 0; i < mel.size() - 1; i++) {
            double[] tmp = new double[mel.get(i).length];
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = mel.get(i + 1)[j] - (mel.get(i)[j] - magdiff.get(i)[j]);
                //tmp[1][j] = mel.get(i + 1)[1][j] - (mel.get(i)[1][j] - magdiff.get(i)[1][j]);
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
                //spec[i][1][j] = maxim;
            }
        }

        List<double[]> magdiff2 = new ArrayList<>(spec.length);
        magdiff2.add(spec[0]);
        for (int i = 0; i < spec.length - 1; i++) {
            double[] tmp = spec[i];
            for (int j = 0; j < tmp.length; j++) {
                //tmp[j] = Math.sqrt(Math.abs(tmp[j] * tmp[j] - spec[i + 1][j] * spec[i + 1][j]));
                tmp[j] = H(spec[i + 1][j] - tmp[j]);
                //tmp[1][j] = H(spec[i + 1][1][j] - tmp[1][j]);
                /*if(tmp[j]>0.01)
                    tmp[j] = spec[i+1][j]/tmp[j];
                else
                    tmp[j] = spec[i+1][j];*/
            }
            magdiff2.add(tmp);
        }


        //magdiff2 = Arrays.asList(spec);
        List<Double> myList = magdiff2.stream()
                .map((doubles -> {
                    double sum = 0;
                    for (int i = 0; i < doubles.length; i++) {
                        sum += doubles[i];
                    }
                    return Math.abs(sum);
                })).sorted(Double::compareTo)
                .collect(Collectors.toList());
        double max = myList.get(myList.size() - 1);
        myList = magdiff2.stream()
                .map((doubles -> {
                    double sum = 0;
                    for (int i = 0; i < doubles.length; i++) {
                        sum += doubles[i];
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
            if (ismax && i - lastonset >= w5 && data[1][i] >= mean + PPTS) onsets1.add(data[0][i]);
        }
        // Maxfiltering
        /*for (int i = 0; i <= data[1].length - 5; i++) {
            double[] tmp = new double[5];
            for (int j = 0; j < 5; j++) {
                tmp[j] = data[1][i + j];
            }
            data[1][i] = Math.max(Math.max(tmp[0], tmp[1]), tmp[2]);
        }*/

        // save onset detection function for beat detection
        onsetDetectionFunction = data;

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
     * HFC + WPD => fmeasure = 77.16 %
     */
    private void onsetDetection2() {
        AudioFile audioFile = new AudioFile(filename, PP2_FFTSIZE, PP2_HOPSIZE);
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

        if (MIX) {
            double wHfc = PP2_HFC_WEIGHT;
            double[] all = new double[audioFile.getNrOfFrames()];
            for (int i = 0; i < audioFile.getNrOfFrames(); i++) {
                all[i] = wHfc * hfc[i] + (1 - wHfc) * wpd[i];
            }
            performPeakPickingFronz(all, audioFile);
        } else {
            performPeakPickingFronz(hfc, audioFile);
            performPeakPickingFronz(wpd, audioFile);

            /**
             * eliminate double onsets
             */
            onsets2 = onsets2.stream().sorted().collect(Collectors.toList());
            List<Double> tmpOnsets = new ArrayList<>();
            double sampleTime = 1d / audioFile.getSampleRate();
            double ignoreOnsetTime = sampleTime * PP2_POST_ONSET_IGNORE;
            for (int i = 0; i < onsets2.size(); i++) {
                double first = onsets2.get(i);

                while (i < onsets2.size() && Math.abs(onsets2.get(i) - first) <= 2 * ignoreOnsetTime) {
                    i++;
                }
                i--;

                double second = onsets2.get(i);
                double mean = (first + second) / 2;

                tmpOnsets.add(mean);
            }
            onsets2 = tmpOnsets;
        }
    }

    private void performPeakPickingFronz(double[] detection, AudioFile audioFile) {
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
                        onsets2.add(data[0][i]);
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
                            onsets2.add(data[0][i]);
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
     * Beat detection => fmeasure = 74.49 %
     */
    private void beatDetection() {
        /**
         * autocorrelation
         */
        // calculate the range of tau values dependeing on tempo range
        double hopTime = HOPSIZE * sampleTime;
        int lowerBound = -1;
        int upperBound = -1;
        for (int i = 1; i < onsetDetectionFunction[0].length; i++) {
            if (upperBound < 0 && (60.0 / (i * hopTime)) < BT_TEMPO_RANGE_LOW) {
                upperBound = i;
            }
            if (lowerBound < 0 && (60.0 / (i * hopTime)) <= BT_TEMPO_RANGE_HIGH) {
                lowerBound = i;
            }
        }
        int tauSize = upperBound - lowerBound;
        // calculate autocorrelation
        double[][] autocorrelation = new double[2][tauSize];
        for (int i = lowerBound; i < upperBound; i++) {
            int index = i - lowerBound;
            autocorrelation[0][index] = i * hopTime;

            double sum = 0;
            for (int t = 0; t < onsetDetectionFunction[1].length; t++) {
                if ((t+i) < 0 || (t+i) >= onsetDetectionFunction[1].length) {
                    continue;
                }

                sum += (onsetDetectionFunction[1][t+i] * onsetDetectionFunction[1][t]);
            }

            autocorrelation[1][index] = sum;
        }

        /**
         * perform peak picking
         */
        double maximum = 0;
        int maxIdx = -1;
        Map<Integer, Double> localMaxima = new HashMap<>();
        for (int i = 0; i < autocorrelation[0].length; i++) {
            double localMax = 0;
            int lmaxIdx = -1;
            for (int j = -BT_LOCAL_WINDOWSIZE/2; j <= BT_LOCAL_WINDOWSIZE/2; j++) {
                if ((i+j) < 0 || (i+j) >= autocorrelation[0].length) {
                    continue;
                }

                double value = autocorrelation[1][i+j];
                if (localMax < value) {
                    localMax = value;
                    lmaxIdx = i+j;
                }
            }
            if (i == lmaxIdx && localMax > maximum) {
                maximum = localMax;
                maxIdx = lmaxIdx;
            }
            if (i == lmaxIdx) {
                // local maximum found
                localMaxima.put(lmaxIdx, localMax);
            }
        }
        double tempoEstimation = autocorrelation[0][maxIdx];
        // pick the right local maximum
        int size = localMaxima.size();
        if (size <= 0) {
            // PROBLEM
            System.out.println("No beat estimation found!");
            return;
        }
        double minLocalMaximum = tempoEstimation;
        int minLocalMaxIdx = maxIdx;
        double mean = mean(localMaxima.values());
        for (int key : localMaxima.keySet()) {
            double value = autocorrelation[0][key];
            if (localMaxima.get(key) >= mean && value < minLocalMaximum && isOctaveOf(value, tempoEstimation)) {
                minLocalMaximum = value;
                minLocalMaxIdx = key;
            }
        }

        // calculate beats
        double beatTime = autocorrelation[0][minLocalMaxIdx];
        final double beatWindow = BT_BEAT_WINDOW_FACTOR * hopTime;
        List<Double> onsets = getOnsets();

        /**
         * find first beat location
         */
        int correlationSize = (int)(BT_FIRST_BEAT_WINDOW / hopTime);
        double pulseLength = BT_TEMPO_PULS_NR_FRAMES * hopTime;
        byte[] tempoTrainPuls = new byte[correlationSize];
        // create train pulse
        tempoTrainPuls[0] = 1;
        for (int i = 1; i < correlationSize; i++) {
            if ((i * hopTime) % minLocalMaximum <= pulseLength) {
                tempoTrainPuls[i] = 1;
            } else {
                tempoTrainPuls[i] = 0;
            }
        }
        // calculate crosscorrelation between train pulse and onsetdetectionfunction for first few seconds
        double[][] crosscorrelation = new double[2][correlationSize];
        for (int i = 0; i < correlationSize; i++) {
            crosscorrelation[0][i] = i * hopTime;

            double sum = 0;
            for (int j = 0; j < correlationSize; j++) {
                if ((j+i) >= onsetDetectionFunction[1].length) {
                    continue;
                }

                sum += (tempoTrainPuls[j] * onsetDetectionFunction[1][j+i]);
            }

            crosscorrelation[1][i] = sum;
        }
        if (PLOT) {
            Plot2DPanel panel = new Plot2DPanel();
            panel.addLinePlot("crosscorrelation", crosscorrelation[0], crosscorrelation[1]);
            JFrame frame = new JFrame(filename + ": crosscorrelation");
            frame.setContentPane(panel);
            frame.setSize(1000, 500);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        }
        // find maximum crosscorrelation
        double maxCc = 0;
        double maxCcOffset = 0;
        for (int i = 0; i < correlationSize; i++) {
            double cc = crosscorrelation[1][i];
            if (cc > maxCc) {
                maxCc = cc;
                maxCcOffset = crosscorrelation[0][i];
            }
        }
        // set first beat location
        double estimatedBeatTime = maxCcOffset;
        beats.add(estimatedBeatTime);

        /**
         * find all other beat locations
         */
        while (estimatedBeatTime <= onsets.get(onsets.size()-1)) {
            estimatedBeatTime += beatTime;

            // find all onsets in range of the current estimated beat
            List<Double> temp = new ArrayList<>();
            for (int i = 0; i < onsets.size(); i++) {
                double onset = onsets.get(i);
                if (onset < estimatedBeatTime - beatWindow) {
                    continue;
                }
                if (onset > estimatedBeatTime + beatWindow) {
                    break;
                }
                temp.add(onset);
            }

            double bestBeatTime = 0;
            if (temp.isEmpty()) {
                // no onset in range -> just take estimation
                bestBeatTime = estimatedBeatTime;
            } else {
                // at least 1 onset in range -> take nearest
                double bestTimeDiff = Double.MAX_VALUE;
                for (double onset : temp) {
                    double diff = Math.abs(estimatedBeatTime - onset);
                    if (diff < bestTimeDiff) {
                        bestTimeDiff = diff;
                        bestBeatTime = onset;
                    }
                }
                // update the estimation to the current onset
                estimatedBeatTime = bestBeatTime;
            }

            beats.add(bestBeatTime);
        }

        // @TODO = Tempo estimation!!!
        tempo.add(60 / tempoEstimation);

        if (PLOT) {
            Plot2DPanel panel = new Plot2DPanel();
            panel.addLinePlot("autocorrelation", autocorrelation[0], autocorrelation[1]);
            JFrame frame = new JFrame(filename + ": autocorrelation");
            frame.setContentPane(panel);
            frame.setSize(1000, 500);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        }
    }

    /**
     * Tempo estimation
     */
    private void tempoEstimation() {
        //System.out.println("Starting Tempo Estimation (NOT IMPLEMENTED!) ...");
        int FFTSIZE = 512, HOPSIZE = 512;
        AudioFile audioFile = new AudioFile(filename, FFTSIZE, HOPSIZE);

        double sampleTime = 1d / audioFile.getSampleRate();

        List<Frame> frames = audioFile.getFrames();
        SpectralTransformator.sampleRate = audioFile.getSampleRate();
        SpectralTransformator.fftSize = FFTSIZE;
        SpectralTransformator.l = lambda;
        List<double[]> magdiff = new ArrayList<>(frames.size());
        List<double[]> mel = frames.stream()
                .map(frame -> {
                    double[] vals = new double[frame.magnitudes.length];
                    for (int i = 0; i < vals.length; i++) {
                        vals[i] = frame.magnitudes[i];
                    }
                    return vals;
                })
                .map(SpectralTransformator::toMel)
                .collect(Collectors.toList());
        /*double[] first = new double[mel.get(0).length];
        for (int i = 0; i < first.length; i++) {
            first[i] = mel.get(0)[i] / 4.0;
        }
        magdiff.add(first);
        for (int i = 0; i < mel.size() - 1; i++) {
            double[] tmp = new double[mel.get(i).length];
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = mel.get(i + 1)[j] - (mel.get(i)[j] - magdiff.get(i)[j]);
            }
            magdiff.add(tmp);
        }*/
        magdiff = mel;
        // Max Filtering
        double[][] spec = new double[magdiff.size()][magdiff.get(0).length];
        for (int i = 0; i < magdiff.size(); i++) {                         //i = Current Time Slot
            for (int j = 0; j < magdiff.get(0).length; j++) {              //j = Current Freq Slot
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
                if (Double.isInfinite(max))
                    System.out.println();
                spec[i][j] = max;//magdiff.get(i)[j];
            }
        }

        List<double[]> magdiff2 = new ArrayList<>(spec.length);
        magdiff2.add(spec[0]);
        for (int i = 0; i < spec.length - 1; i++) {
            double[] tmp = new double[spec[0].length];
            for (int j = 0; j < tmp.length; j++) {
                tmp[j] = H(spec[i + 1][j] - spec[i][j]);
                //tmp[j] = Math.sqrt(Math.abs(tmp[j] * tmp[j] - spec[i + 1][j] * spec[i + 1][j]));
            }
            magdiff2.add(tmp);
        }


        List<Double> myList = magdiff2.stream()
                .map((doubles -> {
                    double sum = 0;
                    for (int i = 0; i < doubles.length; i++) {
                        sum += doubles[i];
                    }
                    return Math.abs(sum);
                })).collect(Collectors.toList());
        //-------------- Same as Onset Detection  till here--------------------


        double[] data = new double[myList.size()];
        for (int i = 0; i < data.length; i++) {
            data[i] = myList.get(i);
        }
        double[] datamed = new double[data.length];
        final int medwinsz = 20;
        final double medamp = 1.5;
        for (int i = 0; i < data.length; i++) {
            double[] tmp = Arrays.copyOfRange(data, i - medwinsz < 0 ? 0 : i - medwinsz, i + medwinsz > data.length - 1 ? data.length - 1 : i + medwinsz);
            Arrays.sort(tmp);
            datamed[i] = H(data[i] - medamp * (tmp.length % 2 == 0 ? (tmp[tmp.length / 2] + tmp[tmp.length / 2 - 1]) / 2 : tmp[tmp.length / 2]));
        }
        data = datamed;
        int maxbeat = Math.round(((60.0f / 35.0f) * audioFile.getSampleRate()) / HOPSIZE);
        int minbeat = Math.round((60.0f / 240.0f) * audioFile.getSampleRate() / HOPSIZE);
        double[] acf = new double[maxbeat - minbeat];
        double[] offsets = new double[acf.length];
        for (int i = minbeat; i < maxbeat; i++) {
            int j = i - minbeat;
            acf[j] = autcorrelate(data, i);
            offsets[j] = i * sampleTime * HOPSIZE;
        }

        acf = weightedacf(acf, minbeat, maxbeat);
        double[] nacf = new double[acf.length - 1];
        for (int i = 1; i < acf.length; i++) {
            //nacf[i - 1] = Math.abs(acf[i] - acf[i - 1]);
        }
        //acf = nacf;
        int maxidx = 0;
        for (int i = 0; i < acf.length; i++) {
            if (acf[maxidx] < acf[i]) {
                maxidx = i;
            }
        }
        int maxidxd2 = (maxidx + minbeat) / 2 - minbeat;
        int maxidxm2 = (maxidx + minbeat) * 2 - minbeat;
        int maxidxm3 = (maxidx + minbeat) * 3 - minbeat;
        int maxidxm4 = (maxidx + minbeat) * 4 - minbeat;
        System.out.println("--------------------" + filename);
        System.out.println(acf[maxidx]);
        if (maxidxd2 >= 0)
            System.out.println("/2 " + acf[maxidxd2]);
        if (maxidxm2 + minbeat < maxbeat)
            System.out.println("*2 " + acf[maxidxm2]);
        if (maxidxm3 + minbeat < maxbeat)
            System.out.println("*3 " + acf[maxidxm3]);
        if (maxidxm4 + minbeat < maxbeat)
            System.out.println("*4 " + acf[maxidxm4]);
        if (maxidxm2 + minbeat < maxbeat && maxidxm3 + minbeat < maxbeat){
            if (acf[maxidxm2] * 0.89 < acf[maxidxm3])
                maxidx = maxidxm3;
        }
        this.tempo = 1 / offsets[maxidx] * 60;
        if (PLOT) {
            Plot2DPanel panel = new Plot2DPanel();
            panel.addLinePlot("data", offsets, acf);
            //panel.addLinePlot("data", data);
            //panel.addLinePlot("data",data[0],newvals);
            JFrame frame = new JFrame(filename);
            frame.setContentPane(panel);
            frame.setSize(1000, 500);
            frame.setVisible(true);
            frame.setDefaultCloseOperation(WindowConstants.DISPOSE_ON_CLOSE);
        }

    }

    private double[] weightedacf(double[] acf, int minframesperbeat, int maxframesperbeat) {
        double[] ret = new double[acf.length];

        for (int i = 0; i < acf.length; i++) {
            ret[i] = 0;
            int currper = minframesperbeat + i;
            /*double period = 2 * Math.PI / (currper);
            for (int j = i; j < acf.length; j++) {
                int framecnt = j + minframesperbeat;
                ret[i] += acf[j] * Math.cos(period * framecnt) * gauss(framecnt, i + minframesperbeat, (maxframesperbeat - minframesperbeat) / 2);
            }
            //ret[i]/= acf.length-i;
            */
            int cnt = 1;
            for (int j = currper; j < acf.length + minframesperbeat; j += currper) {
                ret[i] += acf[j - minframesperbeat] / (cnt * cnt);
                cnt++;
            }
            //ret[i] /= (cnt);
        }
        return ret;
    }

    private double autcorrelate(double[] data, int offset) {
        double sum = 0;
        for (int i = 0; i < data.length - offset; i++) {
            sum += data[i + offset] * data[i];//Math.sqrt(H(data[i + offset]*data[i + offset] - data[i]*data[i]));
        }
        return sum;
    }

    private double gauss(double x, double mu, double sigma) {
        return gauss((x - mu) / sigma) / sigma;
    }

    private double gauss(double x) {
        return Math.exp(-(x * x) / 2) / Math.sqrt(2 * Math.PI);
    }

    private double H(double d) {
        return d < 0 ? 0 : d;
    }

    private boolean isOctaveOf(double v1, double v2) {
        double higher = Math.max(v1, v2);
        double lower = Math.min(v1, v2);

        do {
            lower *= 2;
            double diff = Math.abs(higher - lower);
            if (diff <= higher * BT_OCTAVE_TOLERANCE) {
                return true;
            }
        } while (lower < higher);

        return false;
    }

    private double mean(Collection<Double> m) {
        double sum = 0;
        for (double d : m) {
            sum += d;
        }
        return sum / m.size();
    }
}
