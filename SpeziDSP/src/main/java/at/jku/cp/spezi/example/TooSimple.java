package at.jku.cp.spezi.example;

import java.util.ArrayList;
import java.util.List;

import at.jku.cp.spezi.dsp.AudioFile;
import at.jku.cp.spezi.dsp.Processor;

/**
 *
 * @author andreas arzt
 * @author rainer kelz
 */
public class TooSimple implements Processor {

	private AudioFile audioFile;

	/**
	 * this list contains the results of the onset detection step
	 * 
	 * (time is in seconds)
	 */
	private List<Double> onsets;

	/**
	 * this list contain the results of the beat detection step
	 * 
	 * (beat times in seconds)
	 */
	private List<Double> beats;

	public TooSimple() {
	}

	public void process(String filename) {
		System.out.println("Initializing Processor '" + TooSimple.class.getName() + "'...");
		onsets = new ArrayList<Double>();
		beats = new ArrayList<Double>();
		
		// an AudioFile object is created with the following parameters
		// AudioFile(WAVFILENAME, fftSize, hopSize, OPTIONAL: window function)
		// sizes are in samples

		// the WAV files you were provided with are all sampled at 44100 Hz

		// if you would like to work with multiple DFT resolutions, you would
		// simply create multiple AudioFile objects with different parameters
		System.out.println("Computing STFT ...");
		this.audioFile = new AudioFile(filename, 2048, 1024);

		System.out.println("Running Analysis...");
		onsetDetection();
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
	 * 
	 * TODO: you have to implement *at least* 2 more different onset detection
	 * functions have a look at the class 'Frame' - it contains the magnitudes,
	 * the phases, and more which you can use to implement your detection
	 * function
	 */
	private void onsetDetection() {
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

	/**
	 * TODO: we do not provide any beat detection example implementation. you
	 * ned to implement *at least* two different beat detection functions.
	 */
	private void beatDetection() {
		System.out.println("Starting Beat Detection (NOT IMPLEMENTED!) ...");
	}

	/**
	 * TODO: we do not provide any beat detection example implementation. you
	 * ned to implement *at least* two different tempo estimation functions.
	 */
	private void tempoEstimation() {
		System.out.println("Starting Tempo Estimation (NOT IMPLEMENTED!) ...");
	}
}
