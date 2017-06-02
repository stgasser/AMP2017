/*
 * AudioFile.java
 * Representation of an Audiofile.
 *
 * After Initialization you have to call the process() method to do the feature extraction
 *
 * The most important methods for you are
 *   - List<Frame> getFrames() (returns a list of frames (the STFT))
 *   - List<Double> getSamples() (returns a list of samples (the raw audio))
 *
 * The signal is multiplied by a Hamming window by default.
 *
 */
package at.jku.cp.spezi.dsp;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;

/**
 *
 * @author andreas arzt
 */
public class AudioFile {

	private AudioInputStream rawInputStream;
	private AudioFormat audioFormat;
	private int channels;
	private float sampleRate;
	private AudioInputStream pcmInputStream;
	private byte[] inputBuffer;
	private double[] circBuffer;
	private double[] reBuffer;
	private double[] imBuffer;
	private int cbIndex;
	private double[] window;
	private double frameRMS;

	/**
	 * contains the spectral data (magnitude, phase, unwrapped phase) for each
	 * frame
	 */
	private List<Frame> frames;

	/**
	 * contains the value of each sample
	 */
	private List<Double> samples;

	private int hopSize;
	private int fftSize;

	public AudioFile(String fileName, int fftSize, int hopSize) {
		this(fileName, fftSize, hopSize, FFT.HAMMING);
	}

	public AudioFile(String fileName, int fftSize, int hopSize, int windowFunction) {
		try {
			System.out.println("Reading and preprocessing audiofile " + fileName);
			if (!isPowerOfTwo(fftSize)) {
				throw new IllegalArgumentException("fftSize must be a power of two");
			}

			if (!isPowerOfTwo(hopSize)) {
				throw new IllegalArgumentException("hopSize must be a power of two");
			}

			File audioFile = new File(fileName);
			if (!audioFile.isFile()) {
				throw new FileNotFoundException("Requested file does not exist: " + fileName);
			}

			this.fftSize = fftSize;
			this.hopSize = hopSize;

			rawInputStream = AudioSystem.getAudioInputStream(audioFile);
			audioFormat = rawInputStream.getFormat();
			channels = audioFormat.getChannels();
			sampleRate = audioFormat.getSampleRate();
			pcmInputStream = rawInputStream;

			if ((audioFormat.getEncoding() != AudioFormat.Encoding.PCM_SIGNED)
					|| (audioFormat.getFrameSize() != channels * 2)
					|| audioFormat.isBigEndian()) {

				AudioFormat desiredFormat = new AudioFormat(
					AudioFormat.Encoding.PCM_SIGNED,
					sampleRate,
					16,
					channels,
					channels * 2,
					sampleRate,
					false);

				pcmInputStream = AudioSystem.getAudioInputStream(
					desiredFormat,
					rawInputStream);
				audioFormat = desiredFormat;
			}

			int buffSize = hopSize * channels * 2;
			if ((inputBuffer == null) || (inputBuffer.length != buffSize)) {
				inputBuffer = new byte[buffSize];
			}

			if ((circBuffer == null) || (circBuffer.length != fftSize)) {
				circBuffer = new double[fftSize];
				reBuffer = new double[fftSize];
				imBuffer = new double[fftSize];
				window = FFT.makeWindow(windowFunction, fftSize, fftSize);
				for (int i = 0; i < fftSize; i++) {
					window[i] *= Math.sqrt(fftSize);
				}
			}

			cbIndex = 0;
			frameRMS = 0;

			frames = new ArrayList<Frame>();
			samples = new ArrayList<Double>();

			processFile();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	public List<Frame> getFrames() {
		return frames;
	}

	public List<Double> getSamples() {
		return samples;
	}

	public int getFftSize() {
		return fftSize;
	}

	public int getHopSize() {
		return hopSize;
	}

	public float getSampleRate() {
		return sampleRate;
	}
	
	public int getNrOfSamples() {
		return samples.size();
	}
	
	public int getNrOfFrames() {
		return frames.size();
	}

	/**
	 * checks to see whether the input is a power of two x = 2^i
	 * 
	 * @param x
	 * @return true if x = 2^i, false if it is not
	 */
	private boolean isPowerOfTwo(int x) {
		return (x != 0 && (x & (x - 1)) == 0);
	}

	/**
	 * Reads a frame of input data, averages the channels to mono, scales to a
	 * maximum possible absolute value of 1, and stores the audio data in a
	 * circular input buffer.
	 * 
	 * @return true if a frame (or part of a frame, if it is the final frame) is
	 *         read. If a complete frame cannot be read, the InputStream is set
	 *         to null.
	 */
	private boolean getFrame() {
		if (pcmInputStream == null) {
			return false;
		}
		try {
			int bytesRead = (int) pcmInputStream.read(inputBuffer);

			if (bytesRead < inputBuffer.length) {
				return false;
			}
		} catch (IOException e) {
			return false;
		}
		frameRMS = 0;
		double sample;
		switch (channels) {
		case 1:
			for (int i = 0; i < inputBuffer.length; i += 2) {
				sample = ((inputBuffer[i + 1] << 8)
						| (inputBuffer[i] & 0xff)) / 32768.0;
				frameRMS += sample * sample;
				circBuffer[cbIndex++] = sample;
				samples.add(sample);
				if (cbIndex == fftSize) {
					cbIndex = 0;
				}
			}
			break;
		case 2: // saves ~0.1% of RT (total input overhead ~0.4%) :)
			for (int i = 0; i < inputBuffer.length; i += 4) {
				sample = (((inputBuffer[i + 1] << 8) | (inputBuffer[i] & 0xff))
						+ ((inputBuffer[i + 3] << 8) | (inputBuffer[i + 2] & 0xff))) / 65536.0;
				frameRMS += sample * sample;
				circBuffer[cbIndex++] = sample;
				samples.add(sample);
				if (cbIndex == fftSize) {
					cbIndex = 0;
				}
			}
			break;
		default:
			for (int i = 0; i < inputBuffer.length;) {
				sample = 0;
				for (int j = 0; j < channels; j++, i += 2) {
					sample += (inputBuffer[i + 1] << 8) | (inputBuffer[i] & 0xff);
				}
				sample /= 32768.0 * channels;
				frameRMS += sample * sample;
				circBuffer[cbIndex++] = sample;
				samples.add(sample);
				if (cbIndex == fftSize) {
					cbIndex = 0;
				}
			}
		}
		frameRMS = Math.sqrt(frameRMS / inputBuffer.length);
		return true;
	}

	/**
	 * processes the audiofile; reads frames, computes the STFT
	 */
	private void processFile() {

		while (getFrame()) {
			for (int i = 0; i < fftSize; i++) {
				reBuffer[i] = window[i] * circBuffer[cbIndex];
				// reBuffer[i] = circBuffer[cbIndex];
				if (++cbIndex == fftSize) {
					cbIndex = 0;
				}
			}
			Arrays.fill(imBuffer, 0);
			FFT.magnitudePhaseFFT(reBuffer, imBuffer);
			Frame s = new Frame(reBuffer, imBuffer, fftSize);
			frames.add(s);
		}

		for (int i = 0; i < frames.size(); i++) {
			if (i == 0) {
				Arrays.fill(frames.get(i).unwrappedPhases, 0);
				continue;
			}
			frames.get(i).computeUnwrappedPhases(frames.get(i - 1).unwrappedPhases);
		}
	}

}
