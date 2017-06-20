/*
 * Runner.java
 * - contains the main method
 * - handles the input parameters
 * - writes the results to files
 * - evaluates the results
 * - summarizes the evaluation
 *
 * DO NOT CHANGE ANYTHING in this file!
 *
 */
package at.jku.cp.spezi;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.sql.Time;
import java.util.*;

import at.jku.cp.spezi.beta.Beta;
import at.jku.cp.spezi.dsp.Processor;
import joptsimple.OptionParser;
import joptsimple.OptionSet;

/**
 * @author andreas arzt
 * @author rainer kelz
 */
public class Runner {
	private static final boolean MONTECARLO = false;

	public static String ONSETS = ".onsets";
    public static String BEATS = ".beats";
    public static String TEMPO = ".tempo";

    private static String[] types = {ONSETS, BEATS, TEMPO};
    private static double[] tolerances = {0.05d, 0.07d, 4d};

    private static String AUDIO = ".wav";

    private static String GT = ".gt";
    private static String EV = ".ev";
    private static String PR = ".pr";

    private static Processor processor = null;

    private static void usage() {
        List<String> msg = Arrays.asList(
                "Usage  : <jarfile> -i <inputdirectory> [-n <processorname>] [-p|e|s]",
                "Options: ",
                "   -i input directory (the directory you want to analyze)",
                "   -n name of the processor, defaults to 'example.TooSimple'",
                "   -p predict onsets, beats and tempo for all .wav in the input directory",
                "   -e evaluate all predictions for which a groundtruth exists",
                "   -s summarize all evaluations");

        for (String line : msg) {
            System.out.println(line);
        }
    }

    public static void main(String[] args) throws Exception {
        OptionParser parser = new OptionParser("i:n:pes");
        OptionSet options = parser.parse(args);

        if (!options.has("i")) {
            System.out.println("ERROR: Input directory required! (-i INPUT)");
            usage();
            System.exit(1);
        }

        if (!(options.has("p") || options.has("e") || options.has("s"))) {
            System.out.println("ERROR: None of the options [p|e|s] given!");
            usage();
            System.exit(1);
        }

        String processorName = null;
        if (options.has("n")) {
            processorName = options.valueOf("n").toString();
        } else {
            System.out.println("Defaulting to Processor 'example.TooSimple'");
            processorName = "example.TooSimple";
        }

		if (MONTECARLO) {
			performMontecarloSearchFronz(options, processorName);
		} else {
			performStandardRun(options, processorName);
		}
	}

	private static void performStandardRun(OptionSet options, String processorName) throws IOException {
		String directory = options.valueOf("i").toString();
		processor = findAndInstantiateClass(processorName);

		if (options.has("p")) {
			System.out.println("Predicting ...");
			Files.walk(Paths.get(directory))
					.filter(path -> Files.isRegularFile(path))
					.filter(path -> path.toString().endsWith(AUDIO))
					//.limit(1)	// TODO hack for only first file for testing
					.forEach(path -> predictForFile(path.toString()));
		}

		if (options.has("e")) {
			System.out.println("Evaluating ...");
			Files.walk(Paths.get(directory))
					.filter(path -> Files.isRegularFile(path))
					.filter(path -> path.toString().endsWith(AUDIO))
					.forEach(path -> evalForFile(path.toString()));
		}

		if (options.has("s")) {
			System.out.println("Summarizing ...");
			for (String type : types) {
				Map<String, Integer> summary = new HashMap<>();
				summary.put("tp", 0);
				summary.put("fp", 0);
				summary.put("fn", 0);
				summary.put("error", 0);

				Files.walk(Paths.get(directory))
						.filter(path -> Files.isRegularFile(path))
						.filter(path -> path.toString().endsWith(type + EV))
						.forEach(path -> updateSummaryForFile(summary, path.toString()));

				double precision = 0d;
				double recall = 0d;
				double fmeasure = 0d;

				int tp = summary.get("tp");
				int fp = summary.get("fp");
				int fn = summary.get("fn");
				double error = (double) summary.get("error") / 1e6d;

				if (tp + fp > 0)
					precision = (double) tp / (tp + fp);

				if (tp + fn > 0)
					recall = (double) tp / (tp + fn);

				if (precision + recall > 0)
					fmeasure = (2 * precision * recall) / (precision + recall);

				String summaryfilename = directory + "/summary" + type + EV + ".txt";
				List<String> lines = Arrays.asList(
						"tp " + tp,
						"fp " + fp,
						"fn " + fn,
						"precision " + precision,
						"recall " + recall,
						"fmeasure " + fmeasure,
						"error " + error);
				Utils.writeToFile(summaryfilename, lines);
			}
		}
	}

	private static void performMontecarloSearchStefan(OptionSet options, String processorName) throws IOException {
		int tries = 0;
		int b_mf_tws = 4, b_mffws = 1, b_median = 0, b_avgmax = 0, b_max = 0;
		int b_w1 = 4, b_w2 = 4, b_w3 = 19, b_w4 = 11, b_w5 = 15;
		double b_ppts = 1.097764218995879, b_la = 47;
		double fmmax = 0;
		Random rand = new Random(System.currentTimeMillis());
		while (tries < 1000) {
			tries++;
			if (tries % 10 == 0)
				System.out.println("#" + tries);
			int mf_tws = b_mf_tws + rand.nextInt(3) - 1,
					mffws = b_mffws + rand.nextInt(3) - 1,
					w1 = b_w1 + rand.nextInt(5) - 2,
					w2 = b_w2 + rand.nextInt(5) - 2,
					w3 = b_w3 + rand.nextInt(5) - 2,
					w4 = b_w4 + rand.nextInt(5) - 2,
					w5 = b_w5 + rand.nextInt(5) - 2;
                    /*localmax = rand.nextInt(10) + 1,
                    avgmax = rand.nextInt(10) + 1,
                    median = rand.nextInt(30) + 1;*/

			double la = b_la + rand.nextDouble() * 11 - 5,
					ppts = b_ppts + rand.nextDouble() * b_ppts - b_ppts / 2;
			if (mf_tws < 1) mf_tws = 1;
			if (mffws < 1) mffws = 1;
			if (w1 < 1) w1 = 1;
			if (w2 < 1) w2 = 1;
			if (w3 < 1) w3 = 1;
			if (w4 < 1) w4 = 1;
			if (w5 < 1) w5 = 1;

			Beta.MAXFILTER_FREQ_WINDOWSIZE = mffws;
			Beta.MAXFILTER_TIME_WINDOWSIZE = mf_tws;
			Beta.lambda = la;
			Beta.w1 = w1;
			Beta.w2 = w2;
			Beta.w3 = w3;
			Beta.w4 = w4;
			Beta.w5 = w5;
			//Beta.PEAK_PICKING_MEDIAN = median;
			Beta.PPTS = ppts;
			//Beta.PEAK_PICKING_AVG_MAX = avgmax;
			//Beta.PEAK_PICKING_LOCAL_MAX = localmax;

			String directory = options.valueOf("i").toString();
			processor = findAndInstantiateClass(processorName);

			if (options.has("p")) {
				//System.out.println("Predicting ...");
				Files.walk(Paths.get(directory))
						.filter(path -> Files.isRegularFile(path))
						.filter(path -> path.toString().endsWith(AUDIO))
						//.limit(1)	// TODO hack for only first file for testing
						.forEach(path -> predictForFile(path.toString()));
			}

			if (options.has("e")) {
				//System.out.println("Evaluating ...");
				Files.walk(Paths.get(directory))
						.filter(path -> Files.isRegularFile(path))
						.filter(path -> path.toString().endsWith(AUDIO))
						.forEach(path -> evalForFile(path.toString()));
			}

			if (options.has("s")) {
				//System.out.println("Summarizing ...");
				for (String type : types) {
					Map<String, Integer> summary = new HashMap<>();
					summary.put("tp", 0);
					summary.put("fp", 0);
					summary.put("fn", 0);
					summary.put("error", 0);

					Files.walk(Paths.get(directory))
							.filter(path -> Files.isRegularFile(path))
							.filter(path -> path.toString().endsWith(type + EV))
							.forEach(path -> updateSummaryForFile(summary, path.toString()));

					double precision = 0d;
					double recall = 0d;
					double fmeasure = 0d;

					int tp = summary.get("tp");
					int fp = summary.get("fp");
					int fn = summary.get("fn");
					double error = (double) summary.get("error") / 1e6d;

					if (tp + fp > 0)
						precision = (double) tp / (tp + fp);

					if (tp + fn > 0)
						recall = (double) tp / (tp + fn);

					if (precision + recall > 0)
						fmeasure = (2 * precision * recall) / (precision + recall);
					if (type.equalsIgnoreCase(ONSETS) && fmeasure > fmmax) {
						fmmax = fmeasure;
						b_mf_tws = mf_tws;
						b_mffws = mffws;
						b_la = la;
						//b_avgmax = avgmax;
						//b_max = localmax;
						b_ppts = ppts;
						b_w1 = w1;
						b_w2 = w2;
						b_w3 = w3;
						b_w4 = w4;
						b_w5 = w5;
						//b_median = median;
						System.out.println("new best FMeasure: " + fmeasure);
						System.out.println(b_mf_tws);
						System.out.println(b_mffws);
						System.out.println(b_la);
						System.out.println(b_ppts);
						System.out.println(w1);
						System.out.println(w2);
						System.out.println(w3);
						System.out.println(w4);
						System.out.println(w5);
                        /*System.out.println(b_max);
                        System.out.println(b_avgmax);
                        System.out.println(b_median);*/
					}
					//System.out.println(fmeasure);
					String summaryfilename = directory + "/summary" + type + EV + ".txt";
					List<String> lines = Arrays.asList(
							"tp " + tp,
							"fp " + fp,
							"fn " + fn,
							"precision " + precision,
							"recall " + recall,
							"fmeasure " + fmeasure,
							"error " + error);
					//lines.forEach(System.out::println);
					//Utils.writeToFile(summaryfilename, lines);
				}
			}
		}
	}

	private static void performMontecarloSearchFronz(OptionSet options, String processorName) throws IOException {
		int tries = 0;
		Random rand = new Random(System.currentTimeMillis());
		String directory = options.valueOf("i").toString();
		List<String> lines = new ArrayList<>();
		lines.add("fmeasure;precision;recall;error;tp;fp;fn;"/*TempoRangeLow;TempoRangeHigh;*/+"LocalWindowSize;BeatWindowFactor;OctaveTolerance");

		long startTimeMillis = System.currentTimeMillis();
		Time now = new Time(startTimeMillis);
		System.out.println("Start: " + now.toString());

		double fmmax = 0.6519142957499121;
//		int 	bestTempoRangeLow = 60;
//		int		bestTempoRangeHigh = 200;
		int		bestLocalWindowSize = 5;
		int		bestBeatWindowFactor = 10;
		double	bestOctaveTolerance = 0.30483989214210144;

		while (tries < 7000) {
			tries++;
			if (tries % 100 == 0) {
				now.setTime(System.currentTimeMillis() - startTimeMillis);
				System.out.println(now.toString() + " #" + tries);
			}

//			int deltaTempoRangeLow = rand.nextInt(3) - 1;	// -1 bis +1
//			int deltaTempoRangeHigh = rand.nextInt(5) - 2;	// -2 bis +2
			int deltaLocalWindowSize = rand.nextInt(7) - 3;			// -3 bis +3
			int deltaBeatWindowFactor = rand.nextInt(3) - 1;			// -1 bis +1
			double deltaOctaveTolerance = rand.nextDouble() * 0.5 - 0.025;	// -0.025 bis +0.025

//			int tempoRangeLow = bestTempoRangeLow + deltaTempoRangeLow;
//			int	tempoRangeHigh = bestTempoRangeHigh + deltaTempoRangeHigh;
			int	localWindowSize = bestLocalWindowSize + deltaLocalWindowSize;
			int	beatWindowFactor = bestBeatWindowFactor + deltaBeatWindowFactor;
			double octaveTolerance = bestOctaveTolerance + deltaOctaveTolerance;

//			Beta.BT_TEMPO_RANGE_LOW = tempoRangeLow;
//			Beta.BT_TEMPO_RANGE_HIGH = tempoRangeHigh;
//			Beta.BT_LOCAL_WINDOWSIZE = localWindowSize;
//			Beta.BT_BEAT_WINDOW_FACTOR = beatWindowFactor;
//			Beta.BT_OCTAVE_TOLERANCE = octaveTolerance;

			processor = findAndInstantiateClass(processorName);

			if (options.has("p")) {
				//System.out.println("Predicting ...");
				Files.walk(Paths.get(directory))
						.filter(path -> Files.isRegularFile(path))
						.filter(path -> path.toString().endsWith(AUDIO))
						//.limit(1)	// TODO hack for only first file for testing
						.forEach(path -> predictForFile(path.toString()));
			}

			if (options.has("e")) {
				//System.out.println("Evaluating ...");
				Files.walk(Paths.get(directory))
						.filter(path -> Files.isRegularFile(path))
						.filter(path -> path.toString().endsWith(AUDIO))
						.forEach(path -> evalForFile(path.toString()));
			}

			if (options.has("s")) {
				//System.out.println("Summarizing ...");
				String type = BEATS;
				Map<String, Integer> summary = new HashMap<>();
				summary.put("tp", 0);
				summary.put("fp", 0);
				summary.put("fn", 0);
				summary.put("error", 0);

				Files.walk(Paths.get(directory))
						.filter(path -> Files.isRegularFile(path))
						.filter(path -> path.toString().endsWith(type + EV))
						.forEach(path -> updateSummaryForFile(summary, path.toString()));

				double precision = 0d;
				double recall = 0d;
				double fmeasure = 0d;

				int tp = summary.get("tp");
				int fp = summary.get("fp");
				int fn = summary.get("fn");
				double error = (double) summary.get("error") / 1e6d;

				if (tp + fp > 0)
					precision = (double) tp / (tp + fp);

				if (tp + fn > 0)
					recall = (double) tp / (tp + fn);

				if (precision + recall > 0)
					fmeasure = (2 * precision * recall) / (precision + recall);
				if (type.equalsIgnoreCase(BEATS) && fmeasure > fmmax) {
					fmmax = fmeasure;
					System.out.println("##################################################");
					System.out.println("new best FMeasure: " + fmeasure);
					System.out.println("precision: " + precision);
					System.out.println("recall: " + recall);
					System.out.println("error: " + error);
					System.out.println("parameters:");
//						System.out.println(tempoRangeLow);
//						System.out.println(tempoRangeHigh);
					System.out.println(localWindowSize);
					System.out.println(beatWindowFactor);
					System.out.println(octaveTolerance);
					System.out.println("##################################################");
//						bestTempoRangeLow = tempoRangeLow;
//						bestTempoRangeHigh = tempoRangeHigh;
					bestLocalWindowSize = localWindowSize;
					bestBeatWindowFactor = beatWindowFactor;
					bestOctaveTolerance = octaveTolerance;
				}
//					else if (fmeasure >= 0.995*fmmax) {
//						System.out.println("--------------------------------------------------");
//						System.out.println("very good FMeasure: " + fmeasure);
//						System.out.println("precision: " + precision);
//						System.out.println("recall: " + recall);
//						System.out.println("error: " + error);
//						System.out.println("f:");
//						System.out.println(f1);
//						System.out.println(f2);
//						System.out.println(f3);
//						System.out.println(f4);
//						System.out.println(f5);
//						System.out.println("--------------------------------------------------");
//					}
				else {
					//System.out.println("fmeasure: " + fmeasure);
				}
				lines.add(fmeasure + ";" + precision + ";" + recall + ";" + error + ";" + tp + ";" + fp + ";" + fn + ";" /*+ tempoRangeLow + ";" + tempoRangeHigh + ";"*/ + localWindowSize + ";" + beatWindowFactor + ";" + octaveTolerance);
			}
		}

		now.setTime(System.currentTimeMillis());
		System.out.println("Ende: " + now.toString());

		Utils.writeToFile(directory + "\\montecarlo.csv", lines);
	}

    private static Processor findAndInstantiateClass(String processorName) {
        String fullyQualifiedName = "at.jku.cp.spezi." + processorName;
        try {
            Class<?> clazz = Class.forName(fullyQualifiedName);
            return (Processor) clazz.newInstance();
        } catch (ClassNotFoundException e) {
            throw new RuntimeException(String.format("Processor name '%s' unknown ...", fullyQualifiedName));
        } catch (InstantiationException e) {
            throw new RuntimeException(String.format("Unable to instantiate processor '%s' ...", fullyQualifiedName));
        } catch (IllegalAccessException e) {
            throw new RuntimeException(String.format("Not allowed to instantiate processor '%s' ...", fullyQualifiedName));
        }
    }

    private static void updateSummaryForFile(Map<String, Integer> summary, String filename) {
        Map<String, Integer> dict = Utils.dictFromFile(filename);
        summary.put("tp", summary.get("tp") + dict.getOrDefault("tp", 0));
        summary.put("fp", summary.get("fp") + dict.getOrDefault("fp", 0));
        summary.put("fn", summary.get("fn") + dict.getOrDefault("fn", 0));
        summary.put("error", summary.get("error") + dict.getOrDefault("error", 0));
    }

    private static void predictForFile(String filename) {
        String prefix = filename.substring(0, filename.lastIndexOf("."));

        String onsets = prefix + ONSETS + PR;
        String beats = prefix + BEATS + PR;
        String tempo = prefix + TEMPO + PR;

        processor.process(filename);

        //System.out.println("Outputting Onset Times to " + onsets);
        Utils.writeDataToFile(onsets, processor.getOnsets());

        //System.out.println("Outputting Beat Times to " + beats);
        Utils.writeDataToFile(beats, processor.getBeats());

        //System.out.println("Outputting Tempo to " + tempo);
        Utils.writeDataToFile(tempo, processor.getTempo());
    }

    private static void evalForFile(String filename) {
        String prefix = filename.substring(0, filename.lastIndexOf("."));

        for (int i = 0; i < types.length; i++) {
            String type = types[i];
            double tolerance = tolerances[i];

            String gt = prefix + type + GT;
            String pr = prefix + type + PR;
            if (Files.exists(Paths.get(gt)) && Files.exists(Paths.get(pr))) {
                String ev = prefix + type + EV;
                //System.out.println("Evaluating '" + pr + "' ...");
                evaluateEventFiles(type, gt, pr, ev, tolerance);
            }
        }
    }

    /*
     * general evaluation function, that can be used to evaluate onsets, beats
     * and tempo in the form of event-lists
     */
    private static void evaluateEventFiles(
            String type,
            String groundtruth,
            String predictions,
            String evaluation,
            double tolerance) {

        List<Double> rawGT = Utils.listFromFile(groundtruth);
        List<Double> eventGT = Utils.cleanEventList(rawGT, tolerance);

        List<Double> rawPR = Utils.listFromFile(predictions);
        List<Double> eventPR = Utils.cleanEventList(rawPR, tolerance);

        // tolerance for tempo is in percent of the true tempo!
        if (type.equals(Runner.TEMPO)) {
            if (eventGT.size() != 1)
                throw new RuntimeException("weird tempo file ?!");
            double realTempo = eventGT.get(0);
            tolerance = (realTempo / 100d) * tolerance;
        }

        Map<String, Integer> summary = Utils.evaluateEventList(eventGT, eventPR, tolerance);
        Utils.dictToFile(evaluation, summary);
    }

}
