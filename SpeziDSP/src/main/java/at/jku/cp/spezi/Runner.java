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

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import at.jku.cp.spezi.dsp.Processor;
import joptsimple.OptionParser;
import joptsimple.OptionSet;

/**
 * @author andreas arzt
 * @author rainer kelz
 */
public class Runner {
	public static String ONSETS = ".onsets";
	public static String BEATS = ".beats";
	public static String TEMPO = ".tempo";

	private static String[] types = { ONSETS, BEATS, TEMPO };
	private static double[] tolerances = { 0.05d, 0.07d, 4d };

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

		System.out.println("Outputting Onset Times to " + onsets);
		Utils.writeDataToFile(onsets, processor.getOnsets());

		System.out.println("Outputting Beat Times to " + beats);
		Utils.writeDataToFile(beats, processor.getBeats());

		System.out.println("Outputting Tempo to " + tempo);
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
				System.out.println("Evaluating '" + pr + "' ...");
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
		if(type.equals(Runner.TEMPO))
		{
			if(eventGT.size() != 1)
				throw new RuntimeException("weird tempo file ?!");
			double realTempo = eventGT.get(0);
			tolerance = (realTempo / 100d) * tolerance;
		}
		
		Map<String, Integer> summary = Utils.evaluateEventList(eventGT, eventPR, tolerance);
		Utils.dictToFile(evaluation, summary);
	}

}
