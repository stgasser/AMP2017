package at.jku.cp.spezi.dsp;

import java.util.List;

public interface Processor {
	void process(String filename);
	
	List<Double> getOnsets();
	List<Double> getTempo();
	List<Double> getBeats();
}