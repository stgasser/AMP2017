package at.jku.cp.spezi;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

public class TestEval extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public TestEval( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( TestEval.class );
    }
    
    public void testRoundtripLines() throws IOException
    {
    	File temp = File.createTempFile("testRoundtripLines_", ".txt");
        temp.deleteOnExit();
        
        List<String> expected = new ArrayList<>();
        expected.add("2389746293746076p91827383 20398 ritsn");
        expected.add("2398   9238     esnthrnsth *$^%^&(*&%^&*$%^ ");
        expected.add("a");
        
        Utils.writeToFile(temp.getAbsolutePath(), expected);
        
        List<String> actual = Files.readAllLines(Paths.get(temp.getAbsolutePath()));
        assertEquals(expected, actual);
    }

    public void testRoundtripDict() throws IOException
    {
    	File temp = File.createTempFile("testRoundtripDict_", ".dict");
        temp.deleteOnExit();
        
        Map<String, Integer> expected = new HashMap<>();
        expected.put("x", 0);
        expected.put("!@#$%^", 1);
        expected.put("002182$%^&*()(*&^%^&*()", 2);
     
        Utils.dictToFile(temp.getAbsolutePath(), expected);
        
        Map<String, Integer> actual = Utils.dictFromFile(temp.getAbsolutePath());
        assertEquals(expected, actual);
    }
    
    public void testEvalEvents() {
    	List<Double> gt = new ArrayList<>();
    	gt.add(0d);
    	gt.add(1d);
    	gt.add(2d);
    	
    	
    	List<Double> pr = new ArrayList<>();
    	pr.add(0.09d); // tp
    	pr.add(0.91d); // tp
    	pr.add(2.2d);  // fp
    	// and we missed the last entry from gt, so one fn as well!
    	
    	Map<String, Integer> summary = Utils.evaluateEventList(gt, pr, 0.1d);
    	System.out.println(summary);
    	assertEquals((int) 2, (int) summary.get("tp"));
    	assertEquals((int) 1, (int) summary.get("fp"));
    	assertEquals((int) 1, (int) summary.get("fn"));
    }
    
    public void testEvalEventsProduceCorrectCountsSimple() {
    	List<Double> gt = new ArrayList<>();
    	gt.add(0d);
    	gt.add(1d);
    	gt.add(2d);
    	
    	
    	List<Double> pr = new ArrayList<>();
    	pr.add(0.09d); // tp
    	pr.add(0.91d); // tp
    	pr.add(2.2d);  // fp
    	// and we missed the last entry from gt, so one fn as well!
    	
    	Map<String, Integer> summary = Utils.evaluateEventList(gt, pr, 0.1d);
    	System.out.println(summary);
    	assertEquals((int) 2, (int) summary.get("tp"));
    	assertEquals((int) 1, (int) summary.get("fp"));
    	assertEquals((int) 1, (int) summary.get("fn"));
    }

    public void testEvalEventsProduceCorrectCountsTooCrazy() {
    	List<Double> gt = new ArrayList<>();
    	gt.add(0d);
    	gt.add(1d);
    	gt.add(2d);
    	
    	
    	List<Double> pr = new ArrayList<>();
    	pr.add(0.09d); // tp
    	pr.add(0.91d); // tp
    	pr.add(2.2d);  // fp
    	// and we missed the last entry from gt, so one fn as well!
    	
    	Map<String, Integer> summary = Utils.evaluateEventList(gt, pr, 0.1d);
    	System.out.println(summary);
    	assertEquals((int) 2, (int) summary.get("tp"));
    	assertEquals((int) 1, (int) summary.get("fp"));
    	assertEquals((int) 1, (int) summary.get("fn"));
    }

    
    public void testEventListCleaningComputesMeanOfTooCloseEvents() {
    	List<Double> unclean = new ArrayList<>();
    	unclean.add(0d);
    	unclean.add(1d);
    	unclean.add(1.5d);
    	unclean.add(2d);
    	
    	
    	List<Double> expected_clean = new ArrayList<>();
    	expected_clean.add(0d);
    	expected_clean.add(1.25d);
    	expected_clean.add(2d);
    	
    	List<Double> actual_clean = Utils.cleanEventList(unclean, 0.6d);
    	assertEquals(expected_clean, actual_clean);
    }
    
    public void testEventListCleaningForEmptyEventList() {
    	List<Double> unclean = new ArrayList<>();
    	List<Double> expected_clean = new ArrayList<>();
    	
    	List<Double> actual_clean = Utils.cleanEventList(unclean, 0.6d);
    	assertEquals(expected_clean, actual_clean);
    }
}
