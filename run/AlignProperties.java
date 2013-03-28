package run;

public class AlignProperties 
{
	public final static int numTiles = 9;
	public static double epsilon = 0.2;
	public static double minInlierRatio = 0.5;
	
	public static String tmpName = "tmp_";
	public static String piezoStack = "_piezo.tif";
	public static String piezoProj = "_piezo_avg.tif";
	public static String entropies = "_entropies.txt";
	
	//public static double[] sigma = new double[]{ 0.75, 0.75, 4 };
	public static double[] sigma = new double[]{ 0, 0, 1 };

}
