package run;

import java.util.ArrayList;

import fit.Line;
import fit.LinkedPoint;
import fit.PointFunctionMatch;
import fit.TranslationModel1D;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.models.IllDefinedDataPointsException;
import mpicbg.models.InvertibleBoundable;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.Point;
import mpicbg.models.Tile;

public class MicroscopyPlane extends Tile< TranslationModel1D >
{
	final String baseDir;
	String name;
	final boolean mirror;
	final int tileNumber;
	
	Image< FloatType > image = null, avgProj = null;
	InvertibleBoundable model = null;
	
	public Image< FloatType > referenceStack = null;
	
	// offsets within the same channel
	final ArrayList< PlaneOffset > offsets1 = new ArrayList< PlaneOffset >();
				
	// offsets within the other channel
	final ArrayList< PlaneOffset > offsets2 = new ArrayList< PlaneOffset >();
				
	public MicroscopyPlane( final String baseDir, final String name, final boolean mirror, final int tileNumber )
	{
		super( new TranslationModel1D() );
		
		this.baseDir = baseDir;
		this.name = name;
		this.mirror = mirror;
		this.tileNumber = tileNumber;
	}
	
	public void setImage( final Image< FloatType > image ) { this.image = image; }
	public Image< FloatType > getImage() { return image; }

	public void setAvgProj( final Image< FloatType > avgProj ) { this.avgProj = avgProj; }
	public Image< FloatType > getAvgProj() { return avgProj; }

	public void setXYModel( final InvertibleBoundable model ) { this.model = model; }
	public InvertibleBoundable getXYModel() { return model; }

	public int getTileNumber() { return tileNumber; }
	
	public static int removeOutliers( final ArrayList< PlaneOffset > offsets )
	{
		return removeOutliers( offsets, 0.1, 0.5 );
	}
	
	public static int removeOutliers( final ArrayList< PlaneOffset > offsets, final double epsilon, final double minInlierRatio )
	{
		final ArrayList< PointFunctionMatch > candidates = new ArrayList<PointFunctionMatch>();
		final ArrayList< PointFunctionMatch > inliers = new ArrayList<PointFunctionMatch>();
		
		for ( final PlaneOffset p : offsets )
			candidates.add( new PointFunctionMatch( new LinkedPoint<PlaneOffset>( new float[]{ p.plane.getTileNumber(), p.offset }, p )));//new PointFunctionMatch( p ) );
		
		int numRemoved = 0;
		
		try 
		{
			final Line l = new Line();
			
			l.ransac( candidates, inliers, 100, epsilon, minInlierRatio );
			numRemoved = candidates.size() - inliers.size();
			
			
			l.fit( inliers );
			
			System.out.println( "y = " + l.getM() + " x + " + l.getN() );
			
			offsets.clear();
			
			for ( final PointFunctionMatch p : inliers )
			{
				final LinkedPoint< ? > lp = (LinkedPoint< ? >)p.getP1();
				offsets.add( (PlaneOffset)lp.getLinkedObject() );
			}
			
			//for ( final PointFunctionMatch p : inliers )
			//	System.out.println( l.distanceTo( p.getP1() ) );

		} 
		catch (Exception e) 
		{
			e.printStackTrace();
			numRemoved = offsets.size();
			offsets.clear();
		} 
		
		return numRemoved;
	}
	
}
