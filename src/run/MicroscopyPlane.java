/**
 * License: GPL
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License 2
 * as published by the Free Software Foundation.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * 
 * @author: Stephan Preibisch (stephan.preibisch@gmx.de)
 */
package run;

import ij.io.FileSaver;
import io.ExtractPlane;
import io.OpenPiezoStack;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import loci.formats.FormatException;
import mpicbg.imglib.container.array.ArrayContainerFactory;
import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.io.LOCI;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.imglib.util.Util;
import mpicbg.models.InvertibleBoundable;
import mpicbg.models.Tile;
import process.Mirror;
import fit.Line;
import fit.LinkedPoint;
import fit.PointFunctionMatch;
import fit.TranslationModel1D;

/**
 * Class that represents an individual z plane as acquired by the microscope.
 * 
 * @author preibischs
 *
 */
public class MicroscopyPlane extends Tile< TranslationModel1D >
{
	public static enum Mirroring { DONOT, HORIZONTALLY };
	
	final String baseDir;
	String localDir, tag;
	String darkCountImageName;
	final Mirroring mirror;
	final int tileNumber;
	
	Image< FloatType > image = null, avgProj = null;
	InvertibleBoundable model = null;
	
	public Image< FloatType > referenceStack = null;
	
	// offsets within the same channel
	final ArrayList< PlaneOffset > offsetsSameChannel = new ArrayList< PlaneOffset >();
				
	// offsets within the other channel
	final ArrayList< PlaneOffset > offsetsOtherChannel = new ArrayList< PlaneOffset >();
		
	/**
	 * Instantiate a new {@link MicroscopyPlane} object
	 * 
	 * @param baseDir - base directory of the MFM experiment
	 * @param dirname - the directory that contains all the individual slices of the piezo DNA stack
	 * @param tag - if multiple channels are in the directory a String that selects for the current channel (e.g. green)
	 * @param mirror - mirror the image or not
	 * @param tileNumber - which of the tiles to load (0...8)
	 */
	public MicroscopyPlane( final String baseDir, final String dirname, final String tag, final String darkCountImageName, final Mirroring mirror, final int tileNumber )
	{
		super( new TranslationModel1D() );
		
		this.baseDir = baseDir;
		this.localDir = dirname;
		this.tag = tag;
		this.mirror = mirror;
		this.tileNumber = tileNumber;
		this.darkCountImageName = darkCountImageName;
	}
	
	public String getDarkCountImageName() { return darkCountImageName; }
	public void setDarkCountImageName( final String darkCountImageName ) { this.darkCountImageName = darkCountImageName; }
	public String getBaseDirectory() { return baseDir; }
	public String getFullName() { return tag + "_" + localDir + "_" + tileNumber; }
	public String getLocalDirectory() { return localDir; }
	public String getTagName() { return tag; }
	public Mirroring getMirror() { return mirror; }
	
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
	
	@Override
	public String toString()
	{
		return getFullName();
	}
	
	public Image<FloatType> getImagePiezo() throws FormatException, IOException
	{
		return getImagePiezo( this );
	}
	
	public static Image<FloatType> getImagePiezo( final MicroscopyPlane plane ) throws FormatException, IOException
	{	
		Image<FloatType> image;
		
		final File file = new File( plane.getBaseDirectory(), AlignProperties.tmpName + plane.getFullName() + AlignProperties.piezoStack );
				
		// load or create the 3d-stack
		if ( file.exists() )
		{
			image = LOCI.openLOCIFloatType( file.getAbsolutePath(), new ArrayContainerFactory() );
		}
		else
		{
			image = OpenPiezoStack.openPiezo( new File( plane.getBaseDirectory(), plane.getLocalDirectory() ), plane.getTagName() );
			
			if ( subtractDarkCount( image, plane.getDarkCountImageName() ) )
				System.out.println( "SUBTRACTED darkcount image '" + plane.getDarkCountImageName() + "'" );
			else
				System.out.println( "NOT FOUND Darkcount image '" + plane.getDarkCountImageName() + "'" );
			
			if ( plane.getMirror() == Mirroring.HORIZONTALLY )
				Mirror.horizontal( image );
			
			image = ExtractPlane.extract( image, plane.getTileNumber() );
			
			// save the extracted stack
			FileSaver fs = new FileSaver( ImageJFunctions.copyToImagePlus( image ) );
			fs.saveAsTiffStack( file.getAbsolutePath() );
		}
		
		// load or create the average-projection
		return image;
	}
	
	public static boolean subtractDarkCount( final Image< FloatType > image, final String darkCountFileName )
	{
		if ( darkCountFileName != null && new File( darkCountFileName ).exists() )
		{
			final Image< FloatType > darkCount = LOCI.openLOCIFloatType( new File( darkCountFileName ).getAbsolutePath(), new ArrayContainerFactory() );
			
			if ( darkCount.getDimension( 0 ) != image.getDimension( 0 ) || darkCount.getDimension( 1 ) != image.getDimension( 1 ) )
			{
				System.out.println( "Dimensionality for dark count subtraction does not match:" );
				System.out.println( "image:" + Util.printCoordinates( image.getDimensions() ) );
				System.out.println( "darkcounts:" + Util.printCoordinates( darkCount.getDimensions() ) );
				
				return false;
			}
			
			final Cursor< FloatType > imageCursor = image.createCursor();
			final Cursor< FloatType > darkCountCursor = darkCount.createCursor();
			
			// subtract plane - by - plane
			while ( imageCursor.hasNext() )
			{
				imageCursor.fwd();
				
				if ( !darkCountCursor.hasNext() )
					darkCountCursor.reset();
				
				darkCountCursor.fwd();
				
				imageCursor.getType().set( Math.max( 0, imageCursor.getType().get() - darkCountCursor.getType().get() ) );
			}
			
			return true;			
		}
		else
		{
			return false;
		}		
	}
}
