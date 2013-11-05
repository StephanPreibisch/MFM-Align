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

import ij.CompositeImage;
import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import ij.process.FloatProcessor;
import io.ExtractPlane;
import io.OpenPiezoStack;
import io.TextFileAccess;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import process.AvgProjection3;
import process.CrossCorrelation;
import process.Mirror;
import process.OverlayFusion;
import process.QuantileNormalization;
import run.MicroscopyPlane.Mirroring;
import mpicbg.imglib.container.array.Array;
import mpicbg.imglib.container.array.ArrayContainerFactory;
import mpicbg.imglib.container.basictypecontainer.array.FloatArray;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.ImageFactory;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.interpolation.InterpolatorFactory;
import mpicbg.imglib.interpolation.linear.LinearInterpolatorFactory;
import mpicbg.imglib.io.ImageOpener;
import mpicbg.imglib.io.LOCI;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyFactory;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyMirrorFactory;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyValueFactory;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.models.AbstractAffineModel2D;
import mpicbg.models.IllDefinedDataPointsException;
import mpicbg.models.InvertibleBoundable;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.RigidModel2D;

public class Align 
{
	public static final PrintWriter outAllZ = null;//TextFileAccess.openFileWrite( "/Volumes/TOSHIBA EXT/3D analysis files/alloutZ.txt" );
	public static final PrintWriter outAllXY = null;//TextFileAccess.openFileWrite( "/Volumes/TOSHIBA EXT/3D analysis files/alloutXY.txt" );
	
	/**
	 * Will create a hyperstack with the following properties
	 * 
	 * t = piezo
	 * c = channel
	 * z = plane
	 * 
	 * @param planes
	 */
	public static CompositeImage showAlignedImages( final ArrayList< MicroscopyPlane > planes )
	{
		final ImageStack stack = new ImageStack( planes.get( 0 ).getImage().getDimension( 0 ), planes.get( 0 ).getImage().getDimension( 1 ) );

		// apply the model to the images
		for ( final MicroscopyPlane plane : planes )
		{
			final Image< FloatType > img = plane.getImage();			
			OverlayFusion.fuseChannel( img, img.clone(), new float[ img.getNumDimensions() ], plane.getXYModel(), new LinearInterpolatorFactory<FloatType>( new OutOfBoundsStrategyMirrorFactory<FloatType>() ) );
		}

		// create the hyperstack
		final int numImages = planes.get( 0 ).getImage().getDimension( 2 );
		final int numChannels = 2;
		final int numSlices = AlignProperties.numTiles;
		
		final int[] size = new int[]{ planes.get( 0 ).getImage().getDimension( 0 ), planes.get( 0 ).getImage().getDimension( 1 ) };
		final ImageFactory< FloatType > factory = planes.get( 0 ).getImage().getImageFactory();
		
		final OutOfBoundsStrategyFactory< FloatType > oobs = new OutOfBoundsStrategyValueFactory<FloatType>();
		final InterpolatorFactory< FloatType > interpolatorF = new LinearInterpolatorFactory<FloatType>( oobs );
		
		for ( int t = 0; t < numImages; ++t )
		{
			for ( int z = 0; z < numSlices; ++z )
			{
				// fuse
				for ( int c = 0; c < numChannels; ++c )
				{
					// which plane?
					final int index = z + numSlices*c;
					
					final MicroscopyPlane plane = planes.get( index );
					final Image< FloatType > img = plane.getImage();
					
					final Image< FloatType > planeTmp = factory.createImage( size );
					
					// extract the t'th piezo section
					//CrossCorrelation.extractPlane( img, planeTmp, t );
					CrossCorrelation.extractPlane( img.createInterpolator( interpolatorF ), planeTmp, t - plane.getModel().tx );
					
					final float[] pixels = ((FloatArray)((Array)planeTmp.getContainer()).update( null )).getCurrentStorageArray();
					stack.addSlice( plane.getFullName() + "_tile=" + plane.tileNumber + "_slice=" + t, new FloatProcessor( size[ 0 ], size[ 1 ], pixels ) );
				}
			}
		}
		
		//convertXYZCT ...
		ImagePlus result = new ImagePlus( "registered", stack );
				
		result.setDimensions( numChannels, numSlices, numImages );		
		//result = OverlayFusion.switchZCinXYCZT( result );
		return new CompositeImage( result, CompositeImage.COMPOSITE );
	}

	public static CompositeImage showAlignedProjections( final ArrayList< MicroscopyPlane > planes )
	{
		ImageStack stack = new ImageStack( planes.get( 0 ).getAvgProj().getDimension( 0 ), planes.get( 0 ).getAvgProj().getDimension( 1 ) );
		final int numChannels = 2;
		final int numSlices = AlignProperties.numTiles;

		for ( int z = 0; z < numSlices; ++z )
		{
			// fuse
			for ( int c = 0; c < numChannels; ++c )
			{
				// which plane?
				final int index = z + numSlices*c;
				final MicroscopyPlane plane = planes.get( index );
				
				final Image< FloatType > proj = plane.getAvgProj();
				OverlayFusion.fuseChannel( proj, proj.clone(), new float[ proj.getNumDimensions() ], plane.getXYModel(), new LinearInterpolatorFactory<FloatType>( new OutOfBoundsStrategyMirrorFactory<FloatType>() ) );			
				stack.addSlice( plane.getFullName(), ImageJFunctions.copyToImagePlus( proj ).getProcessor() );
			}
		}

		ImagePlus result = new ImagePlus( "registered", stack );
		result.setDimensions( 2, AlignProperties.numTiles, 1 );
		
		return new CompositeImage( result, CompositeImage.COMPOSITE );
	}
	
	public static CompositeImage createFinalImages( final ArrayList< MicroscopyPlane > planes, final String baseDir, final String[] target, final Mirroring[] mirror,
						final boolean adjust, final boolean quantile ) throws Exception, IOException
	{
		final File t1 = new File( baseDir, target[ 0 ] );
		final File t2 = new File( baseDir, target[ 1 ] );
		final File t3 = new File( baseDir, target[ 2 ] );
		
		System.out.println( t1.getAbsolutePath() );
		System.out.println( t2.getAbsolutePath() );
		System.out.println( t3.getAbsolutePath() );
			
		final Image< FloatType > img1 = new ImageOpener().openImage( t1.getAbsolutePath(), new ImageFactory<FloatType>( new FloatType(), new ArrayContainerFactory() ) );
		final Image< FloatType > img2 = new ImageOpener().openImage( t2.getAbsolutePath(), new ImageFactory<FloatType>( new FloatType(), new ArrayContainerFactory() ) );
		final Image< FloatType > img3 = new ImageOpener().openImage( t3.getAbsolutePath(), new ImageFactory<FloatType>( new FloatType(), new ArrayContainerFactory() ) );
		
		if ( mirror[ 0 ] == Mirroring.HORIZONTALLY )
			Mirror.horizontal( img1 );
		if ( mirror[ 1 ] == Mirroring.HORIZONTALLY )
			Mirror.horizontal( img2 );
		if ( mirror[ 2 ] == Mirroring.HORIZONTALLY )
			Mirror.horizontal( img3 );
		
		// extract all planes and overwrite them in the planes ArrayList
		for ( int i = 0; i < 9; ++i )
		{
			final MicroscopyPlane plane1 = planes.get( i );
			final MicroscopyPlane plane2 = planes.get( i + 9 );
			
			final MicroscopyPlane plane3 = new MicroscopyPlane( plane1.getBaseDirectory(), target[ 2 ], plane1.tag, mirror[ 2 ], i );					
			//new MicroscopyPlane( plane1.baseDir, target[ 2 ], mirror[ 2 ], i );
			
			plane1.tag = target[ 0 ];
			plane2.tag = target[ 1 ];
			plane1.localDir = "";
			plane2.localDir = "";
			
			plane1.setImage( ExtractPlane.extract( img1, i ) );
			plane2.setImage( ExtractPlane.extract( img2, i ) );
			plane3.setImage( ExtractPlane.extract( img3, i ) );
			
			// set the alignment of plane1
			plane3.setXYModel( ((InvertibleBoundable)((AbstractAffineModel2D)plane1.getXYModel()).copy()) );
			
			// apply the model to the images
			Image< FloatType > img = plane1.getImage();
			OverlayFusion.fuseChannel( img, img.clone(), new float[ img.getNumDimensions() ], plane1.getXYModel(), new LinearInterpolatorFactory<FloatType>( new OutOfBoundsStrategyMirrorFactory<FloatType>() ) );
			img = plane2.getImage();
			OverlayFusion.fuseChannel( img, img.clone(), new float[ img.getNumDimensions() ], plane2.getXYModel(), new LinearInterpolatorFactory<FloatType>( new OutOfBoundsStrategyMirrorFactory<FloatType>() ) );
			img = plane3.getImage();
			OverlayFusion.fuseChannel( img, img.clone(), new float[ img.getNumDimensions() ], plane3.getXYModel(), new LinearInterpolatorFactory<FloatType>( new OutOfBoundsStrategyMirrorFactory<FloatType>() ) );
			
			// compute the average projection of plane3
			final Image< FloatType > planeImg = AvgProjection3.project( plane3.getImage() );
			plane3.setAvgProj( planeImg );
			
			planes.add( plane3 );
		}
		
		if ( adjust )
		{
			// extract all planes and overwrite them in the planes ArrayList
			final double avg1 = CrossCorrelation.avg10( planes.get( 4 ).getImage() );
			final double avg2 = CrossCorrelation.avg10( planes.get( 4 + 9 ).getImage() );
			final double avg3 = CrossCorrelation.avg10( planes.get( 4 + 18 ).getImage() );
			
			for ( int i = 0; i < 9; ++i )
			{
				final MicroscopyPlane plane1 = planes.get( i );
				final MicroscopyPlane plane2 = planes.get( i + 9 );			
				final MicroscopyPlane plane3 = planes.get( i + 18 );
				
				final double diff1 = avg1 / CrossCorrelation.avg10( plane1.getImage() );
				final double diff2 = avg2 / CrossCorrelation.avg10( plane2.getImage() );
				final double diff3 = avg3 / CrossCorrelation.avg10( plane3.getImage() );
				
				for ( final FloatType t : plane1.getImage() )
					t.set( t.get() * (float)diff1 );
				
				for ( final FloatType t : plane2.getImage() )
					t.set( t.get() * (float)diff2 );
				
				for ( final FloatType t : plane3.getImage() )
					t.set( t.get() * (float)diff3 );
			}
		}
		
		final int numTimepoints = img1.getDimension( 2 );
		final int numSlices = 9;
		final int numChannels = 3; // mrna, npc, mrna dna
		final ImageStack stack = new ImageStack( planes.get( 0 ).getImage().getDimension( 0 ), planes.get( 0 ).getImage().getDimension( 1 ) );
		
		final int[] size = new int[]{ planes.get( 0 ).getImage().getDimension( 0 ), planes.get( 0 ).getImage().getDimension( 1 ) };
		final ImageFactory< FloatType > factory = planes.get( 0 ).getImage().getImageFactory();
		
		
		Image< FloatType >[] ref = null;
		
		if ( quantile )
		{
			ref = new Image[ numSlices * 2 ];
			
			for ( int z = 0; z < numSlices; ++z )
			{
				// fuse
				for ( int c = 0; c < 2; ++c )
				{
					// which plane?
					final int index = z + numSlices*c;
					
					final MicroscopyPlane plane = planes.get( index );
					final Image< FloatType > img = plane.getImage();
					
					// extract the 5'th piezo section
					final Image< FloatType > planeTmp = factory.createImage( size );
					CrossCorrelation.extractPlane( img, planeTmp, 4 );
					
					ref[ index ] = planeTmp;
				}
			}
		}
		
		
		for ( int t = 0; t < numTimepoints; ++t )
		{
			for ( int z = 0; z < numSlices; ++z )
			{
				// fuse
				for ( int c = 0; c < numChannels; ++c )
				{
					// which plane?
					final int index = z + numSlices*c;
					
					final MicroscopyPlane plane = planes.get( index );
					final Image< FloatType > img = plane.getImage();
					
					final Image< FloatType > planeTmp;
					
					// extract the t'th piezo section
					if ( c == 2 )
					{
						// the 3rd channel is always just copied
						planeTmp = plane.getAvgProj().clone();
					}
					else
					{
						planeTmp = factory.createImage( size );
						CrossCorrelation.extractPlane( img, planeTmp, t );
						
						if ( quantile )
							QuantileNormalization.normalizeTo( ref[ index ], planeTmp );
					}
					
					final float[] pixels = ((FloatArray)((Array)planeTmp.getContainer()).update( null )).getCurrentStorageArray();
					stack.addSlice( plane.getFullName() + "_time=" + t, new FloatProcessor( size[ 0 ], size[ 1 ], pixels ) );
				}
			}
		}
		//ref = ExtractPlane.extract( refGlobal, refIndex );
		
		ImagePlus result = new ImagePlus( "registered", stack );		
		result.setDimensions( numChannels, numSlices, numTimepoints );
		return new CompositeImage( result, CompositeImage.COMPOSITE );
	}

	public static void findAllDataDirs( final String base, final ArrayList< String > data )
	{
		final File f = new File( base );
		
		String[] entries = f.list();
		
		for ( final String entry : entries )
		{
			if ( entry.contains( "mRNAreg.tif") )
			{
				data.add( f.getAbsolutePath() );
				return;
			}
		}
		
		for ( final String entry : entries )
		{
			final File newF = new File( base, entry );
			
			if ( newF.isDirectory() )
				findAllDataDirs( newF.getAbsolutePath(), data );
		}
	}
	
	public static void alignAll( final String baseDir, final AbstractAffineModel2D< ? > model2d ) throws Exception
	{		
		final String localDir = "DNA stack";
			
		final String[] tags = new String[] { "2464" /*green*/, "4283" /*red*/ };
		final Mirroring[] mirror = new Mirroring[]{ Mirroring.HORIZONTALLY, Mirroring.DONOT };
		
		//
		// set up the planes
		// 	
		final ArrayList< MicroscopyPlane > planes = new ArrayList< MicroscopyPlane >();
		
		for ( int c = 0; c < tags.length; ++c )
			for ( int t = 0; t < AlignProperties.numTiles; ++t )
				planes.add( new MicroscopyPlane( baseDir, localDir, tags[ c ], mirror[ c ], t ) );
		
		AlignZ alignZ = new AlignZ( planes );
		AlignXY alignXY = new AlignXY( alignZ.getPlanes(), model2d );
		
		// apply models to the avg projections
		showAlignedProjections( alignZ.getPlanes() ).show();
		
		// apply to the images
		// showAlignedImages( alignZ.getPlanes() ).show();

		final String greenChannelLarge = findImageFile( new File( baseDir ), "2464", 200*1024*1024, 500*1024*1024 );
		final String greenChannelDNA = findImageFile( new File( baseDir ), "2464", 5*1024*1024, 15*1024*1024 );
		final String redChannelLarge = findImageFile( new File( baseDir ), "4283", 200*1024*1024, 500*1024*1024 );

		final String[] target = new String[]{ greenChannelLarge, redChannelLarge, greenChannelDNA };
		final Mirroring[] mirrorTarget = new Mirroring[]{ Mirroring.HORIZONTALLY, Mirroring.DONOT, Mirroring.HORIZONTALLY };
		
		CompositeImage ci = createFinalImages( alignZ.getPlanes(), baseDir, target, mirrorTarget, false, false );
		FileSaver fs = new FileSaver( ci );
		fs.saveAsTiffStack( new File( baseDir, "raw_aligned.tif" ).getAbsolutePath() );
		ci.close();
		
		ci = createFinalImages( alignZ.getPlanes(), baseDir, target, mirrorTarget, true, false );
		fs = new FileSaver( ci );
		fs.saveAsTiffStack( new File( baseDir, "avgcorrected_aligned.tif" ).getAbsolutePath() );
		ci.close();

		ci = createFinalImages( alignZ.getPlanes(), baseDir, target, mirrorTarget, true, true );
		fs = new FileSaver( ci );
		fs.saveAsTiffStack( new File( baseDir, "avgcorrected_quantile_aligned.tif" ).getAbsolutePath() );
		ci.close();
	}
	
	/**
	 * Find file that corresponds to these criterian
	 * 
	 * @param baseDir - the directory
	 * @param tag - a part of the filename
	 * @param minSize - min size in byte
	 * @param maxSize - max size in byte
	 * @return
	 */
	public static String findImageFile( final File baseDir, final String tag, final int minSize, final int maxSize )
	{
		final String[] list = baseDir.list(
				new FilenameFilter() 
				{	
					@Override
					public boolean accept( final File dir, final String name ) 
					{
						final File newFile = new File( dir, name );
						
						// ignore directories and hidden files
						if ( newFile.isHidden() || newFile.isDirectory() )
							return false;
						else if ( name.contains( tag ) && newFile.length() > minSize && newFile.length() < maxSize )
							return true;
						else
							return false;
					}
				}
				);
	
		if ( list.length == 0 )
		{
			System.out.println( "No file found in '" + baseDir + "' that contains '" + tag + "' and is larger than '" + minSize + "' and smaller than '" + maxSize + "'" );
			return null;
		}
		else if ( list.length > 1 )
		{
			System.out.println( "MORE that one file found in '" + baseDir + "' that contains '" + tag + "' and is larger than '" + minSize + "' and smaller than '" + maxSize + "'" );

			for ( final String s : list )
			{
				final File f = new File( baseDir, s );
				System.out.println( f + " >> " + f.length() );
			}
			return null;
		}
		else
		{
			return list[ 0 ];
		}
	}
	
	public static void main( String[] args ) throws Exception
	{
		final AbstractAffineModel2D< ? > model2d = new RigidModel2D();
		
		final String root = "/home/stephanpreibisch/Desktop/stephan/";
		final String experimentDir = "1 (20110525, dish 2, cell 22)";

		new ImageJ();
		
		final ArrayList< String > allDirs = new ArrayList<String>();
		findAllDataDirs( root, allDirs );
	
		allDirs.clear();
		//allDirs.add( "/Users/preibischs/Documents/Microscopy/david/sample 1/" );
		//allDirs.add( "/Volumes/TOSHIBA EXT/3D analysis files/2011-05-25/Dish 03/cell 06" );
		allDirs.add( root + experimentDir );
		
		for ( final String baseDir : allDirs )
		{
			System.out.println( baseDir );
			
			try
			{
				alignAll( baseDir, model2d.copy() );
			}
			catch (Exception e)
			{
				System.err.println( "Failed on: " + baseDir );
			}
			
			System.out.println( "done." );
		}
		
		if ( outAllZ != null )
			outAllZ.close();
		
		if ( outAllXY != null )
			outAllXY.close();
		
		System.out.println( "All done. exiting." );
		
		System.exit( 0 );
	}

}
