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
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import process.AvgProjection3;
import process.CrossCorrelation;
import process.Mirror;
import process.OverlayFusion;
import process.QuantileNormalization;

import loci.formats.FormatException;
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
					stack.addSlice( plane.name + "_tile=" + plane.tileNumber + "_slice=" + t, new FloatProcessor( size[ 0 ], size[ 1 ], pixels ) );
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
				stack.addSlice( plane.name + "_tile=" + plane.tileNumber, ImageJFunctions.copyToImagePlus( proj ).getProcessor() );
			}
		}

		ImagePlus result = new ImagePlus( "registered", stack );
		result.setDimensions( 2, AlignProperties.numTiles, 1 );
		
		return new CompositeImage( result, CompositeImage.COMPOSITE );
	}
	
	public static CompositeImage createFinalImages( final ArrayList< MicroscopyPlane > planes, final String baseDir, final String[] target, final boolean[] mirror, final boolean adjust, final boolean quantile ) throws FormatException, IOException
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
		
		if ( mirror[ 0 ] )
			Mirror.horizontal( img1 );
		if ( mirror[ 1 ] )
			Mirror.horizontal( img2 );
		if ( mirror[ 2 ] )
			Mirror.horizontal( img3 );
		
		// extract all planes and overwrite them in the planes ArrayList
		for ( int i = 0; i < 9; ++i )
		{
			final MicroscopyPlane plane1 = planes.get( i );
			final MicroscopyPlane plane2 = planes.get( i + 9 );
			
			final MicroscopyPlane plane3 = new MicroscopyPlane( plane1.baseDir, target[ 2 ], mirror[ 2 ], i );
			
			plane1.name = target[ 0 ];
			plane2.name = target[ 1 ];
			
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
					stack.addSlice( plane.name + "_tile=" + plane.tileNumber + "_time=" + t, new FloatProcessor( size[ 0 ], size[ 1 ], pixels ) );
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
	
	public static void alignAll(String baseDir) throws FormatException, IOException, NotEnoughDataPointsException, IllDefinedDataPointsException
	{
		final String[] names = new String[]{ "DNA stack mRNA", "DNA stack NPC" };
		final boolean[] mirror = new boolean[]{ true, false };
		final String[] target = new String[]{ "mRNA.tif", "NPC.tif", "mRNAreg.tif" };
		final boolean[] mirrorTarget = new boolean[]{ true, false, true };
		
		AlignZ alignZ = new AlignZ( baseDir, names, mirror );
		AlignXY alignXY = new AlignXY( alignZ.getPlanes() );
		
		// apply models to the avg projections
		showAlignedProjections( alignZ.getPlanes() ).show();
		
		// apply to the images
		// showAlignedImages( alignZ.getPlanes() ).show();
		
		/*
		CompositeImage ci = createFinalImages( alignZ.getPlanes(), baseDir, target, mirrorTarget, false );
		FileSaver fs = new FileSaver( ci );
		fs.saveAsTiffStack( new File( baseDir, "raw_aligned.tif" ).getAbsolutePath() );
		ci.close();
		*/
		/*
		CompositeImage ci = createFinalImages( alignZ.getPlanes(), baseDir, target, mirrorTarget, true, false );
		FileSaver fs = new FileSaver( ci );
		fs.saveAsTiffStack( new File( baseDir, "avgcorrected_aligned.tif" ).getAbsolutePath() );
		ci.close();
		*/
		/*
		CompositeImage ci = createFinalImages( alignZ.getPlanes(), baseDir, target, mirrorTarget, true, true );
		FileSaver fs = new FileSaver( ci );
		fs.saveAsTiffStack( new File( baseDir, "avgcorrected_quantile_aligned.tif" ).getAbsolutePath() );
		ci.close();
		*/	
	}
	
	public static void main( String[] args ) throws FormatException, IOException, NotEnoughDataPointsException, IllDefinedDataPointsException
	{

		new ImageJ();
		
		final ArrayList< String > allDirs = new ArrayList<String>();
		//findAllDataDirs( "/Volumes/TOSHIBA EXT/3D analysis files", allDirs );
	
		allDirs.clear();
		//allDirs.add( "/Users/preibischs/Documents/Microscopy/david/sample 1/" );
		//allDirs.add( "/Volumes/TOSHIBA EXT/3D analysis files/2011-05-25/Dish 03/cell 06" );
		allDirs.add( "/Users/preibischs/Documents/Microscopy/david/sample 1/" );
		
		for ( final String baseDir : allDirs )
		{
			System.out.println( baseDir );
			
			try
			{
				alignAll( baseDir );
			}
			catch (Exception e)
			{
				System.err.println( "Failed on: " + baseDir );
			}
			
			System.out.println( "done." );
		}
		
		outAllZ.close();
		outAllXY.close();
		
		System.exit( 0 );
	}

}
