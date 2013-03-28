package run;

import ij.ImageJ;
import ij.io.FileSaver;
import io.ExtractPlane;
import io.OpenPiezoStack;
import io.TextFileAccess;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;

import loci.formats.FormatException;
import mpicbg.imglib.algorithm.gauss.GaussianConvolutionReal;
import mpicbg.imglib.container.array.ArrayContainerFactory;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.io.LOCI;
import mpicbg.imglib.multithreading.SimpleMultiThreading;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyMirrorFactory;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.models.IllDefinedDataPointsException;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.Point;
import mpicbg.models.PointMatch;
import mpicbg.models.TileConfiguration;
import process.Alignment;
import process.AutoFocus;
import process.ComputeEntropy;
import process.CrossCorrelation;
import process.Mirror;

public class AlignZ 
{
	final ArrayList< MicroscopyPlane > planes;
	final HashMap< String, Image<FloatType> > piezoStacks = new HashMap<String, Image<FloatType>>();	
	
	public AlignZ( final String baseDir, final String[] names, final boolean[] mirror ) throws FormatException, IOException, NotEnoughDataPointsException, IllDefinedDataPointsException
	{
		//
		// set up the planes
		// 	
		planes = new ArrayList< MicroscopyPlane >();
		
		for ( int c = 0; c < names.length; ++c )
			for ( int t = 0; t < AlignProperties.numTiles; ++t )
				planes.add( new MicroscopyPlane( baseDir, names[ c ], mirror[ c ], t ) );
		
		//
		// Compute all pairwise matches
		//		
		for ( final MicroscopyPlane planeA : planes )
		{
			// get all corresponding plane offsets
			for ( final MicroscopyPlane planeB : planes )
			{
				if ( planeA != planeB )
				{
					final float offset = computePairwiseAlignment( baseDir, planeA, planeB );
					
					if ( planeB.name.equals( planeA.name ) )
						planeA.offsets1.add( new PlaneOffset( planeB, offset ) );
					else
						planeA.offsets2.add( new PlaneOffset( planeB, offset ) );					
				}
			}
			
			// compute RANSAC to filter the correct ones
		
			//for ( final PlaneOffset planeOffset : planeA.offsets1 )
			//	System.out.println( "same " + planeA.name + " " + planeA.tileNumber + " <-> " + planeOffset.plane.name + " " + planeOffset.plane.tileNumber + " = " + planeOffset.offset );

			
			
			int numRemoved = MicroscopyPlane.removeOutliers( planeA.offsets1, AlignProperties.epsilon, AlignProperties.minInlierRatio );
			System.out.println( "same " + planeA.name + " " + planeA.tileNumber + " -> removed " + numRemoved );
			
			//for ( final PlaneOffset planeOffset : planeA.offsets2 )
			//	System.out.println( "diff " + planeA.name + " " + planeA.tileNumber + " <-> " + planeOffset.plane.name + " " + planeOffset.plane.tileNumber + " = " + planeOffset.offset );

			numRemoved = MicroscopyPlane.removeOutliers( planeA.offsets2, AlignProperties.epsilon, AlignProperties.minInlierRatio );
			System.out.println( "diff " + planeA.name + " " + planeA.tileNumber + " -> removed " + numRemoved );
			
			//for ( final PlaneOffset planeOffset : planeA.offsets2 )
			//	System.out.println( "diff " + planeA.name + " " + planeA.tileNumber + " <-> " + planeOffset.plane.name + " " + planeOffset.plane.tileNumber + " = " + planeOffset.offset );

			final ArrayList< PlaneOffset > allPlanes = new ArrayList<PlaneOffset>();
			allPlanes.addAll( planeA.offsets1 );
			allPlanes.addAll( planeA.offsets2 );
			
			// add them to the tile configuration
			for ( final PlaneOffset offset : allPlanes )
			{
				final Point p1 = new Point( new float[]{ 0 } );
				
				planeA.addMatch( new PointMatch( p1, offset ) );
				offset.plane.addMatch( new PointMatch( offset, p1 ) );

				planeA.addConnectedTile( offset.plane );
				offset.plane.addConnectedTile( planeA );
			}
			
			//System.exit( 0 );
		}
		
		final TileConfiguration tc = new TileConfiguration();
		tc.addTiles( planes );

		tc.fixTile( planes.get( 4 ) );
		
		tc.optimize( 10, 1000, 200 );
		
		/*
		double avgError = tc.getError();
		double maxError = tc.getMaxError();			
		System.out.println( avgError  + " " + maxError );
		*/
		
		PrintWriter out = TextFileAccess.openFileWrite( new File( baseDir, "_zPositions.txt" ) );
		
		for ( final MicroscopyPlane plane : planes )
		{
			System.out.println( "Plane " + plane.tileNumber + " of " + plane.name + " <-\t" + plane.getModel().tx );
			out.println( "Plane " + plane.tileNumber + " of " + plane.name + " <-\t" + plane.getModel().tx );
			
			if ( Align.outAllZ != null )
				Align.outAllZ.println( baseDir + "\t" + plane.name + "\t" + plane.tileNumber + "\t" + plane.getModel().tx );
		}

		if ( Align.outAllZ != null )
			Align.outAllZ.flush();

		out.close();
	}
	
	public ArrayList< MicroscopyPlane > getPlanes() { return planes; }
	
	public float computePairwiseAlignment( final String baseDir, final MicroscopyPlane refPlane, final MicroscopyPlane templatePlane ) throws FormatException, IOException
	{
		final String refDir = refPlane.name;
		final int refIndex = refPlane.tileNumber;
		final boolean mirrorRef = refPlane.mirror;
		
		final String templateDir = templatePlane.name;
		final int templateIndex = templatePlane.tileNumber; 
		final boolean mirrorTemplate = templatePlane.mirror;
		
		Image<FloatType> ref, template;
		
		float[] entropiesReference, entropiesTemplate;
		
		final File refFile = new File( baseDir, AlignProperties.tmpName + refDir + "_" + refIndex + AlignProperties.piezoStack );
		final File refEntropiesFile = new File( baseDir, AlignProperties.tmpName + refDir + "_" + refIndex + AlignProperties.entropies );
		final File templateFile = new File( baseDir, AlignProperties.tmpName + templateDir + "_"  + templateIndex + AlignProperties.piezoStack );
		final File templateEntropiesFile = new File( baseDir, AlignProperties.tmpName + templateDir + "_" + templateIndex + AlignProperties.entropies );

		if ( !refEntropiesFile.exists() )
		{
			if ( refFile.exists() )
			{
				ref = LOCI.openLOCIFloatType( refFile.getAbsolutePath(), new ArrayContainerFactory() );
			}
			else
			{
				final Image< FloatType > tmp = piezoStacks.get( refDir );
				if ( tmp == null )
				{
					ref = OpenPiezoStack.openPiezo( new File( baseDir, refDir ).getAbsolutePath() );
					if ( mirrorRef )
						Mirror.horizontal( ref );
					
					piezoStacks.put( refDir, ref.clone() );
				}
				else
				{
					ref = tmp;
				}

				ref = ExtractPlane.extract( ref, refIndex );
				
				// save the extracted stack
				FileSaver fs = new FileSaver( ImageJFunctions.copyToImagePlus( ref ) );
				fs.saveAsTiffStack( refFile.getAbsolutePath() );
			}
			
			if ( AlignProperties.sigma != null )
			{
				final GaussianConvolutionReal< FloatType > gauss = new GaussianConvolutionReal<FloatType>( ref, new OutOfBoundsStrategyMirrorFactory<FloatType>(), AlignProperties.sigma );
				gauss.process();
				ref = gauss.getResult();
			}
			
			//ImageJFunctions.show( ref );
			
			final Image< FloatType > focusStackReference = AutoFocus.focus( ref );
			
			//ImageJFunctions.show( focusStackReference );
			//SimpleMultiThreading.threadHaltUnClean();
			
			entropiesReference = ComputeEntropy.computeEntropyForSlices( focusStackReference, 256 );

			// normalize by avg and stdev
			CrossCorrelation.normalize( entropiesReference );
			
			entropiesReference = CrossCorrelation.median3( entropiesReference );
			
			// save the entropies
			final PrintWriter out = TextFileAccess.openFileWrite( refEntropiesFile );
			out.println( "entries\t" + entropiesReference.length );
			for ( final float e : entropiesReference )
				out.println( e );
			out.close();
		}
		else
		{
			// load entropies
			entropiesReference = loadEntropies( refEntropiesFile );
		}
		
		if ( !templateEntropiesFile.exists() )
		{
			if ( templateFile.exists() )
			{
				template = LOCI.openLOCIFloatType( templateFile.getAbsolutePath(), new ArrayContainerFactory() );
			}
			else
			{
				final Image< FloatType > tmp = piezoStacks.get( templateDir );
				if ( tmp == null )
				{
					template = OpenPiezoStack.openPiezo( new File( baseDir, templateDir ).getAbsolutePath() );
						
					if ( mirrorTemplate )
						Mirror.horizontal( template );
					
					piezoStacks.put( templateDir, template.clone() );
				}
				else
				{
					template = tmp;
				}
				template = ExtractPlane.extract( template, templateIndex );
				
				// save the stack
				FileSaver fs = new FileSaver( ImageJFunctions.copyToImagePlus( template ) );
				fs.saveAsTiffStack( templateFile.getAbsolutePath() );
			}

			if ( AlignProperties.sigma != null )
			{
				final GaussianConvolutionReal< FloatType > gauss = new GaussianConvolutionReal<FloatType>( template, new OutOfBoundsStrategyMirrorFactory<FloatType>(), AlignProperties.sigma );
				gauss.process();
				template = gauss.getResult();
			}
			
			final Image< FloatType > focusStackTemplate = AutoFocus.focus( template );
			entropiesTemplate = ComputeEntropy.computeEntropyForSlices( focusStackTemplate, 256 );
			
			// normalize by avg and stdev
			CrossCorrelation.normalize( entropiesTemplate );

			entropiesTemplate = CrossCorrelation.median3( entropiesTemplate );
			
			// save the entropies
			final PrintWriter out = TextFileAccess.openFileWrite( templateEntropiesFile );
			out.println( "entries\t" + entropiesTemplate.length );
			for ( final float e : entropiesTemplate )
				out.println( e );
			out.close();
		}
		else
		{
			// load entries
			entropiesTemplate = loadEntropies( templateEntropiesFile );
		}		
					
		PrintWriter out = null;//TextFileAccess.openFileWrite( new File( baseDir,"debug_z_registration_" + templateDir + "_" + templateIndex + "-onto-" + refDir + "_" + refIndex + ".txt" ) );
		final float offset = Alignment.align1d( entropiesReference, entropiesTemplate, 1.4, 0.1, out );
		//out.close();
		
		//IJ.log( "offset [px]\t" + offset );
		out = TextFileAccess.appendFileWrite( new File( baseDir,"z_registration.txt" ) );
		out.println( refDir + "_" + refIndex + "\t" + templateDir + "_" + templateIndex + "\t" + offset );
		out.close();
		
		
		final Image< FloatType > alignedTemplate = Alignment.getAlignedSeries( Alignment.createImageFromArray( entropiesTemplate, new int[]{ entropiesTemplate.length } ), offset );
		
		// write a small log file
		out = TextFileAccess.openFileWrite( new File( baseDir,"values_z_registration_" + templateDir + "_" + templateIndex + "-onto-" + refDir + "_" + refIndex + ".txt" ) );
		out.println( "ref_entropy" + "\t" + "template_entropy" + "\t" + "template_adjusted" );
					
		int i = 0;
		for ( final FloatType t : alignedTemplate )
			out.println( entropiesReference[ i ] + "\t" + entropiesTemplate[ i++ ] + "\t" + t );
		out.close();
		
		
		return offset;
	}
	
	protected float[] loadEntropies( final File file ) throws IOException
	{
		final BufferedReader in = TextFileAccess.openFileRead( file );
		
		final int numEntries = Integer.parseInt( in.readLine().split( "\t" )[ 1 ] );
		final float[] entropies = new float[ numEntries ];
		
		int i = 0;
		while ( in.ready() )
			entropies[ i++ ] = Float.parseFloat( in.readLine() );
		
		return entropies;
	}
	
	protected int findMax( final float[] values )
	{
		int max = 0;
		float maxV = values[ 0 ];
		
		for ( int i = 1; i < values.length; ++i )
		{
			if ( values[ i ] > maxV )
			{
				maxV = values[ i ];
				max = i;
			}
		}
		
		return max;
	}
	
	public static void main( String[] args ) throws FormatException, IOException, NotEnoughDataPointsException, IllDefinedDataPointsException
	{
		new ImageJ();
		
		final String[] names = new String[]{ "DNA stack mRNA", "DNA stack NPC" };
		final boolean[] mirror = new boolean[]{ true, false };
		
		new AlignZ( "/Users/preibischs/Documents/Microscopy/david/sample 1/", names, mirror );
	}

}
