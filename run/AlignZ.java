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
import run.MicroscopyPlane.Mirroring;

public class AlignZ 
{
	final ArrayList< MicroscopyPlane > planes;
	final HashMap< String, Image<FloatType> > allPiezoStacks = new HashMap<String, Image<FloatType>>();	
	
	public AlignZ( final ArrayList< MicroscopyPlane > planes ) throws FormatException, IOException, NotEnoughDataPointsException, IllDefinedDataPointsException
	//final String baseDir, final String[] names, final boolean[] mirror ) throws FormatException, IOException, NotEnoughDataPointsException, IllDefinedDataPointsException
	{
		//
		// set up the planes
		// 	
		this.planes = planes;
		
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
					final float offset = computePairwiseAlignment( planeA, planeB );
					
					// if it is from the same channel
					if ( planeB.getTagName().equals( planeA.getTagName() ) )
						planeA.offsetsSameChannel.add( new PlaneOffset( planeB, offset ) );
					else
						planeA.offsetsOtherChannel.add( new PlaneOffset( planeB, offset ) );					
				}
			}
			
			// compute RANSAC to filter the correct ones
		
			//for ( final PlaneOffset planeOffset : planeA.offsets1 )
			//	System.out.println( "same " + planeA.name + " " + planeA.tileNumber + " <-> " + planeOffset.plane.name + " " + planeOffset.plane.tileNumber + " = " + planeOffset.offset );

			int numRemoved = MicroscopyPlane.removeOutliers( planeA.offsetsSameChannel, AlignProperties.epsilon, AlignProperties.minInlierRatio );
			System.out.println( "same " + planeA.getFullName() + " " + planeA.tileNumber + " -> removed " + numRemoved );
			
			//for ( final PlaneOffset planeOffset : planeA.offsets2 )
			//	System.out.println( "diff " + planeA.name + " " + planeA.tileNumber + " <-> " + planeOffset.plane.name + " " + planeOffset.plane.tileNumber + " = " + planeOffset.offset );

			numRemoved = MicroscopyPlane.removeOutliers( planeA.offsetsOtherChannel, AlignProperties.epsilon, AlignProperties.minInlierRatio );
			System.out.println( "diff " + planeA.getFullName() + " " + planeA.tileNumber + " -> removed " + numRemoved );
			
			//for ( final PlaneOffset planeOffset : planeA.offsets2 )
			//	System.out.println( "diff " + planeA.name + " " + planeA.tileNumber + " <-> " + planeOffset.plane.name + " " + planeOffset.plane.tileNumber + " = " + planeOffset.offset );

			final ArrayList< PlaneOffset > allPlanes = new ArrayList<PlaneOffset>();
			allPlanes.addAll( planeA.offsetsSameChannel );
			allPlanes.addAll( planeA.offsetsOtherChannel );
			
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
		
		PrintWriter out = TextFileAccess.openFileWrite( new File( planes.get( 0 ).getBaseDirectory(), "_zPositions.txt" ) );
		
		for ( final MicroscopyPlane plane : planes )
		{
			System.out.println( "Plane " + plane.getFullName() + " <-\t" + plane.getModel().tx );
			out.println( "Plane " + plane.getFullName() + " <-\t" + plane.getModel().tx );
			
			if ( Align.outAllZ != null )
				Align.outAllZ.println( planes.get( 0 ).getBaseDirectory() + "\t" + plane.getFullName() + "\t" + plane.tileNumber + "\t" + plane.getModel().tx );
		}

		if ( Align.outAllZ != null )
			Align.outAllZ.flush();

		out.close();
	}
	
	public ArrayList< MicroscopyPlane > getPlanes() { return planes; }
	
	public float computePairwiseAlignment( final MicroscopyPlane refPlane, final MicroscopyPlane templatePlane ) throws FormatException, IOException
	{
		Image<FloatType> ref, template;
		
		float[] entropiesReference, entropiesTemplate;
		
		final File refFile = new File( refPlane.getBaseDirectory(), AlignProperties.tmpName + refPlane.getFullName() + AlignProperties.piezoStack );
		final File refEntropiesFile = new File( refPlane.getBaseDirectory(), AlignProperties.tmpName + refPlane.getFullName() + AlignProperties.entropies );
		final File templateFile = new File( templatePlane.getBaseDirectory(), AlignProperties.tmpName + templatePlane.getFullName() + AlignProperties.piezoStack );
		final File templateEntropiesFile = new File( templatePlane.getBaseDirectory(), AlignProperties.tmpName + templatePlane.getFullName() + AlignProperties.entropies );

		// try to load the entropies
		if ( !refEntropiesFile.exists() )
		{
			// try to load the raw dna stack per plane 
			if ( refFile.exists() )
			{
				ref = LOCI.openLOCIFloatType( refFile.getAbsolutePath(), new ArrayContainerFactory() );
			}
			else
			{
				// do not open the huge image (containing 9 planes) for every new tile that we extract
				if ( ( ref = allPiezoStacks.get( refPlane.getTagName() ) ) == null )
				{
					ref = OpenPiezoStack.openPiezo( new File( refPlane.getBaseDirectory(), refPlane.getLocalDirectory() ), refPlane.getTagName() );

					if ( refPlane.getMirror() == Mirroring.HORIZONTALLY )
						Mirror.horizontal( ref );

					allPiezoStacks.put( refPlane.getTagName(), ref.clone() );
				}

				ref = ExtractPlane.extract( ref, refPlane.getTileNumber() );
								
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
				if ( ( template = allPiezoStacks.get( templatePlane.getTagName() ) ) == null )
				{
					template = OpenPiezoStack.openPiezo( new File( templatePlane.getBaseDirectory(), templatePlane.getLocalDirectory() ), templatePlane.getTagName() );
						
					if ( templatePlane.getMirror() == Mirroring.HORIZONTALLY )
						Mirror.horizontal( template );
					
					allPiezoStacks.put( templatePlane.getTagName(), template.clone() );
				}

				template = ExtractPlane.extract( template, templatePlane.getTileNumber() );
				
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
		out = TextFileAccess.appendFileWrite( new File( refPlane.getBaseDirectory() ,"z_registration.txt" ) );
		out.println( refPlane.getFullName() + "\t" + templatePlane.getFullName() + "\t" + offset );
		out.close();
		
		
		final Image< FloatType > alignedTemplate = Alignment.getAlignedSeries( Alignment.createImageFromArray( entropiesTemplate, new int[]{ entropiesTemplate.length } ), offset );
		
		// write a small log file
		out = TextFileAccess.openFileWrite( new File( refPlane.getBaseDirectory(), "values_z_registration_" + templatePlane.getFullName() + "-onto-" + refPlane.getFullName() + ".txt" ) );
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
		
		final String root = "/home/stephanpreibisch/Desktop/stephan/";
		final String experimentDir = "1 (20110525, dish 2, cell 22)";
		
		final String localDir = "DNA stack";
		
		final String[] tags = new String[] { "green", "red" };
		final Mirroring[] mirror = new Mirroring[]{ Mirroring.HORIZONTALLY, Mirroring.DONOT };
		
		//
		// set up the planes
		// 	
		final ArrayList< MicroscopyPlane > planes = new ArrayList< MicroscopyPlane >();
		
		for ( int c = 0; c < tags.length; ++c )
			for ( int t = 0; t < AlignProperties.numTiles; ++t )
				planes.add( new MicroscopyPlane( root + experimentDir, localDir, tags[ c ], mirror[ c ], t ) );
		
		new AlignZ( planes );
	}

}
