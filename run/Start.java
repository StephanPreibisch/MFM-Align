package run;

import ij.IJ;
import ij.ImageJ;
import ij.io.FileSaver;
import io.ExtractPlane;
import io.OpenPiezoStack;
import io.TextFileAccess;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import loci.formats.FormatException;
import mpicbg.imglib.algorithm.gauss.GaussianConvolutionReal;
import mpicbg.imglib.container.array.ArrayContainerFactory;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.interpolation.linear.LinearInterpolatorFactory;
import mpicbg.imglib.io.LOCI;
import mpicbg.imglib.multithreading.SimpleMultiThreading;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyMirrorFactory;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.imglib.util.DevUtil;
import mpicbg.models.RigidModel2D;
import plugin.DescriptorParameters;
import plugin.Descriptor_based_registration;
import process.Alignment;
import process.AutoFocus;
import process.AvgProjection3;
import process.ComputeEntropy;
import process.CrossCorrelation;
import process.Matching;
import process.Mirror;
import process.OverlayFusion;

public class Start 
{
	public static int minNumInliers = 10;
	
	final String tmpName = "tmp_";
	final String piezoProj = "_piezo_avg.tif";
	final String piezoStack = "_piezo.tif";
	
	public Start( final String baseDir, final String refDir, final int refIndex, final String templateDir, final int templateIndex, final boolean mirror ) throws FormatException, IOException
	{
		Image<FloatType> ref, template, refProjection, templateProjection;
		
		final File refFile = new File( baseDir, tmpName + refDir + "_" + refIndex + piezoStack );
		final File refProjFile = new File( baseDir, tmpName + refDir + "_" + refIndex + piezoProj );
		
		final File templateFile = new File( baseDir, tmpName + templateDir + "_" + templateIndex + "-onto-" +  refDir + "_" + refIndex + piezoStack );
		final File templateProjFile = new File( baseDir, tmpName + templateDir + "_" + templateIndex + "-onto-" +  refDir + "_" + refIndex + piezoProj );
		
		if ( refFile.exists() && templateFile.exists() && refProjFile.exists() && templateProjFile.exists() )
		{
			ref = LOCI.openLOCIFloatType( refFile.getAbsolutePath(), new ArrayContainerFactory() );
			template = LOCI.openLOCIFloatType( templateFile.getAbsolutePath(), new ArrayContainerFactory() );
			refProjection = LOCI.openLOCIFloatType( refProjFile.getAbsolutePath(), new ArrayContainerFactory() );
			templateProjection = LOCI.openLOCIFloatType( templateProjFile.getAbsolutePath(), new ArrayContainerFactory() );
		}
		else
		{	
			// load the files and extract the center planes
			if ( refFile.exists() )
				ref = LOCI.openLOCIFloatType( refFile.getAbsolutePath(), new ArrayContainerFactory() );
			else
				ref = ExtractPlane.extract( OpenPiezoStack.openPiezo( new File( baseDir, refDir ).getAbsolutePath() ), refIndex );
			
			template = ExtractPlane.extract( OpenPiezoStack.openPiezo( new File( baseDir, templateDir ).getAbsolutePath() ), templateIndex );
	
			// mirror the npc
			if ( mirror )
				Mirror.horizontal( template );
			
			// do the average intensity projections
			if ( refFile.exists() )
				refProjection = LOCI.openLOCIFloatType( refProjFile.getAbsolutePath(), new ArrayContainerFactory() );
			else
				refProjection = AvgProjection3.project( ref );
			
			templateProjection = AvgProjection3.project( template );
	
			// compute the per-plane registration
			// of NPC and mRNA
			final DescriptorParameters params = getParametersForProjection();
			final int numInliers = Matching.descriptorBasedRegistration( ImageJFunctions.copyToImagePlus( templateProjection ), ImageJFunctions.copyToImagePlus( refProjection ), params );
			
			if ( numInliers < minNumInliers )
			{
				IJ.log( "num inliers: " + numInliers );
				IJ.log( "insufficient inliers (at least " + minNumInliers + " required), please check manually what is happening." );
				ImageJFunctions.show( refProjection ).setTitle( "reference = " + refDir + "_" + refIndex + piezoStack );
				ImageJFunctions.show( templateProjection ).setTitle( "template = " + templateDir + "_" + templateIndex + piezoProj );
				
				SimpleMultiThreading.threadHaltUnClean();
				
				return;
			}	

			// compute cross-correlation before
			final float r1 = CrossCorrelation.corrlate( ref, template );

			// re-arrange the mRNA channel
			OverlayFusion.fuseChannel( templateProjection, templateProjection.clone(), new float[ templateProjection.getNumDimensions() ], Descriptor_based_registration.lastModel1, new LinearInterpolatorFactory<FloatType>( new OutOfBoundsStrategyMirrorFactory<FloatType>() ) );
			OverlayFusion.fuseChannel( template, template.clone(), new float[ template.getNumDimensions() ], Descriptor_based_registration.lastModel1, new LinearInterpolatorFactory<FloatType>( new OutOfBoundsStrategyMirrorFactory<FloatType>() ) );

			// compute cross-correlation
			final float r2 = CrossCorrelation.corrlate( ref, template );

			// write some debugging info
			IJ.log( "num inliers:\t" + numInliers );
			IJ.log( "cross correlation before alignment:\t" + r1 );
			IJ.log( "cross correlation after alignment:\t" + r2 );
			IJ.log( "model mapping template (" + templateDir + "_" + templateIndex + ") onto reference (" + refDir + "_" + refIndex + "):\t" + Descriptor_based_registration.lastModel1 );

			// write a small log file
			PrintWriter out = TextFileAccess.openFileWrite( new File( baseDir,"log_image_registration_" + templateDir + "_" + templateIndex + "-onto-" + refDir + "_" + refIndex + ".txt" ) );
			
			out.println( "num inliers:\t" + numInliers );
			out.println( "cross correlation before alignment:\t" + r1 );
			out.println( "cross correlation after alignment:\t" + r2 );
			out.println( "model mapping template (" + templateDir + "_" + templateIndex + ") onto reference (" + refDir + "_" + refIndex + "):\t" + Descriptor_based_registration.lastModel1 );
			
			out.close();

			// save the projections
			FileSaver fs = new FileSaver( ImageJFunctions.copyToImagePlus( refProjection ) );
			fs.saveAsTiff( refProjFile.getAbsolutePath() );
	
			fs = new FileSaver( ImageJFunctions.copyToImagePlus( templateProjection ) );
			fs.saveAsTiff( templateProjFile.getAbsolutePath() );
			
			// save the extracted stacks
			fs = new FileSaver( ImageJFunctions.copyToImagePlus( ref ) );
			fs.saveAsTiffStack( refFile.getAbsolutePath() );
	
			fs = new FileSaver( ImageJFunctions.copyToImagePlus( template ) );
			fs.saveAsTiffStack( templateFile.getAbsolutePath() );
		}		
		
		GaussianConvolutionReal< FloatType > gauss = new GaussianConvolutionReal<FloatType>( ref, new OutOfBoundsStrategyMirrorFactory<FloatType>(), new double[]{0.75, 0.75, 4} );
		gauss.process();
		ref = gauss.getResult();
		
		gauss = new GaussianConvolutionReal<FloatType>( template, new OutOfBoundsStrategyMirrorFactory<FloatType>(), new double[]{0.75, 0.75, 4} );
		gauss.process();
		template = gauss.getResult();
		
		
		final Image< FloatType > focusStackReference = AutoFocus.focus( ref );
		final Image< FloatType > focusStackTemplate = AutoFocus.focus( template );
		
		final float[] entropiesReference = ComputeEntropy.computeEntropyForSlices( focusStackReference, 256 );
		final float[] entropiesTemplate = ComputeEntropy.computeEntropyForSlices( focusStackTemplate, 256 );

		//ImageJFunctions.show( focusStackReference );
		//ImageJFunctions.show( focusStackTemplate );
		
		// normalize by avg and stdev
		CrossCorrelation.normalize( entropiesReference );
		CrossCorrelation.normalize( entropiesTemplate );
				
		PrintWriter out = TextFileAccess.openFileWrite( new File( baseDir,"debug_z_registration_" + templateDir + "_" + templateIndex + "-onto-" + refDir + "_" + refIndex + ".txt" ) );
		final float offset = Alignment.align1d( entropiesReference, entropiesTemplate, 1.4, 0.1, out );
		out.close();
		
		out = TextFileAccess.appendFileWrite( new File( baseDir,"z_registration.txt" ) );
		IJ.log( "offset [px]\t" + offset );
		out.println( refDir + "_" + refIndex + "\t" + templateDir + "_" + templateIndex + "\t" + offset );
		out.close();
		/*
		final Image< FloatType > alignedTemplate = Alignment.getAlignedSeries( DevUtil.createImageFromArray( entropiesTemplate, new int[]{ entropiesTemplate.length } ), offset );
		
		// write a small log file
		out = TextFileAccess.openFileWrite( new File( baseDir,"values_z_registration_" + templateDir + "_" + templateIndex + "-onto-" + refDir + "_" + refIndex + ".txt" ) );
		out.println( "ref_entropy" + "\t" + "template_entropy" + "\t" + "template_adjusted" );
					
		int i = 0;
		for ( final FloatType t : alignedTemplate )
			out.println( entropiesReference[ i ] + "\t" + entropiesTemplate[ i++ ] + "\t" + t );
		out.close();
	*/
		

		/*
		final int planesStart = 85;
		final int numPlanes = 80;
		final int range = 10;
		final float[][] r = new float[ numPlanes ][];
		
		for ( int i = 0; i < numPlanes; ++i )
			r[ i ] = CrossCorrelation.corrlatePlane( mRNA, npc, planesStart + i, range );
		
		for ( int i = 0; i < numPlanes; ++i )
		{
			int max = findMax( r[ i ] );
			
			System.out.println( (planesStart+i) + "\t"  + (max-range) );
		}
		
		
		
		for ( int i = 0; i < 2 * range + 1; ++i )
		{
			System.out.print(  (i-range) );
			for ( int j = 0; j < numPlanes; ++j )
				System.out.print( "\t" + r[ j ][ i ] );
			System.out.println();
		}
		*/
	}
	
	public static int findMax( final float[] values )
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
	
	protected DescriptorParameters getParametersForProjection()
	{
		final DescriptorParameters params = new DescriptorParameters();
		
		params.dimensionality = 2;
		params.sigma1 = 2.099f;
		params.sigma2 = 2.4961457f;
		params.threshold = 0.020566484f;
		params.lookForMaxima = true;
		params.lookForMinima = true;
		params.model = new RigidModel2D();
		params.similarOrientation = true;
		params.numNeighbors = 3;
		params.redundancy = 1;
		params.significance = 3;
		params.ransacThreshold = 5;
		params.channel1 = 0;
		params.channel2 = 0;
		
		// for stack-registration
		params.globalOpt = 0; // 0=all-to-all; 1=all-to-all-withrange; 2=all-to-1; 3=Consecutive
		params.range = 5;	
		params.directory = "";
		
		params.reApply = false;
		params.roi1 = null;
		params.roi2 = null;
		
		params.setPointsRois = false;
		
		params.silent = true;
		
		// 0 == fuse in memory, 1 == write to disk, 2 == nothing
		params.fuse = 2;
		
		return params;
	}
	
	public static void main( String[] args ) throws FormatException, IOException
	{
		new ImageJ();
		final int reference = 4;
		
		for ( int i = 0; i < 9; ++i )
		{
			if ( i != reference )
				new Start( "/Users/preibischs/Documents/Microscopy/david/sample 1/", "DNA stack mRNA", 4, "DNA stack mRNA", i, false );
		}

		for ( int i = 0; i < 9; ++i )
			new Start( "/Users/preibischs/Documents/Microscopy/david/sample 1/", "DNA stack mRNA", 4, "DNA stack NPC", i, true );

	}

}
