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

import loci.formats.FormatException;
import mpicbg.imglib.container.array.Array;
import mpicbg.imglib.container.array.ArrayContainerFactory;
import mpicbg.imglib.container.basictypecontainer.array.FloatArray;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.ImageFactory;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.interpolation.linear.LinearInterpolatorFactory;
import mpicbg.imglib.io.LOCI;
import mpicbg.imglib.multithreading.SimpleMultiThreading;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyMirrorFactory;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.models.IllDefinedDataPointsException;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.RigidModel2D;
import plugin.DescriptorParameters;
import plugin.Descriptor_based_series_registration;
import process.AvgProjection3;
import process.CrossCorrelation;
import process.Matching;
import process.Mirror;
import process.OverlayFusion;

public class AlignXY 
{
	final ArrayList< MicroscopyPlane > planes;
	
	public AlignXY( final ArrayList< MicroscopyPlane > planes ) throws FormatException, IOException
	{
		this.planes = planes;
		
		align();
	}
	
	public AlignXY( final String baseDir, final String[] names, final boolean[] mirror ) throws FormatException, IOException
	{
		//
		// set up the planes
		// 	
		planes = new ArrayList< MicroscopyPlane >();
		
		for ( int c = 0; c < names.length; ++c )
			for ( int t = 0; t < AlignProperties.numTiles; ++t )
				planes.add( new MicroscopyPlane( baseDir, names[ c ], mirror[ c ], t ) );
		
		align();
	}
	
	public void align() throws FormatException, IOException
	{
		//
		// extract a stack of z-avg-projections
		//
		ImageStack planeStack = null;
		for ( final MicroscopyPlane plane : planes )
		{
			final Image< FloatType > image = getImage( plane.baseDir, plane.name, plane.tileNumber, plane.mirror );
			final Image< FloatType > planeImg = AvgProjection3.project( image );
			
			// store the images
			plane.setImage( image );
			plane.setAvgProj( planeImg );
			
			if ( planeStack == null )
				planeStack = new ImageStack( planeImg.getDimension( 0 ), planeImg.getDimension( 1 ) );
			planeStack.addSlice( plane.name + "_" + plane.tileNumber, ImageJFunctions.copyToImagePlus( planeImg ).getProcessor() );
		}
		
		ImagePlus stack = new ImagePlus( "stack of avg proj", planeStack );
		
		// make it a timelapse and not a stack
		stack = OverlayFusion.switchZTinXYCZT( stack );

		//stack.show();
		//SimpleMultiThreading.threadHaltUnClean();
		
		// compute the per-plane registration
		// of NPC and mRNA
		final DescriptorParameters params = getParametersForProjection();
		Matching.descriptorBasedStackRegistration( stack, params );

		// set the models
		PrintWriter out = TextFileAccess.openFileWrite( new File( planes.get( 0 ).baseDir, "_xyAlignment.txt" ) );
		for ( int i = 0; i < planes.size(); ++i )
		{
			planes.get( i ).setXYModel( Descriptor_based_series_registration.lastModels.get( i ) );
			out.println( "Plane " + planes.get( i ).tileNumber + " of " + planes.get( i ).name + " <-\t" + Descriptor_based_series_registration.lastModels.get( i ) );
			
			if ( Align.outAllXY != null )
				Align.outAllXY.println( planes.get( i ).baseDir + "\t" + planes.get( i ).name + "\t" + planes.get( i ).tileNumber + "\t" + Descriptor_based_series_registration.lastModels.get( i ) );

		}
		if ( Align.outAllXY != null )
			Align.outAllXY.flush();
		
		out.close();
	}

	
	public Image<FloatType> getImage( final String baseDir, final String dir, final int index, final boolean mirror ) throws FormatException, IOException
	{	
		Image<FloatType> image;
		
		final File file = new File( baseDir, AlignProperties.tmpName + dir + "_" + index + AlignProperties.piezoStack );
		
		// load or create the 3d-stack
		if ( file.exists() )
		{
			image = LOCI.openLOCIFloatType( file.getAbsolutePath(), new ArrayContainerFactory() );
		}
		else
		{
			image = OpenPiezoStack.openPiezo( new File( baseDir, dir ).getAbsolutePath() );
			if ( mirror )
				Mirror.horizontal( image );
			image = ExtractPlane.extract( image, index );
			
			// save the extracted stack
			FileSaver fs = new FileSaver( ImageJFunctions.copyToImagePlus( image ) );
			fs.saveAsTiffStack( file.getAbsolutePath() );
		}
		
		// load or create the average-projection
		return image;
	}

	protected DescriptorParameters getParametersForProjection()
	{
		final DescriptorParameters params = new DescriptorParameters();
		
		params.dimensionality = 2;
		params.sigma1 = 2.99f;
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
		
		params.setPointsRois = true;
		
		params.silent = false;
		
		// 0 == fuse in memory, 1 == write to disk, 2 == nothing
		params.fuse = 2;
		
		return params;
	}

	public static void main( String[] args ) throws FormatException, IOException, NotEnoughDataPointsException, IllDefinedDataPointsException
	{
		new ImageJ();
		
		final String[] names = new String[]{ "DNA stack mRNA", "DNA stack NPC" };
		final boolean[] mirror = new boolean[]{ true, false };
		
		new AlignXY( "/Users/preibischs/Documents/Microscopy/david/sample 1/", names, mirror );
	}
	
}
