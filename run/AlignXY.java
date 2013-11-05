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
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.FileSaver;
import io.ExtractPlane;
import io.OpenPiezoStack;
import io.TextFileAccess;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import loci.formats.FormatException;
import mpicbg.imglib.container.array.ArrayContainerFactory;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.io.LOCI;
import mpicbg.imglib.multithreading.SimpleMultiThreading;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.models.IllDefinedDataPointsException;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.RigidModel2D;
import plugin.DescriptorParameters;
import plugin.Descriptor_based_series_registration;
import process.AvgProjection3;
import process.Matching;
import process.Mirror;
import process.OverlayFusion;
import run.MicroscopyPlane.Mirroring;

public class AlignXY 
{
	final ArrayList< MicroscopyPlane > planes;
	
	public AlignXY( final ArrayList< MicroscopyPlane > planes ) throws FormatException, IOException
	{
		this.planes = planes;
		
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
			final Image< FloatType > image = getImagePiezo( plane );
			final Image< FloatType > planeImg = AvgProjection3.project( image );
			
			// store the images
			plane.setImage( image );
			plane.setAvgProj( planeImg );
			
			if ( planeStack == null )
				planeStack = new ImageStack( planeImg.getDimension( 0 ), planeImg.getDimension( 1 ) );
			planeStack.addSlice( plane.getFullName(), ImageJFunctions.copyToImagePlus( planeImg ).getProcessor() );
		}
		
		ImagePlus stack = new ImagePlus( "stack of avg proj", planeStack );
		
		// make it a timelapse and not a stack
		stack = OverlayFusion.switchZTinXYCZT( stack );

		stack.getProcessor().resetMinAndMax();
		stack.show();
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
			out.println( "Plane " + planes.get( i ).getFullName() + " <-\t" + Descriptor_based_series_registration.lastModels.get( i ) );
			
			if ( Align.outAllXY != null )
				Align.outAllXY.println( planes.get( i ).baseDir + "\t" + planes.get( i ).getFullName() + "\t" + Descriptor_based_series_registration.lastModels.get( i ) );

		}
		if ( Align.outAllXY != null )
			Align.outAllXY.flush();
		
		out.close();
	}

	/**
	 * It looks if this file has been already converted from the slices in the directory into its own file, if so it just loads it
	 * 
	 * @param baseDir - base directory of the MFM experiment
	 * @param dir - the directory that contains all the individual slices of the piezo DNA stack
	 * @param fileNameTag - if multiple channels are in the directory a String that selects for the current channel (e.g. green)
	 * @param index - which of the tiles to load (0...8)
	 * @param mirror - mirror the image or not
	 * @return
	 * @throws FormatException
	 * @throws IOException
	 */
	//public Image<FloatType> getImagePiezo( final String baseDir, final String dir, final String fileNameTag, final int index, final Mirroring mirror ) throws FormatException, IOException
	public Image<FloatType> getImagePiezo( final MicroscopyPlane plane ) throws FormatException, IOException
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

	protected DescriptorParameters getParametersForProjection()
	{
		final DescriptorParameters params = new DescriptorParameters();
		
		params.dimensionality = 2;
		params.sigma1 = 2.99f; // before 2.099f
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

		new AlignXY( planes );
	}
	
}
