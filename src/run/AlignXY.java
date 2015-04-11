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
import io.TextFileAccess;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import loci.formats.FormatException;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.models.AbstractModel;
import mpicbg.models.IllDefinedDataPointsException;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.RigidModel2D;
import plugin.DescriptorParameters;
import plugin.Descriptor_based_series_registration;
import process.AvgProjection3;
import process.Matching;
import process.OverlayFusion;
import run.MicroscopyPlane.Mirroring;

/**
 * Class that runs the XY alignment using Descriptor-based registration package.
 * 
 * @author preibischs
 *
 */
public class AlignXY 
{
	final ArrayList< MicroscopyPlane > planes;
	
	/**
	 * Makes the max projections for all {@link MicroscopyPlane}s, alignes them and stores the transformation model in the
	 * {@link MicroscopyPlane} object.
	 * 
	 * @param planes - the {@link MicroscopyPlane}s
	 * @param model - which model to use, e.g. RigidModel2d, AffineModel2d
	 * @throws FormatException
	 * @throws IOException
	 */
	public AlignXY( final ArrayList< MicroscopyPlane > planes, final AbstractModel< ? > model ) throws FormatException, IOException
	{
		this.planes = planes;
		
		align( model );
	}
	
	public void align( final AbstractModel< ? > model ) throws FormatException, IOException
	{
		//
		// extract a stack of z-avg-projections
		//
		ImageStack planeStack = null;
		for ( final MicroscopyPlane plane : planes )
		{
			final Image< FloatType > image = plane.getImagePiezo();
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

		//stack.getProcessor().resetMinAndMax();
		//stack.show();
		//SimpleMultiThreading.threadHaltUnClean();
		
		// compute the per-plane registration
		// of NPC and mRNA
		final DescriptorParameters params = getParametersForProjection( model );
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

	protected DescriptorParameters getParametersForProjection( final AbstractModel< ? > model )
	{
		final DescriptorParameters params = new DescriptorParameters();
		
		params.dimensionality = 2;
		params.sigma1 = 2.99f; // before 2.099f
		params.sigma2 = 3.55f; // before 2.4961457f;
		params.threshold = 0.010566484f;
		params.lookForMaxima = true;
		params.lookForMinima = true;
		
		if ( model == null )
			params.model = new RigidModel2D();
		else
			params.model = model;
			
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
		
		final String root = "/media/f3df52e9-9b20-4b66-a4cc-1fc741ff481e/stephan/";
		final String experimentDir = "1 (20110525, dish 2, cell 22)";
		
		final String localDir = "DNA stack";
		
		final String[] tags = new String[] { "2464"/*"green"*/, "4283"/*"red"*/ };
		final Mirroring[] mirror = new Mirroring[]{ Mirroring.HORIZONTALLY, Mirroring.DONOT };
		
		final String darkCounts[] = new String[]{ root + "Dark Counts/MED_avgstack_DNA_2464 green.tif",
												   root + "Dark Counts/MED_avgstack_DNA_4283 red.tif" }; // can be null

		//
		// set up the planes
		// 	
		final ArrayList< MicroscopyPlane > planes = new ArrayList< MicroscopyPlane >();
		
		for ( int c = 0; c < tags.length; ++c )
			for ( int t = 0; t < AlignProperties.numTiles; ++t )
				planes.add( new MicroscopyPlane( root + experimentDir, localDir, tags[ c ], darkCounts[ c ], mirror[ c ], t ) );

		new AlignXY( planes, new RigidModel2D() );
	}
	
}
