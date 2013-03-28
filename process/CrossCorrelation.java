package process;

import java.util.Arrays;

import mpicbg.imglib.container.array.Array;
import mpicbg.imglib.container.basictypecontainer.array.FloatArray;
import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.cursor.LocalizableByDimCursor;
import mpicbg.imglib.cursor.LocalizableCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.interpolation.Interpolator;
import mpicbg.imglib.multithreading.SimpleMultiThreading;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.imglib.util.Util;

public class CrossCorrelation 
{
	public static float[] corrlatePlane( final Image< FloatType > image1, final Image< FloatType > image2, final int plane1, final int range )
	{
		final float[] r = new float[ 2* range + 1 ];
		
		final int[] size = new int[ 2 ];
		size[ 0 ] = image1.getDimension( 0 );
		size[ 1 ] = image1.getDimension( 1 );
		
		final Image< FloatType > planeA = image1.createNewImage( size );
		final Image< FloatType > planeB = image1.createNewImage( size );
			
		extractPlane( image1, planeA, plane1 );
		
		//planeA.getDisplay().setMinMax();
		//ImageJFunctions.copyToImagePlus(planeA ).show();
		
		int j = 0;
		
		for ( int i = plane1 - range; i <= plane1 + range; i = i + 1 )
		{
			extractPlane( image2, planeB, i );
			
			QuantileNormalization.normalizeTo( planeA, planeB );
			
/*
			planeB.getDisplay().setMinMax();
			ImageJFunctions.copyToImagePlus(planeB ).show();
			
			QuantileNormalization.normalizeTo( planeA, planeB );

			planeB.getDisplay().setMinMax();
			ImageJFunctions.copyToImagePlus(planeB ).show();
			
			SimpleMultiThreading.threadHaltUnClean();
	*/		
			r[ j++ ] = corrlate( planeA, planeB );
			
			//System.out.println( (plane1-i) + "\t" + corrlate( planeA, planeB ) );// + "\t" + difference( planeA, planeB ) + "\t" + squareDifference( planeA, planeB ) + "\t" + avg( planeA ) + "\t" + avg( planeB ) );
		}
		
		return r;
	}
	
	public static void extractPlane( final Image< FloatType > image, final Image< FloatType > plane, final int planeNo )
	{
		final int pos[] = new int[ 3 ];
		
		final LocalizableCursor< FloatType > cursor = plane.createLocalizableCursor();
		final LocalizableByDimCursor< FloatType > randomAccess = image.createLocalizableByDimCursor();
		
		while ( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.getPosition( pos );
			pos[ 2 ] = planeNo;
			randomAccess.setPosition( pos );
			
			cursor.getType().set( randomAccess.getType() );
		}
	}

	public static void extractPlane( final Interpolator< FloatType > interpolator, final Image< FloatType > plane, final float planeNo )
	{
		final int[] tmp = new int[ 2 ];
		final float pos[] = new float[ 3 ];
		
		final LocalizableCursor< FloatType > cursor = plane.createLocalizableCursor();
		
		while ( cursor.hasNext() )
		{
			cursor.fwd();
			cursor.getPosition( tmp );
			pos[ 0 ] = tmp[ 0 ];
			pos[ 1 ] = tmp[ 1 ];
			pos[ 2 ] = planeNo;
			interpolator.setPosition( pos );
			
			cursor.getType().set( interpolator.getType() );
		}
	}

	public static float difference( final Image< FloatType > image1, final Image< FloatType > image2 )
	{
		final int numPx = image1.getNumPixels();
		
		final Cursor< FloatType > cursor1 = image1.createCursor();
		final Cursor< FloatType > cursor2 = image2.createCursor();
		
		double diff = 0;
		
		while ( cursor1.hasNext() )
			diff += Math.abs( cursor1.next().getRealFloat() - cursor2.next().getRealFloat() );

		return (float)(diff/numPx);
	}

	public static float avg( final Image< FloatType > image1 )
	{
		final int numPx = image1.getNumPixels();
		
		final Cursor< FloatType > cursor1 = image1.createCursor();
		
		//
		// compute average
		//
		double avg1 = 0;
		
		
		while ( cursor1.hasNext() )
		{
			avg1 += cursor1.next().getRealFloat();
			
		}

		avg1 /= (double) numPx;
		
		return (float)avg1;
	}

	public static float avg10( final Image< FloatType > image1 )
	{
		final int numPx = image1.getNumPixels();
		
		final float[] sort = new float[ numPx ];
		
		int i = 0;
		for ( final FloatType t : image1 )
			sort[ i++ ] = t.get();
		
		Arrays.sort( sort );
		
		final int inside = numPx / 10;
		
		final float[] highest = new float[ inside ];
		
		for ( int j = 0; j < inside; ++j )
			highest[ j ] = sort[ sort.length - 1 - j ];
		
		//
		// compute average
		//
		double avg1 = 0;
		
		for ( final float v: highest )
			avg1 += v;
		
		avg1 /= (double) inside;
		
		return (float)avg1;
	}

	public static float squareDifference( final Image< FloatType > image1, final Image< FloatType > image2 )
	{
		final int numPx = image1.getNumPixels();
		
		final Cursor< FloatType > cursor1 = image1.createCursor();
		final Cursor< FloatType > cursor2 = image2.createCursor();
		
		double diff = 0;
		
		while ( cursor1.hasNext() )
			diff += Math.pow( cursor1.next().getRealFloat() - cursor2.next().getRealFloat(), 2 );

		return (float)(Math.sqrt( diff )/numPx);
	}

	public static float[] median3( final float[] image1 )
	{
		float[] filtered = new float[ image1.length ];
		float[] tmp = new float[ 3 ];
		
		for ( int i = 1; i < image1.length - 1; ++i )
		{
			tmp[ 0 ] = image1[ i - 1 ];
			tmp[ 1 ] = image1[ i ];
			tmp[ 2 ] = image1[ i + 1 ];
			
			Arrays.sort( tmp );
			
			filtered[ i ] = tmp[ 1 ];
		}
		
		filtered[ 0 ] = Math.min( image1[ 0 ], image1[ 1 ] );
		filtered[ image1.length - 1 ] = Math.min( image1[ image1.length - 1 ], image1[ image1.length - 2 ] );
		
		return filtered;
	}
	public static void normalize( final float[] image1 )
	{
		final int numPx = image1.length;
		
		//
		// compute average
		//
		double avg1 = 0;

		for ( int i = 0; i < numPx; ++i )
			avg1 += image1[ i ];

		avg1 /= (double) numPx;
				
		//
		// compute cross correlation
		//
		double var1 = 0;
		
		for ( int i = 0; i < numPx; ++i )
		{
			final float pixel1 = image1[ i ];
			final double dist1 = pixel1 - avg1;
			var1 += dist1 * dist1;
		}		
		
		var1 /= (double) numPx;
		double stDev1 = Math.sqrt(var1);
		
		for ( int i = 0; i < numPx; ++i )
			image1[ i ] = (float)((image1[ i ] - avg1)/stDev1);
	}

	public static float corrlate( final float[] image1, final float[] image2 )
	{
		final int numPx = image1.length;
		
		//
		// compute average
		//
		double avg1 = 0;
		double avg2 = 0;

		for ( int i = 0; i < numPx; ++i )
		{			
			avg1 += image1[ i ];
			avg2 += image2[ i ];
		}

		avg1 /= (double) numPx;
		avg2 /= (double) numPx;
				
		//
		// compute cross correlation
		//
		double var1 = 0, var2 = 0;
		double coVar = 0;
		
		for ( int i = 0; i < numPx; ++i )
		{
			final float pixel1 = image1[ i ];
			final float pixel2 = image2[ i ];
			
			final double dist1 = pixel1 - avg1;
			final double dist2 = pixel2 - avg2;

			coVar += dist1 * dist2;
			var1 += dist1 * dist1;
			var2 += dist2 * dist2;
		}		
		
		var1 /= (double) numPx;
		var2 /= (double) numPx;
		coVar /= (double) numPx;

		double stDev1 = Math.sqrt(var1);
		double stDev2 = Math.sqrt(var2);

		// all pixels had the same color....
		if (stDev1 == 0 || stDev2 == 0)
		{
			if ( stDev1 == stDev2 && avg1 == avg2 )
				return 1;
			else
				return 0;
		}

		// compute correlation coeffienct
		return (float)(coVar / (stDev1 * stDev2));	
	}

	public static float corrlate( final Image< FloatType > image1, final Image< FloatType > image2 )
	{
		final int numPx = image1.getNumPixels();
		
		final Cursor< FloatType > cursor1 = image1.createCursor();
		final Cursor< FloatType > cursor2 = image2.createCursor();
		
		//
		// compute average
		//
		double avg1 = 0;
		double avg2 = 0;
		
		while ( cursor1.hasNext() )
		{
			avg1 += cursor1.next().getRealFloat();
			avg2 += cursor2.next().getRealFloat();
		}

		avg1 /= (double) numPx;
		avg2 /= (double) numPx;
				
		//
		// compute cross correlation
		//
		cursor1.reset();
		cursor2.reset();
				
		double var1 = 0, var2 = 0;
		double coVar = 0;
		
		while ( cursor1.hasNext() )
		{
			final float pixel1 = cursor1.next().getRealFloat();
			final float pixel2 = cursor2.next().getRealFloat();
			
			final double dist1 = pixel1 - avg1;
			final double dist2 = pixel2 - avg2;

			coVar += dist1 * dist2;
			var1 += dist1 * dist1;
			var2 += dist2 * dist2;
		}		
		
		var1 /= (double) numPx;
		var2 /= (double) numPx;
		coVar /= (double) numPx;

		double stDev1 = Math.sqrt(var1);
		double stDev2 = Math.sqrt(var2);

		// all pixels had the same color....
		if (stDev1 == 0 || stDev2 == 0)
		{
			if ( stDev1 == stDev2 && avg1 == avg2 )
				return 1;
			else
				return 0;
		}

		// compute correlation coeffienct
		return (float)(coVar / (stDev1 * stDev2));
	}
}
