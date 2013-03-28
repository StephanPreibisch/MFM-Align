package process;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;

import mpicbg.imglib.container.array.Array;
import mpicbg.imglib.container.array.ArrayContainerFactory;
import mpicbg.imglib.container.basictypecontainer.FloatAccess;
import mpicbg.imglib.container.basictypecontainer.array.FloatArray;
import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.cursor.LocalizableCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.interpolation.Interpolator;
import mpicbg.imglib.interpolation.InterpolatorFactory;
import mpicbg.imglib.interpolation.lanczos.LanczosInterpolatorFactory;
import mpicbg.imglib.interpolation.linear.LinearInterpolatorFactory;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyValueFactory;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.imglib.util.DevUtil;

public class Alignment 
{
	public static float align1d( final float[] reference, final float[] image2, final double stepSize, final double minPrecision )
	{
		return align1d( reference, image2, stepSize, minPrecision, null );
	}
	
	public static float align1d( final float[] reference, final float[] image2, final double stepSize, final double minPrecision, PrintWriter out )
	{
		final Image< FloatType > refImg = createImageFromArray( reference, new int[]{ reference.length } );
		final Image< FloatType > img2 = createImageFromArray( image2, new int[]{ reference.length } );
	
		final float[] steps = computeSteps( reference.length, stepSize, minPrecision );

		final float lowestValue = get2ndPercentile( img2 );
		final InterpolatorFactory< FloatType > factory = new LanczosInterpolatorFactory<FloatType>( new OutOfBoundsStrategyValueFactory<FloatType>( new FloatType( lowestValue ) ), 5, false );
		final Interpolator< FloatType > interpolator = factory.createInterpolator( img2 );
		
		float offset = 0;
		float bestDifference = computeDifference( refImg, img2.getDimension( 0 ), interpolator, 0 ); 
		
		for ( final float step : steps )
		{
			boolean foundBetter = false;
			
			do
			{
				foundBetter = false;
				
				// try + and - current step
				final float d1 = computeDifference( refImg, img2.getDimension( 0 ), interpolator, offset - step );
				final float d2 = computeDifference( refImg, img2.getDimension( 0 ), interpolator, offset + step );
	
				if ( out != null )
				{
					out.print( -step + " (" + offset + "): " + bestDifference + " <<< " + d1 );
					
					if ( d1 < bestDifference && d1 < d2 )
						out.println( " <--" );
					else
						out.println( "" );
					
					out.print( "+" + step + " (" + offset + "): " + bestDifference + " >>> " + d2 );
		
					if ( d2 < bestDifference && d1 < d2 )
						out.println( " <--" );
					else
						out.println( "" );
				}
				
				if ( d1 < bestDifference || d2 < bestDifference )
				{
					foundBetter = true;
					
					if ( d1 < d2 )
					{
						offset -= step;
						bestDifference = d1;
					}
					else
					{
						offset += step;
						bestDifference = d2;						
					}
				}			
			}
			while( foundBetter );
		}
		
		if ( out != null )
			out.println( "Final offset: " + offset + " (" + bestDifference + ")" );
		
		return offset;
	}
	
	public static float get2ndPercentile( final Image< FloatType > image )
	{
		final float[] tmp = new float[ image.getDimension( 0 ) ];

		int i = 0;
		
		for ( final FloatType t : image )
			tmp[ i++ ] = t.get();
		
		Arrays.sort( tmp );
		
		return tmp[ Math.round( (tmp.length/100.0f) * 5.0f ) ];
	}
	
	public static Image< FloatType > getAlignedSeries( final Image< FloatType > imgA, final float offset )
	{
		final Image< FloatType > imgOut = imgA.createNewImage();
		final LocalizableCursor< FloatType > refCursor = imgOut.createLocalizableCursor();
		
		final float lowestValue = get2ndPercentile( imgA );
		final InterpolatorFactory< FloatType > factory = new LanczosInterpolatorFactory<FloatType>( new OutOfBoundsStrategyValueFactory<FloatType>( new FloatType( lowestValue ) ), 5, false );
		//final InterpolatorFactory< FloatType > factory = new LinearInterpolatorFactory<FloatType>( new OutOfBoundsStrategyValueFactory<FloatType>( new FloatType( lowestValue ) ) );
		final Interpolator< FloatType > interpolator = factory.createInterpolator( imgA );

		final float[] tmp = new float[ 1 ];

		while ( refCursor.hasNext() )
		{
			refCursor.fwd();
			tmp[ 0 ] = refCursor.getPosition( 0 ) + offset;			
			interpolator.setPosition( tmp );
			
			refCursor.getType().set( interpolator.getType().get() );
		}
		
		return imgOut;
	}
	
	public static float computeDifference( final Image< FloatType > imgA, final int sizeImgB, final Interpolator< FloatType > interpolator, final float offset )
	{
		final LocalizableCursor< FloatType > refCursor = imgA.createLocalizableCursor();		
		
		double difference = 0;
		int countInside = 0;
		
		final float[] tmp = new float[ 1 ];
		
		while ( refCursor.hasNext() )
		{
			final float v1 = refCursor.next().get();
			tmp[ 0 ] = refCursor.getPosition( 0 ) + offset;
			
			if ( tmp[ 0 ] >= 0 && tmp[ 0 ] < sizeImgB )
			{
				++countInside;
				interpolator.setPosition( tmp );
				difference += Math.abs( v1 - interpolator.getType().get() );
			}
		}
		
		return (float)(difference / countInside);
	}
	
	protected static float[] computeSteps( final double size, final double stepSize, final double minPrecision )
	{
		// the step size we will use
		final ArrayList< Float > stepList = new ArrayList<Float>();
		
		float s = (float)size/2;
		
		do
		{
			s /= stepSize;
			stepList.add( s );
		}
		while ( s > minPrecision );
				
		final float[] steps = new float[ stepList.size() ];
		
		for ( int i = 0; i < steps.length; ++i )
			steps[ i ] = stepList.get( i );
		
		return steps;
	}
	
	final public static Image<FloatType> createImageFromArray( final float[] data, final int[] dim )
	{
		final FloatAccess access = new FloatArray( data );
		final Array<FloatType, FloatAccess> array = 
			new Array<FloatType, FloatAccess>(new ArrayContainerFactory(), access, dim, 1 );
			
		// create a Type that is linked to the container
		final FloatType linkedType = new FloatType( array );
		
		// pass it to the DirectAccessContainer
		array.setLinkedType( linkedType );
		
		return new Image<FloatType>(array, new FloatType());
	}

}
