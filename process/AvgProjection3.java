package process;

import ij.IJ;
import mpicbg.imglib.cursor.LocalizableByDimCursor;
import mpicbg.imglib.cursor.LocalizableCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.util.RealSum;

public class AvgProjection3 
{
	public static Image< FloatType > project( final Image< FloatType > image )
	{
		if ( image.getNumDimensions() != 3 )
		{
			IJ.error( "number of dimensions is not 3, quitting." );
			return null;
		}
		
		final int[] dim = new int[ 2 ];
		dim[ 0 ] = image.getDimension( 0 );
		dim[ 1 ] = image.getDimension( 1 );
		final int d = image.getDimension( 2 );
		
		final Image< FloatType > projection = image.createNewImage( dim );
		
		final LocalizableCursor< FloatType > cursor = projection.createLocalizableCursor();
		final LocalizableByDimCursor< FloatType > randomAccess = image.createLocalizableByDimCursor();
		
		final int[] tmp = new int[ 3 ];
		
		while ( cursor.hasNext() )
		{
			final FloatType target = cursor.next();
			
			cursor.getPosition( tmp );
			tmp[ 2 ] = 0;
			
			randomAccess.setPosition( tmp );
			
			final RealSum sum = new RealSum( d );
			
			for ( int i = 0; i < d; ++i )
			{
				sum.add( randomAccess.getType().get() );
				
				if ( d != i - 1)
					randomAccess.fwd( 2 );
			}
			
			target.set( (float)(sum.getSum() / (double)d) );
		}
		return projection;
	}
}