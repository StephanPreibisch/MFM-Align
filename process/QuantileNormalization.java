package process;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import mpicbg.imglib.cursor.LocalizableByDimCursor;
import mpicbg.imglib.cursor.LocalizableCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.type.numeric.real.FloatType;

public class QuantileNormalization 
{
	public static void normalizeTo( final Image< FloatType > reference, final Image< FloatType > template )
	{
		final float[] values = new float[ reference.getNumPixels() ];
		
		int i = 0;
		for ( final FloatType t : reference )
			values[ i++ ] = t.get();
		
		Arrays.sort( values );
		
		final ArrayList< ValuePosition > templateList = new ArrayList< ValuePosition >();	
		final LocalizableCursor< FloatType > cursor = template.createLocalizableCursor();
		
		while ( cursor.hasNext() )
		{
			cursor.fwd();
			templateList.add( new ValuePosition( cursor.getType().get(), cursor.getPosition() ) );
		}
		
		Collections.sort( templateList );
		
		final LocalizableByDimCursor< FloatType > randomAccess = template.createLocalizableByDimCursor();
		
		for ( i = 0; i < values.length; ++i )
		{
			randomAccess.setPosition( templateList.get( i ).position );
			randomAccess.getType().set( values[ i ] );
		}
	}
}
