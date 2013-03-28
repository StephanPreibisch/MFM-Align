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
