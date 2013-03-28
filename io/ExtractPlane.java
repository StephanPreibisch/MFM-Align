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
package io;

import ij.IJ;
import mpicbg.imglib.cursor.LocalizableByDimCursor;
import mpicbg.imglib.cursor.LocalizableCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.type.numeric.real.FloatType;

public class ExtractPlane 
{
	public static Image< FloatType > extract( final Image< FloatType > stack, final int index )
	{
		return extract( stack, index % 3, index / 3 );
	}
	
	public static Image< FloatType > extract( final Image< FloatType > stack, final int x, final int y )
	{
		IJ.log( "Extracting tile (" + x + ", " + y + ")" );
		final int[] dimensions = stack.getDimensions();
		
		dimensions[ 0 ] /= 3;
		dimensions[ 1 ] /= 3;
		
		final int offsetX = x * dimensions[ 0 ];
		final int offsetY = y * dimensions[ 1 ];

		//TODO:REMOVE
		//dimensions[ 0 ] = 190;
		//dimensions[ 1 ] = 190;

		final Image< FloatType > extractedPlane = stack.createNewImage( dimensions );
		
		final LocalizableCursor< FloatType > cursor = extractedPlane.createLocalizableCursor();
		final LocalizableByDimCursor< FloatType > randomAccess = stack.createLocalizableByDimCursor();
		
		final int[] tmp = dimensions.clone();
		
		while ( cursor.hasNext() )
		{
			final FloatType t = cursor.next();
			
			cursor.getPosition( tmp );
			tmp[ 0 ] += offsetX;// + 65;
			tmp[ 1 ] += offsetY;
			
			randomAccess.setPosition( tmp );
			
			t.set( randomAccess.getType() );
		}
		
		return extractedPlane;
	}
}
