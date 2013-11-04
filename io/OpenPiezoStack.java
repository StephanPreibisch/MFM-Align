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
import ij.ImageJ;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.Arrays;

import loci.formats.FormatException;
import mpicbg.imglib.container.array.ArrayContainerFactory;
import mpicbg.imglib.cursor.LocalizableByDimCursor;
import mpicbg.imglib.cursor.LocalizableCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.ImageFactory;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.io.ImageOpener;
import mpicbg.imglib.type.numeric.real.FloatType;
import mpicbg.imglib.util.Util;


public class OpenPiezoStack 
{
	public static Image< FloatType > openPiezo( final String folder, final String tag ) throws FormatException, IOException
	{
		final File dir = new File( folder );
		
		// read many 2d-images if it is a directory
		if ( dir.isDirectory() )
		{
			final String[] files = dir.list( new FilenameFilter() 
			{	
				@Override
				public boolean accept( final File dir, final String name ) 
				{
					final File newFile = new File( dir, name );
					
					// ignore directories and hidden files
					if ( newFile.isHidden() || newFile.isDirectory() )
						return false;
					else if ( name.contains( tag ) )
						return true;
					else
						return false;
				}
			});
			
			Arrays.sort( files );
			
			final int depth = files.length;
			
			final ImageOpener opener = new ImageOpener();			
			final ImageFactory< FloatType > factory = new ImageFactory< FloatType>( new FloatType(), new ArrayContainerFactory() );
			
			final int[] dimIndvidual = opener.openImage( dir.getAbsolutePath() + File.separator + files[ 0 ], factory ).getDimensions();
			final int[] dim = dimIndvidual.clone();
			dim[ 2 ] *= depth;
			
			IJ.log( depth + " files, should all be '" + files[ 0 ] + "' [" + dimIndvidual[ 0 ] + "x" + dimIndvidual[ 1 ] + "x" + dimIndvidual[ 2 ] + " image=Image<FloatType>]" );
			IJ.log( "trying to load " + dim[2] + " planes." );
			
			final Image< FloatType > piezoStack = factory.createImage( dim );
			final LocalizableByDimCursor< FloatType > piezoCursor = piezoStack.createLocalizableByDimCursor();
			final int[] pos = new int[ 3 ];
			
			for ( int i = 0; i < depth; ++i )
			{
				final Image< FloatType > tmp = opener.openImage( dir.getAbsolutePath() + File.separator + files[ i ], new ImageFactory< FloatType>( new FloatType(), new ArrayContainerFactory()));
				//IJ.log( Util.printCoordinates( tmp.getDimensions() )  + " - " + files[ i ] );
				
				if ( tmp.getDimension( 0 ) != dimIndvidual[ 0 ] || tmp.getDimension( 1 ) != dimIndvidual[ 1 ] || tmp.getDimension( 2 ) != dimIndvidual[ 2 ] )
				{
					IJ.log( "--- these dimensions are incompatible. Please check file.");
					return null;
				}
				
				final LocalizableCursor< FloatType > cursor = tmp.createLocalizableCursor();
				
				while ( cursor.hasNext() )
				{
					cursor.fwd();
					cursor.getPosition( pos );
					pos[ 2 ] += i * dimIndvidual[ 2 ];
					
					piezoCursor.setPosition( pos );
					
					piezoCursor.getType().set( cursor.getType() );
				}
			}
			
			return piezoStack;
		}
		else
		{
			System.out.println( "'" + dir.getAbsolutePath() + "' is no directory." );
			return null;
		}
	}	
}
