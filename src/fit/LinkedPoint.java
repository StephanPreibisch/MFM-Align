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
package fit;

import mpicbg.imglib.util.Util;
import mpicbg.models.Point;

/**
 * Extends the Point class (mpicbg package) by linking of an arbitrary Object.
 * 
 * @author preibischs
 *
 * @param <P>
 */
public class LinkedPoint<P> extends Point
{
	private static final long serialVersionUID = 1L;

	final P link;
	
	public LinkedPoint( final float[] l, final P link )
	{
		super( l.clone() );		
		this.link = link;
	}

	public LinkedPoint( final float[] l, final float[] w, final P link )
	{
		super( l.clone(), w.clone() );
		this.link = link;
	}

	public P getLinkedObject() { return link; }
	
	public String toString() { return "LinkedPoint " + Util.printCoordinates( l ); }
}
