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

import java.util.ArrayList;
import java.util.Collection;

import mpicbg.models.AbstractModel;
import mpicbg.models.IllDefinedDataPointsException;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.Point;
import mpicbg.models.PointMatch;

/**
 * Implementation of the Function interface for function fitting.
 * 
 * @author preibischs
 *
 * @param <M>
 */
public abstract class AbstractFunction< M extends AbstractFunction< M > > extends AbstractModel< M > implements Function< Point >
{
	private static final long serialVersionUID = -5628897768456805162L;

	@Override
	public int getMinNumMatches() { return getMinNumPoints(); }

	@Deprecated
	@Override
	public <P extends PointMatch> void fit( final Collection<P> matches ) throws NotEnoughDataPointsException, IllDefinedDataPointsException
	{
		final ArrayList<Point> list = new ArrayList<Point>();
		
		for ( final P pm : matches )
			list.add( pm.getP1() );
		
		fitFunction( list );
	}

	@Override
	public float[] apply( final float[] location ) { return null; }

	@Override
	public void applyInPlace( final float[] location ) {}	
}
