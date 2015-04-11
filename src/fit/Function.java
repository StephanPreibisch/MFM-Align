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

import java.util.Collection;

import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.Point;

/**
 * Interface for a function that can be fitted to a list of Points (mpicbg package).
 * 
 * @author Stephan Preibisch
 *
 * @param <P> - if a special extension of {@link Point} is necessary, otherwise just implement Function<Point>
 */
public interface Function< P extends Point >
{
	/**
	 * @return - how many points are at least necessary to fit the function
	 */
	public int getMinNumPoints();

	/**
	 * Fits this Function to the set of {@link Point}s.

	 * @param points - {@link Collection} of {@link Point}s
	 * @throws NotEnoughDataPointsException - thrown if not enough {@link Point}s are in the {@link Collection}
	 */
	public void fitFunction( final Collection<P> points ) throws NotEnoughDataPointsException;
	
	/**
	 * Computes the minimal distance of a {@link Point} to this function
	 *  
	 * @param point - the {@link Point}
	 * @return - distance to the {@link Function}
	 */
	public double distanceTo( final P point );

}
