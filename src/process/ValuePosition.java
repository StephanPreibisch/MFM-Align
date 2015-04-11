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

/**
 * Utility class for quantile normalization.
 * 
 * @author preibischs
 *
 */
public class ValuePosition implements Comparable< ValuePosition >
{
	public float value;
	public int[] position;
	
	public ValuePosition( final float value, final int[] position )
	{
		this.value = value;
		this.position = position;
	}

	@Override
	public int compareTo( final ValuePosition o )
	{
		if ( o.value < value )
			return 1;
		else if ( o.value == value )
			return 0;
		else
			return -1;
	}
}