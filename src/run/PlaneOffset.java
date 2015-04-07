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
package run;

import mpicbg.models.Point;

public class PlaneOffset extends Point
{
	private static final long serialVersionUID = 1L;
	
	final MicroscopyPlane plane;
	final float offset;
	
	public PlaneOffset( final MicroscopyPlane plane, final float offset )
	{
		super( new float[]{ offset } );
		this.plane = plane;
		this.offset = offset;
	}
}
