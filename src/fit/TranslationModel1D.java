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

import mpicbg.models.AbstractModel;
import mpicbg.models.InvertibleCoordinateTransform;
import mpicbg.models.NotEnoughDataPointsException;
import mpicbg.models.PointMatch;

/**
 * 1d-translation model used to calculate the location of the z-planes.
 * 
 * @version 0.2b
 */
public class TranslationModel1D extends AbstractModel< TranslationModel1D > implements InvertibleCoordinateTransform
{
	private static final long serialVersionUID = 7193913463278025602L;

	static final protected int MIN_NUM_MATCHES = 1;
	
	public float tx = 0;
	
	@Override
	final public int getMinNumMatches(){ return MIN_NUM_MATCHES; }
	
	@Override
	final public float[] apply( final float[] l )
	{
		assert l.length == 1 : "1d translation transformations can be applied to 1d points only.";
		
		return new float[]{ l[ 0 ] + tx };
	}
	
	@Override
	final public void applyInPlace( final float[] l )
	{
		assert l.length == 1 : "1d translation transformations can be applied to 1d points only.";
		
		l[ 0 ] += tx;
	}
	
	@Override
	final public float[] applyInverse( final float[] l )
	{
		assert l.length == 1 : "1d translation transformations can be applied to 1d points only.";
		
		return new float[]{ l[ 0 ] - tx };
	}

	@Override
	final public void applyInverseInPlace( final float[] l )
	{
		assert l.length == 1 : "1d translation transformations can be applied to 1d points only.";
		
		l[ 0 ] -= tx;
	}
		
	@Override
	final public < P extends PointMatch >void fit( final Collection< P > matches ) throws NotEnoughDataPointsException
	{
		if ( matches.size() < MIN_NUM_MATCHES ) throw new NotEnoughDataPointsException( matches.size() + " data points are not enough to estimate a 1d translation model, at least " + MIN_NUM_MATCHES + " data points required." );
		
		// center of mass:
		double pcx = 0;
		double qcx = 0;
		
		double ws = 0.0f;
		
		for ( final P m : matches )
		{
			final float[] p = m.getP1().getL(); 
			final float[] q = m.getP2().getW(); 
			
			final double w = m.getWeight();
			ws += w;
			
			pcx += w * p[ 0 ];
			qcx += w * q[ 0 ];
		}
		pcx /= ws;
		qcx /= ws;

		tx = ( float )( qcx - pcx );
	}

	@Override
	public TranslationModel1D copy()
	{
		final TranslationModel1D m = new TranslationModel1D();
		m.tx = tx;
		m.cost = cost;
		return m;
	}
	
	@Override
	final public void set( final TranslationModel1D m )
	{
		tx = m.tx;
		cost = m.getCost();
	}
	
	/**
	 * Initialize the model such that the respective affine transform is:
	 * 
	 * 1 tx
	 * 0 1
	 * 
	 * @param tx
	 */
	final public void set( final float tx )
	{
		this.tx = tx;
	}
	
	@Override
	public TranslationModel1D createInverse()
	{
		final TranslationModel1D ict = new TranslationModel1D();
		
		ict.tx = -tx;	
		ict.cost = cost;
		
		return ict;
	}
}