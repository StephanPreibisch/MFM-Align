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

import java.util.concurrent.atomic.AtomicInteger;

import mpicbg.imglib.image.Image;
import mpicbg.imglib.multithreading.SimpleMultiThreading;
import mpicbg.imglib.type.numeric.real.FloatType;

/**
 * Computes the entropy from the autofocus function.
 * 
 * @author preibischs
 *
 */
public class ComputeEntropy 
{
	public static float[] computeEntropyForSlices( final Image< FloatType > image, final int bins )
	{
		// get min and max of the stack
		image.getDisplay().setMinMax();
		final float min = (float)image.getDisplay().getMin();
		final float max = (float)image.getDisplay().getMax();
		
		// create a plane to extract
		final int[] size = new int[ 2 ];
		size[ 0 ] = image.getDimension( 0 );
		size[ 1 ] = image.getDimension( 1 );
		
		final int numSlices = image.getDimension( 2 );
		final float[] entropies = new float[ numSlices ];
		
		final AtomicInteger ai = new AtomicInteger(0);
		final int numThreads = Runtime.getRuntime().availableProcessors();
        final Thread[] threads = SimpleMultiThreading.newThreads( numThreads );
        for ( int ithread = 0; ithread < threads.length; ++ithread )
            threads[ithread] = new Thread(new Runnable()
            {
                public void run()
                {
                	// Thread ID
                	final int myNumber = ai.getAndIncrement();
                	
                	final Image< FloatType > plane = image.createNewImage( size );
                	
					for ( int i = myNumber; i < numSlices; i = i + numThreads )
					{
						CrossCorrelation.extractPlane( image, plane, i );
						entropies[ i ] = computeEntropy( plane, bins, min, max );
					}
					
					plane.close();
                }
		});

        SimpleMultiThreading.startAndJoin( threads );
        
		/*
		for ( int i = 0; i < numSlices; ++i )
		{
			CrossCorrelation.extractPlane( image, plane, i );
			entropies[ i ] = computeEntropy( plane, bins, min, max );
		}
		 */
		

		return entropies;
	}
	
	public static float computeEntropy( final Image< FloatType > image, final int bins, final float min, final float max )
	{
		final int[] hist = new int[ bins ];
		
		final float diff = max - min;
		
		for ( final FloatType t : image )
			if ( t.get() != 0 )
				++hist[ Math.round( ((t.get() - min)/diff) * (bins-1) ) ];
		
		float entropy = 0;
		
		long sum = 0;
		for ( final int entry : hist )
			sum += entry;
		
		for ( final float entry : hist )
		{
			if ( entry > 0 )
			{
				final float prob = entry / sum;				
				entropy -= prob * (Math.log(prob) / Math.log(2)); 
			}
		}
		
		return entropy;
	}
}
