package process;

import java.util.concurrent.atomic.AtomicInteger;

import mpicbg.imglib.image.Image;
import mpicbg.imglib.multithreading.SimpleMultiThreading;
import mpicbg.imglib.type.numeric.real.FloatType;

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
