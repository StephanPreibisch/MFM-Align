package process;

import ij.ImagePlus;

import java.util.concurrent.atomic.AtomicInteger;

import mpicbg.imglib.algorithm.fft.Bandpass;
import mpicbg.imglib.algorithm.fft.FourierTransform;
import mpicbg.imglib.algorithm.fft.FourierTransform.PreProcessing;
import mpicbg.imglib.cursor.LocalizableByDimCursor;
import mpicbg.imglib.cursor.LocalizableCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.display.Display;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.multithreading.SimpleMultiThreading;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyValueFactory;
import mpicbg.imglib.type.numeric.complex.ComplexFloatType;
import mpicbg.imglib.type.numeric.real.FloatType;

public class AutoFocus 
{
	public static Image< FloatType > focus( final Image< FloatType > image )
	{
		// create a plane to extract
		final int[] size = new int[ 2 ];
		size[ 0 ] = image.getDimension( 0 );
		size[ 1 ] = image.getDimension( 1 );
		
		final Image< FloatType > plane1 = image.createNewImage( size );
		
		final int numSlices = image.getDimension( 2 );
		
		// compute the first fft to get the size
		final FourierTransform< FloatType, ComplexFloatType > fft1 = new FourierTransform<FloatType, ComplexFloatType>( plane1, new ComplexFloatType() );
		fft1.process();

		final int size2[] = image.getDimensions();
		size2[ 0 ] = fft1.getResult().getDimension( 0 );
		size2[ 1 ] = fft1.getResult().getDimension( 1 );
		
		final Image< FloatType > powerSpectrum = image.createNewImage( size2 );
		
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
												
						final FourierTransform< FloatType, ComplexFloatType > fft = new FourierTransform<FloatType, ComplexFloatType>( plane, new ComplexFloatType() );
						//fft.setPreProcessing( PreProcessing.USE_GIVEN_OUTOFBOUNDSSTRATEGY );
						//fft.setCustomOutOfBoundsStrategy( new OutOfBoundsStrategyValueFactory<FloatType>() );
						fft.process();
						Image< ComplexFloatType> fourier = fft.getResult();
									
						final Bandpass< ComplexFloatType > bandpass = new Bandpass<ComplexFloatType>( fourier, 0, 90 );
						bandpass.process();
						fourier = bandpass.getResult();
						
						final Display< ComplexFloatType > display = fourier.getDisplay();
						final LocalizableByDimCursor< FloatType > randomAccess = powerSpectrum.createLocalizableByDimCursor();
						final LocalizableCursor< ComplexFloatType > cursor = fourier.createLocalizableCursor();
						final int[] tmp = new int[ 3 ];
						
						while ( cursor.hasNext() )
						{
							cursor.fwd();
							cursor.getPosition( tmp );
							tmp[ 2 ] = i;
							
							randomAccess.setPosition( tmp );
							randomAccess.getType().set( display.get32Bit( cursor.getType() ) );
						}
					}
					
					plane.close();
                }
		});

        SimpleMultiThreading.startAndJoin( threads );
		
		return powerSpectrum;
	}
	
	
}
