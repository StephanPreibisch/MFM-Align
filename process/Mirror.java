package process;

import mpicbg.imglib.algorithm.mirror.MirrorImage;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.type.numeric.real.FloatType;

public class Mirror 
{
	public static void horizontal( final Image< FloatType > image )
	{
		final MirrorImage< FloatType > mirror = new MirrorImage<FloatType>( image, 0 );		
		mirror.process();
	}
}
