package run;

import ij.ImageJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.io.Opener;
import ij.process.FloatProcessor;
import io.TextFileAccess;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

import mpicbg.imglib.algorithm.gauss.GaussianConvolutionReal;
import mpicbg.imglib.container.array.ArrayContainerFactory;
import mpicbg.imglib.cursor.Cursor;
import mpicbg.imglib.cursor.LocalizableByDimCursor;
import mpicbg.imglib.image.Image;
import mpicbg.imglib.image.ImageFactory;
import mpicbg.imglib.image.display.imagej.ImageJFunctions;
import mpicbg.imglib.multithreading.SimpleMultiThreading;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyMirrorFactory;
import mpicbg.imglib.outofbounds.OutOfBoundsStrategyValueFactory;
import mpicbg.imglib.type.numeric.real.FloatType;

public class AnalyzeTracks 
{
	public static void findAllDataDirs( final String base, final ArrayList< String > data )
	{
		final File f = new File( base );
	
		String[] entries = f.list();
		
		for ( final String entry : entries )
		{
			if ( entry.contains( "track") )
			{
				data.add( f.getAbsolutePath() );
				return;
			}
		}
		
		for ( final String entry : entries )
		{
			final File newF = new File( base, entry );
			
			if ( newF.isDirectory() )
				findAllDataDirs( newF.getAbsolutePath(), data );
		}
	}
	
	public static void writeTracks( ArrayList< float[][] > tracks, final String file, final Image<FloatType> dna )
	{
		PrintWriter out = TextFileAccess.openFileWrite( file );
		final LocalizableByDimCursor< FloatType > cursor = dna.createLocalizableByDimCursor( new OutOfBoundsStrategyValueFactory<FloatType>() );
		
		out.println( "track" + "\t" + "pos" + "\t" + "x"  + "\t" + "y" + "\t" + "z" + "\t" + "t" + "\t" +"DNA");
		
		int trackId = 0;
		for ( final float[][] track : tracks )
		{
			for ( int i = 0; i < track.length; ++i )
			{
				cursor.setPosition( Math.round( track[ i ][ 0 ] ), 0 );
				cursor.setPosition( Math.round( track[ i ][ 1 ] ), 1 );
				cursor.setPosition( Math.round( track[ i ][ 2 ] ), 2 );
				
				out.println( trackId + "\t" + i + "\t" + track[ i ][ 0 ]  + "\t" + track[ i ][ 1 ] + "\t" + track[ i ][ 2 ] + "\t" + track[ i ][ 3 ] + "\t" + cursor.getType().get() );
			}
			++trackId;
		}
		
		out.close();
	}
	
	public static ArrayList< float[] > readData( String filename ) throws IOException
	{
		final BufferedReader in = TextFileAccess.openFileRead( filename );
		
		final ArrayList< float[] > file = new ArrayList< float[] >();
		
		// discard header
		in.readLine();
		
		while ( in.ready() )
		{
			String line = in.readLine().trim();
			
			line = line.replace( "\t", " " );
			
			line = line.replace( "         ", " " );
			
			while ( line.contains( "  ") )
				line = line.replace( "  ", " " );
			
			final String[] entries = line.split( " " );
			final float[] values = new float[ entries.length ];
			
			for ( int i = 0; i < entries.length; ++i )
				values[ i ] = Float.parseFloat( entries[ i ] );
			
			//System.out.println( values.length );
			
			file.add( values );
		}
		
		in.close();
		
		return file;
	}
	
	public static ArrayList< float[][] > getTracks( final ArrayList< float[] > data, final int minLength )
	{
		final ArrayList< float[][] > tracks = new ArrayList< float[][] >();
		
		final int numCols = data.get( 0 ).length / 3;
		final float[] line0 = data.get( 0 );
		
		int valid = 0;
		
		for ( int i = 0; i < numCols; ++i )
		{
			final int time = Math.round( line0[ i*3 ] );
			
			double x, y, z, t;
			int j = 1;
			
			int consec = 0;
			int maxConsc = 0;
			int conscOffset = 0;
			
			do
			{
				final float[] line = data.get( j );
				y = line[ i*3 ];
				x = line[ i*3 + 1 ];
				z = line[ i*3 + 2 ];
				t = time + j - 1;
				
				if ( Math.abs( y - Math.floor( y ) ) > 0 && Math.abs( x - Math.floor( x ) ) > 0 )
				{
					++consec;
				}
				else
				{
					if ( consec > maxConsc )
					{
						maxConsc = consec;
						conscOffset = j - maxConsc;
					}
					consec = 0;
				}
				
				++j;
 			} 
			while ( x != 0 );
			
			if ( maxConsc >= minLength )
			{
				valid++;
				final float[][] track = new float[ maxConsc ][];			
					
				for ( int k = conscOffset; k < conscOffset + track.length; ++k )
				{
					final float[] line = data.get( k );
					y = line[ i*3 ];
					x = line[ i*3 + 1 ];
					z = line[ i*3 + 2 ];
					t = time + k - 1;
				
					track[ k - conscOffset ] = new float[] { (float)x-1, (float)y-1, (float)z-1, (float)t-1 };
					
					//System.out.println( (k - conscOffset) + "/" + track.length + " "+ x + " " + y + " " + z  + " " + t );
	 			}
				
				tracks.add( track );
			}
		}
		
		System.out.println( "Found " +  valid + " valid tracks out of " + numCols );
		
		return tracks;
	}
	
	public static Image< FloatType > extractDNA( final ImagePlus imp )
	{
		final ImageFactory< FloatType > f = new ImageFactory<FloatType>( new FloatType(), new ArrayContainerFactory() );
		final Image< FloatType > img = f.createImage( new int[] {imp.getWidth(), imp.getHeight(), 9 } );
		
		final ImageStack stack = imp.getImageStack();
		
		final Cursor< FloatType > cursor = img.createCursor();
		
		for ( int i = 0; i < 9; ++i )
		{
			final FloatProcessor fp = (FloatProcessor)stack.getProcessor( i * 3 + 1 + 2 );
		
			for ( final float value : (float[])fp.getPixels() )
				cursor.next().set( value );
		}
		
		return img;
	}

	public static Image< FloatType > copy( final ImagePlus imp )
	{
		final ImageFactory< FloatType > f = new ImageFactory<FloatType>( new FloatType(), new ArrayContainerFactory() );
		final Image< FloatType > img = f.createImage( new int[] {imp.getWidth(), imp.getHeight(), 3, 9, 300 } );
		
		final ImageStack stack = imp.getImageStack();
		
		final Cursor< FloatType > cursor = img.createCursor();
		
		for ( int i = 0; i < stack.getSize(); ++i )
		{
			final FloatProcessor fp = (FloatProcessor)stack.getProcessor( i + 1 );
		
			for ( final float value : (float[])fp.getPixels() )
				cursor.next().set( value );
		}
		
		return img;
	}
	
	public static void writeBack( final Image< FloatType > img, final ImagePlus imp )
	{
		final ImageStack stack = imp.getImageStack();
		
		final Cursor< FloatType > cursor = img.createCursor();
		
		for ( int i = 0; i < stack.getSize(); ++i )
		{
			final FloatProcessor fp = (FloatProcessor)stack.getProcessor( i + 1 );
			final float[] pixels = (float[])fp.getPixels();
			
			for ( int j = 0; j < pixels.length; ++j )
				pixels[ j ] = cursor.next().get();
			
			//for ( final float value : (float[])fp.getPixels() )
		//		cursor.next().set( value );
		}
		
	}

	public static void norm( final Image< FloatType > img )
	{
		img.getDisplay().setMinMax();
		final float min  = (float)img.getDisplay().getMin();
		final float max  = (float)img.getDisplay().getMax();
		
		for ( final FloatType t : img )
			t.set( (t.get() - min)/(max-min) );
	}

	public static void main( String[] args ) throws IOException
	{
		new ImageJ();
		
		final ArrayList< String > allDirs = new ArrayList<String>();
		//findAllDataDirs( "/Volumes/TOSHIBA EXT/3D analysis files", allDirs );
		
		findAllDataDirs( "/Volumes/Time Machine Travel/data/cell22/", allDirs ) ;
		
		//allDirs.clear();
		//allDirs.add( "/Users/preibischs/Documents/Microscopy/david/sample 1/" );
		
		for ( final String baseDir : allDirs )
		{
			System.out.println( baseDir );
			
			final File f = new File( baseDir );
			
			String[] entries = f.list();
			ArrayList< String > textFiles = new ArrayList<String>();
			
			for ( final String entry : entries )
			{
				if ( entry.contains( "track") )
				{
					textFiles.add( new File( baseDir, entry ).getAbsolutePath() );
				}
			}
			
			ImagePlus imp = new Opener().openImage( new File( baseDir, "avgcorrected_quantile_aligned.tif").getAbsolutePath() );
			//imp.show();
			
			//final Image<FloatType> img2 = copy( imp );
			//GaussianConvolutionReal< FloatType > conv = new GaussianConvolutionReal<FloatType>( img2, new OutOfBoundsStrategyMirrorFactory<FloatType>(), new double[]{ 1.7, 1.7, 0, 1.7, 0 } );
			//conv.process();
			//writeBack( conv.getResult(), imp );
			//imp.show();
			
			final Image< FloatType > dna = extractDNA( imp );
			norm( dna );
			//ImageJFunctions.show( dna );

			for ( final String file : textFiles )
			{
				System.out.println( "\t" + file );
				
				//final ArrayList< float[] > data = readData( new File( baseDir, "_tracks.txt").getAbsolutePath() );
				final ArrayList< float[] > data = readData( file );
				final ArrayList< float[][] > tracks3 = getTracks( data, 3 );
				final ArrayList< float[][] > tracks6 = getTracks( data, 6 );
				
				writeTracks( tracks3, file + ".3.txt", dna );
				writeTracks( tracks6, file + ".6.txt", dna );
				//writeTracks( tracks, new File( baseDir, "_tracksextracted.txt").getAbsolutePath(), dna );
			}
			
			imp.close();
			dna.close();

			System.out.println( "done." );
		}
		
	}
}
