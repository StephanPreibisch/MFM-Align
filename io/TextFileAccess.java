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
package io;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class TextFileAccess 
{
	public static BufferedReader openFileRead(final File file)
	{
		BufferedReader inputFile;
		try
		{
			inputFile = new BufferedReader(new FileReader(file));
		}
		catch (IOException e)
		{
			System.out.println("TextFileAccess.openFileRead(): " + e);
			inputFile = null;
		}
		return (inputFile);
	}

	public static BufferedReader openFileRead(final String fileName)
	{
		BufferedReader inputFile;
		try
		{
			inputFile = new BufferedReader(new FileReader(fileName));
		}
		catch (IOException e)
		{
			System.out.println("TextFileAccess.openFileRead(): " + e);
			inputFile = null;
		}
		return (inputFile);
	}

	public static PrintWriter openFileWrite(final File file)
	{
		PrintWriter outputFile;
		try
		{
			outputFile = new PrintWriter(new FileWriter(file));
		}
		catch (IOException e)
		{
			System.out.println("TextFileAccess.openFileWrite(): " + e);
			outputFile = null;
		}
		return (outputFile);
	}
	
	public static PrintWriter appendFileWrite( final File file )
	{
		PrintWriter outputFile;
		try
		{
			outputFile = new PrintWriter(new FileWriter(file, true));
		}
		catch (IOException e)
		{
			System.out.println("TextFileAccess.appendFileWrite(): " + e);
			outputFile = null;
		}
		return (outputFile);		
	}
	
	public static PrintWriter openFileWrite(final String fileName)
	{
		PrintWriter outputFile;
		try
		{
			outputFile = new PrintWriter(new FileWriter(fileName));
		}
		catch (IOException e)
		{
			System.out.println("TextFileAccess.openFileWrite(): " + e);
			outputFile = null;
		}
		return (outputFile);
	}
	
	public static PrintWriter openFileWriteEx(final File file) throws IOException
	{
		return new PrintWriter(new FileWriter(file));
	}

	public static BufferedReader openFileReadEx(final File file) throws IOException
	{
		return new BufferedReader(new FileReader(file));
	}

	public static PrintWriter openFileWriteEx(final String fileName) throws IOException
	{
		return new PrintWriter(new FileWriter(fileName));
	}

	public static BufferedReader openFileReadEx(final String fileName) throws IOException
	{
		return new BufferedReader(new FileReader(fileName));
	}
}
