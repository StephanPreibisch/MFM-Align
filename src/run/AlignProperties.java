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

/**
 * Defines parameters for the alignment in XY and Z.
 * 
 * @author preibischs
 *
 */
public class AlignProperties 
{
	public final static int numTiles = 9;
	public static double epsilon = 0.2;
	public static double minInlierRatio = 0.5;
	
	public static String tmpName = "tmp_";
	public static String piezoStack = "_piezo.tif";
	public static String piezoProj = "_piezo_avg.tif";
	public static String entropies = "_entropies.txt";
	
	//public static double[] sigma = new double[]{ 0.75, 0.75, 4 };
	public static double[] sigma = new double[]{ 0, 0, 1 };
}
