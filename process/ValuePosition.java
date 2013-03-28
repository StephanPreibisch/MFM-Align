package process;

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