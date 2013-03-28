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
