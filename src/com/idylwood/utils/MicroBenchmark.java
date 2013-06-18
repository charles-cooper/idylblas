
public class MicroBenchmark
{
	private long time = System.currentTimeMillis();
	public void logTime(final String msg)
	{
		final long old_time = time;
		time = System.currentTimeMillis();
		System.out.println(msg+" "+(time-old_time));
	}
	public static void main(final String args[])
	{
		final MicroBenchmark mb = new MicroBenchmark();
		mb.logTime("start");
		mb.logTime("end");
	}
}

