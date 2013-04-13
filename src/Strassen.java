// Copyright Omri Floman. Downloaded from http://www.cs.huji.ac.il/~omrif01/ on April 14, 2013
// by Charles Cooper of Idylwood Technologies.
public class Strassen 
{
	public static double[][] multiply(double[][] A, double[][] B)
	{
		try
		{
			checkInputStrassen(A,B);
		}
		catch (RuntimeException e)
		{
			throw e;
		}
		return strassenRecursive(A,B);
	}

	private static double[][] reconstructAnswer(double[][] r, double[][] s,
			double[][] t, double[][] u) 
	{
		int n = r.length*2;
		double[][] C = new double[n][n];
		copyBack(C,r,0,0);
		copyBack(C,s,0,n/2);
		copyBack(C,t,n/2,0);
		copyBack(C,u,n/2,n/2);
		return C;
	}

	private static void copyBack(double[][] C, double[][] r, int x, int y) 
	{
		for (int i=0; i<r.length; i++)
		{
			for (int j=0; j<r.length; j++)
			{
				C[x+i][y+j] = r[i][j];
			}
		}
	}

	private static void copy(double[][] a, double[][] A, int x, int y) 
	{
		for (int i=0; i<a.length; i++)
		{
			for (int j=0; j<a.length; j++)
			{
				a[i][j] = A[x+i][y+j]; 
			}
		}
	}

	private static double[][] strassenRecursive(double[][] A, double[][] B) 
	{
		int n = A.length;
		if (n==1)
		{
			double[][] C = new double[1][1];
			C[0][0]=A[0][0]*B[0][0];
			return C;
		}

		double[][] r,s,t,u, a,b,c,d,e,f,g,h, P1,P2,P3,P4,P5,P6,P7;
		r = new double[n/2][n/2];		s = new double[n/2][n/2];		t = new double[n/2][n/2];		
		u = new double[n/2][n/2];		a = new double[n/2][n/2];		b = new double[n/2][n/2];		
		c = new double[n/2][n/2];		d = new double[n/2][n/2];		e = new double[n/2][n/2];		
		f = new double[n/2][n/2];		g = new double[n/2][n/2];		h = new double[n/2][n/2];
		P1 = new double[n/2][n/2];		P2 = new double[n/2][n/2];		P3 = new double[n/2][n/2];		
		P4 = new double[n/2][n/2];		P5 = new double[n/2][n/2];		P6 = new double[n/2][n/2];		
		P7 = new double[n/2][n/2];		
		copy(a,A,0,0);
		copy(b,A,0,n/2);
		copy(c,A,n/2,0);
		copy(d,A,n/2,n/2);
		copy(e,B,0,0);
		copy(f,B,0,n/2);
		copy(g,B,n/2,0);
		copy(h,B,n/2,n/2);

		P1= strassenRecursive(a, add(f,h,-1)); // P1 = a(f-h) = af-ah
		P2= strassenRecursive(add(a,b,1), h); // P2 = (a+b)h = ah+bh
		P3= strassenRecursive(add(c,d,1), e); // P3 = (c+d)e = ce+de
		P4= strassenRecursive(d, add(g,e,-1)); // P4 = d(g-e) = dg-de
		P5= strassenRecursive(add(a,d,1), add(e,h,1)); // P5 = (a+d)(e+h)=ae+de+ah+dh
		P6= strassenRecursive(add(b,d,-1), add(g,h,1)); // P6 = (b-d)(g+h)=bg-dg+bh-dh
		P7= strassenRecursive(add(a,c,-1), add(e,f,1)); // P7 = (a-c)(e+f)=ae-ce+af-cf

		r = add(add(P5,P4,1),add(P2,P6,-1),-1); // r = P5+P4-P2+P6 = ae+bg
		s = add(P1,P2,1); // s = P1+P2 = af+bh
		t = add(P3,P4,1); // t = P3+P4 = ce+dg
		u = add(add(P5,P1,1),add(P3,P7,1),-1); //u = P5+P1-P3-P7 = cf+dh
		return reconstructAnswer(r,s,t,u);
	}

	private static double[][] add(double[][] A, double[][] B, int signofB)
	{
		int n = A.length;
		double[][] C = new double[n][n];
		for (int i=0; i<n; i++)
		{
			for (int j=0; j<n; j++)
			{
				C[i][j] = A[i][j] + signofB*B[i][j];
			}
		}
		return C;

	}

	private static void checkInputStrassen(double[][] A, double[][] B) 
	{
		int p = A.length;
		if (p==0)
		{
			throw new IllegalArgumentException("Null matrix");
		}
		int n=p;
		while (n>1)
		{
			if (n%2 != 0)
			{
				throw new IllegalArgumentException("Non power of 2 matrix");
			}
			n/=2;
		}

		int q = A[0].length;
		if (q==0)
		{
			throw new IllegalArgumentException("Null matrix");
		}

		if (q!=p)
		{
			throw new IllegalArgumentException("Nonsquare Matrix");
		}

		for (int i=1; i<p; i++)
		{
			if (A[i].length != q)
			{
				throw new IllegalArgumentException("Inconsistent matrix");
			}
		}
		if (B.length != q)
		{
			throw new IllegalArgumentException("Inconsistent dimensions");
		}

		int r = B[0].length;
		if (r!=p)
		{
			throw new IllegalArgumentException("Nonsquare Matrix");
		}
		for (int i=1; i<q; i++)
		{
			if (B[i].length != r)
			{
				throw new IllegalArgumentException("Inconsistent matrix");
			}
		}
	}
}
