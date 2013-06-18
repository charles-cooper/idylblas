package com.idylwood;

import com.idylwood.matrix.DoubleMatrix1D;
import com.idylwood.matrix.DoubleMatrix2D;
import com.idylwood.matrix.DoubleMatrix2D.Pointer;
import com.idylwood.utils.MathUtils;
import com.idylwood.Blas;

public class FastBlas implements Blas
{

	/**
	  Returns the sum of absolute values; <tt>|x[0]| + |x[1]| + ... </tt>.
	  In fact equivalent to <tt>x.aggregate(cern.jet.math.Functions.plus, cern.jet.math.Functions.abs)</tt>.
	  @param x the first vector.
	  */
	@Override public double dasum(DoubleMatrix1D x){



		double sum = 0;
		double err = 0;
		final int unroll = 6; // empirically it doesn't get much better than this
		final int len = x.length - x.length%unroll;

		// unroll the loop. due to IEEE 754 restrictions
		// the JIT shouldn't be allowed to unroll it dynamically, so it's
		// up to us to do it by hand ;)

		int i = 0;
		double temp[] = new double[x.length];

		for (; i < x.length; i++)
		{
			temp[i] = MathUtils.abs(x.get(i));
		}

		i = 0;
		for (; i < len; i+=unroll)
		{
			final double val = temp[i] + temp[i+1]
				+ temp[i+2] + temp[i+3]
				+ temp[i+4] + temp[i+5];
			final double partial = val - err;
			final double hi = sum + val;
			err = (hi - sum) - partial;
			sum = hi;
		}
		for (; i < x.length; i++)
		{
			final double val = temp[i];
			final double partial = val - err;
			final double hi = sum + val;
			err = (hi - sum) - partial;
			sum = hi;
		}
		return sum;
	}
	/**
	  Combined vector scaling; <tt>y = y + alpha*x</tt>.
	  In fact equivalent to <tt>y.assign(x,cern.jet.math.Functions.plusMult(alpha))</tt>.

	  @param alpha a scale factor.
	  @param x the first source vector.
	  @param y the second source vector, this is also the vector where results are stored.

	  @throws IllegalArgumentException <tt>x.size() != y.size()</tt>..
	  */
	@Override public void daxpy(double alpha, DoubleMatrix1D x, DoubleMatrix1D y){
		if ( x.length != y.length){
			throw new ArrayIndexOutOfBoundsException("Vectors are of different dimension.");
		}
		else{
			int i = 0;
			for(; i < x.length; i++)
			{
				y.set(i, (y.get(i)+ (alpha*(x.get(i)))));
			}
		}
	}
	/**
	  Combined matrix scaling; <tt>B = B + alpha*A</tt>.
	  In fact equivalent to <tt>B.assign(A,cern.jet.math.Functions.plusMult(alpha))</tt>.

	  @param alpha a scale factor.
	  @param A the first source matrix.
	  @param B the second source matrix, this is also the matrix where results are stored.

	  @throws IllegalArgumentException if <tt>A.columns() != B.columns() || A.rows() != B.rows()</tt>.
	  */
	@Override public void daxpy(double alpha, DoubleMatrix2D A, DoubleMatrix2D B){
		if(A.rows() != B.rows() || A.cols() != B.cols())
			throw new IllegalArgumentException ("Matricies are of different dimension.");
		Pointer Aptr = A.ptr();
		Pointer Bptr = B.ptr();
		for(; Aptr.hasNext(); Aptr.increment(), Bptr.increment())
			Bptr.setValue(Bptr.getValue() + (Aptr.getValue()*alpha));
	}
	/**
	  Vector assignment (copying); <tt>y = x</tt>.
	  In fact equivalent to <tt>y.assign(x)</tt>.

	  @param x the source vector.
	  @param y the destination vector.

	  @throws IllegalArgumentException <tt>x.size() != y.size()</tt>.
	  */
	@Override public void dcopy(DoubleMatrix1D x, DoubleMatrix1D y){
		if ( x.length != y.length){
			throw new ArrayIndexOutOfBoundsException("Vectors are of different dimension.");
		}
		else{
			int i = 0;
			for(; i < x.length; i++)
			{
				y.set(i, (x.get(i)));
			}
		}
	}
	/**
	  Matrix assignment (copying); <tt>B = A</tt>.
	  In fact equivalent to <tt>B.assign(A)</tt>.

	  @param A the source matrix.
	  @param B the destination matrix.

	  @throws IllegalArgumentException if <tt>A.columns() != B.columns() || A.rows() != B.rows()</tt>.
	  */
	@Override public void dcopy(DoubleMatrix2D A, DoubleMatrix2D B){
		if(A.rows() != B.rows() || A.cols() != B.cols())
			throw new IllegalArgumentException ("Matricies are of different dimension.");
		DoubleMatrix2D.copy(A, B);
	}
	/**
	  Returns the dot product of two vectors x and y, which is <tt>Sum(x[i]*y[i])</tt>.
	  In fact equivalent to <tt>x.zDotProduct(y)</tt>.
	  @param x the first vector.
	  @param y the second vector.
	  @return the sum of products.

	  @throws IllegalArgumentException if <tt>x.size() != y.size()</tt>.
	  */
	@Override public double ddot(DoubleMatrix1D x, DoubleMatrix1D y){
		if ( x.length != y.length){
			throw new ArrayIndexOutOfBoundsException("Vectors are of different dimension.");
		}
		else{
			int i = 0;
			double product = 0;
			for(; i < x.length; i++)
			{
				product += (y.get(i)*x.get(i));
			}
			return product;
		}
	}
	/**
	  Generalized linear algebraic matrix-matrix multiply; <tt>C = alpha*A*B + beta*C</tt>.
	  In fact equivalent to <tt>A.zMult(B,C,alpha,beta,transposeA,transposeB)</tt>.
Note: Matrix shape conformance is checked <i>after</i> potential transpositions.

@param transposeA set this flag to indicate that the multiplication shall be performed on A'.
@param transposeB set this flag to indicate that the multiplication shall be performed on B'.
@param alpha a scale factor.
@param A the first source matrix.
@param B the second source matrix.
@param beta a scale factor.
@param C the third source matrix, this is also the matrix where results are stored.

@throws IllegalArgumentException if <tt>B.rows() != A.columns()</tt>.
@throws IllegalArgumentException if <tt>C.rows() != A.rows() || C.columns() != B.columns()</tt>.
@throws IllegalArgumentException if <tt>A == C || B == C</tt>.
*/
	@Override public void dgemm(boolean transposeA, boolean transposeB, double alpha, DoubleMatrix2D A, DoubleMatrix2D B, double beta, DoubleMatrix2D C){
		if (B.rows() != A.cols())
			throw new IllegalArgumentException ("Multipication undefined.");
		if (C.rows() != A.rows() || C.cols() != B.cols())
			throw new IllegalArgumentException ("C has wrong dimensions .");
		if (A == C || B == C)
			throw new IllegalArgumentException ("A or B is bad.");

		if (transposeA)
			A = MathUtils.transpose(A);
		if (transposeB)
			B = MathUtils.transpose(B);
		for ( final Pointer ptr = C.ptr(); ptr.hasNext(); ptr.increment())
			ptr.setValue(ptr.getValue() * beta);
		daxpy(alpha, MathUtils.matrixMultiply(A, B), C);
	}

	/**
	 * Generalized linear algebraic matrix-vector multiply; <tt>y = alpha*A*x + beta*y</tt>.
	 * In fact equivalent to <tt>A.zMult(x,y,alpha,beta,transposeA)</tt>.
	 * Note: Matrix shape conformance is checked <i>after</i> potential transpositions.
	 *
	 * @param transposeA set this flag to indicate that the multiplication shall be performed on A'.
	 * @param alpha a scale factor.
	 * @param A the source matrix.
	 * @param x the first source vector.
	 * @param beta a scale factor.
	 * @param y the second source vector, this is also the vector where results are stored.
	 * 
	 * @throws IllegalArgumentException <tt>A.columns() != x.size() || A.rows() != y.size())</tt>..
	 */
	@Override public void dgemv(boolean transposeA, double alpha, DoubleMatrix2D A, DoubleMatrix1D x, double beta, DoubleMatrix1D y){
		if (A.cols() != x.length || A.rows() != y.length)
			throw new IllegalArgumentException ("Bad dimensions.");
		if(transposeA)
			A = MathUtils.transpose(A);
		A.multiply(x).scale(alpha) // alpha * A * x
			.plus(y.scale(beta)) // + y*beta
			.copy(y); // copy in
	}
	/**
	 * Performs a rank 1 update; <tt>A = A + alpha*x*y'</tt>.
	 * Example:
	 * <pre>
	 * A = { {6,5}, {7,6} }, x = {1,2}, y = {3,4}, alpha = 1 -->
	 * A = { {9,9}, {13,14} }
	 * </pre>

	 * @param alpha a scalar.
	 * @param x an m element vector.
	 * @param y an n element vector.
	 * @param A an m by n matrix.
	 */
	@Override public void dger(double alpha, DoubleMatrix1D x, DoubleMatrix1D y, DoubleMatrix2D A){
		if (x.length != A.rows())
			throw new IllegalArgumentException("Bad rows");
		if (y.length != A.cols())
			throw new IllegalArgumentException("Bad columns");
		for (final Pointer ptr = A.ptr(); ptr.hasNext(); ptr.increment())
			ptr.setValue(ptr.getValue() + alpha*x.get(ptr.currentColumn())*y.get(ptr.currentRow()));
	}
	/**
	  Return the 2-norm; <tt>sqrt(x[0]^2 + x[1]^2 + ...)</tt>.
	  In fact equivalent to <tt>Math.sqrt(Algebra.DEFAULT.norm2(x))</tt>.

	  @param x the vector.
	  */
	@Override public double dnrm2(DoubleMatrix1D x){

		double num = 0; 
		for(int i = 0; i < x.length ; i++){
			num += x.get(i);
		}
		return num;
	}
	/**
	  Applies a givens plane rotation to (x,y); <tt>x = c*x + s*y; y = c*y - s*x</tt>.
	  @param x the first vector.
	  @param y the second vector.
	  @param c the cosine of the angle of rotation.
	  @param s the sine of the angle of rotation.
	  */
	@Override public void drot(DoubleMatrix1D x, DoubleMatrix1D y, double c, double s){

	}
	/**
	  Constructs a Givens plane rotation for <tt>(a,b)</tt>.
	  Taken from the LINPACK translation from FORTRAN to Java, interface slightly modified.
	  In the LINPACK listing DROTG is attributed to Jack Dongarra

	  @param  a  rotational elimination parameter a.
	  @param  b  rotational elimination parameter b.
	  @param  rotvec[]  Must be at least of length 4. On output contains the values <tt>{a,b,c,s}</tt>.
	  */
	@Override public void drotg(double a, double b, double rotvec[]){

	}
	/**
	  Vector scaling; <tt>x = alpha*x</tt>.
	  In fact equivalent to <tt>x.assign(cern.jet.math.Functions.mult(alpha))</tt>.

	  @param alpha a scale factor.
	  @param x the first vector.
	  */
	@Override public void dscal(double alpha, DoubleMatrix1D x){
		for(int i = 0; i < x.length; i++){
			x.set(i, x.get(i)*alpha);
		}
	}
	@Override public void dscal(double alpha, DoubleMatrix2D A)
	{
		final DoubleMatrix2D.Pointer ptr = A.new Pointer();
		for (; ptr.hasNext(); ptr.increment())
			ptr.setValue(ptr.getValue() * alpha);
	}
	/**
	  Matrix scaling; <tt>A = alpha*A</tt>.
	  In fact equivalent to <tt>A.assign(cern.jet.math.Functions.mult(alpha))</tt>.

	  @param alpha a scale factor.
	  @param A the matrix.
	  */
	@Override public void dswap(DoubleMatrix1D x, DoubleMatrix1D y){
		double holder = 0;
		for (int i = 0; i < x.length; i++){
			holder = y.get(i);
			y.set(i, x.get(i));
			x.set(i, holder);
		}
	}
	/**
	  Swaps the elements of two matrices; <tt>B <==> A</tt>.

	  @param A the first matrix.
	  @param B the second matrix.

	  @throws IllegalArgumentException if <tt>A.columns() != B.columns() || A.rows() != B.rows()</tt>.
	  */
	@Override public void dswap(DoubleMatrix2D x, DoubleMatrix2D y){
		if (x.cols()!=y.cols() || x.rows()!=y.rows())
			throw new IllegalArgumentException("Bad dimensions!");
		double tmp;
		Pointer x_ptr = x.ptr();
		Pointer y_ptr = y.ptr();
		for (; x_ptr.hasNext(); x_ptr.increment(), y_ptr.increment())
		{
			tmp = y_ptr.getValue();
			y_ptr.setValue(x_ptr.getValue()); 
			x_ptr.setValue(tmp);
		}
	}
	/**
	  Symmetric matrix-vector multiplication; <tt>y = alpha*A*x + beta*y</tt>.
	  Where alpha and beta are scalars, x and y are n element vectors and

	  A is an n by n symmetric matrix.
	  A can be in upper or lower triangular format.
	  @param isUpperTriangular is A upper triangular or lower triangular part to be used?
	  @param alpha scaling factor.
	  @param A the source matrix.
	  @param x the first source vector.
	  @param beta scaling factor.
	  @param y the second vector holding source and destination.
	  */
	@Override public void dsymv(boolean isUpperTriangular, double alpha, DoubleMatrix2D A, DoubleMatrix1D x, double beta, DoubleMatrix1D y)
	{

	}
	/**
	  Triangular matrix-vector multiplication; <tt>x = A*x</tt> or <tt>x = A'*x</tt>.
	  Where x is an n element vector and A is an n by n unit, or non-unit,
	  upper or lower triangular matrix.
	  @param isUpperTriangular is A upper triangular or lower triangular?
	  @param transposeA set this flag to indicate that the multiplication shall be performed on A'.
	  @param isUnitTriangular true --> A is assumed to be unit triangular; false --> A is not assumed to be unit triangular
	  @param A the source matrix.
	  @param x the vector holding source and destination.
	  */
	@Override public void dtrmv(boolean isUpperTriangular, boolean transposeA, boolean isUnitTriangular, DoubleMatrix2D A, DoubleMatrix1D x){

	}
	/**
	  Returns the index of largest absolute value; <tt>i such that |x[i]| == max(|x[0]|,|x[1]|,...).</tt>.

	  @param x the vector to search through.
	  @return the index of largest absolute value (-1 if x is empty).
	  */
	@Override public int idamax(DoubleMatrix1D x)
	{
		double check = MathUtils.abs(x.get(0));
		int idx = 0; 
		for (int i = 0; i < x.length; i++){
			if(MathUtils.abs(x.get(i)) > check){
				check = MathUtils.abs(x.get(i));
				idx = i;
			}
		}
		return idx; 
	}

	public static void main(String[] args)
	{
		System.out.println("hello!");
		double testdata[] = {1,1,0}; 
		double testdataz[][] = {{1,1,1},{1,1,1},{1,1,1}};
		double testdatax[] = {0,0,0};



		DoubleMatrix1D test = new DoubleMatrix1D(testdata); 
		DoubleMatrix2D test2 = new DoubleMatrix2D(testdataz);
		DoubleMatrix1D test3 = new DoubleMatrix1D(testdatax );


		new FastBlas().dgemv(false,1, test2, test, 1, test3);



		for(int i = 0; i < test3.length; i++){
			System.out.println(test3.get(i));
		}
	}
}


