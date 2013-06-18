/*
 * ====================================================
 * Copyright (C) 2013 by Idylwood Technologies, LLC. All rights reserved.
 *
 * Developed at Idylwood Technologies, LLC.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * The License should have been distributed to you with the source tree.
 * If not, it can be found at
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Author: Charles Cooper
 * Date: 2013
 * ====================================================
 */
package com.idylwood.utils;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;
import java.util.Random;

import com.idylwood.matrix.DoubleMatrix2D;

// This class contains a bunch of numerical methods which are intended to be fast and precise.
// Several of the methods come in 'fast', 'normal' and 'slow' versions to give the end user
// more control over the trade-off between accuracy and speed.
// Most of the methods are marked final as a matter of style
// (and also intended to avoid vtable lookups and let the JIT inline better).
// TODO write unit tests!

public final class MathUtils {
	// remove public instantiation
	private MathUtils() {}

	// TODO put this thing somewhere else
	public static class LinearRegression {
		public final double intercept, slope;
		public LinearRegression(final double intercept, final double slope) {
			this.intercept = intercept; this.slope = slope;
		}
	}

	// TODO get rid of garbage
	public static final LinearRegression regress(final double [] x, final double [] y)
	{
		final double intercept, slope;
		final double xMean = mean(x);
		final double yMean = mean(y);
		final double [] xCentered = shift(x,-xMean);
		final double [] yCentered = shift(y,-yMean);
		slope = sum( multiply(xCentered,yCentered) ) / sum (pow(xCentered,2));
		intercept = yMean - slope * xMean;
		return new LinearRegression(intercept,slope);
	}

	// Faster and more correct (round away from zero instead of round towards +infinity) than java.lang.Math.round.
	private static final long ONE_HALF = Double.doubleToRawLongBits(0.5); // bits of 0.5
	/**
	 * Faster and more correct version of round than {@link java.lang.Math#round(double)}.
	 * Rounds away from zero instead of round towards infinity.
	 * In most cases it will have the same return as Math.round,
	 * but for numbers like -0.5, -1.5, this will round down to -1, -2,
	 * whereas JDK specification is to round up to 0, -1.
	 * It also has no branching and so is faster.
	 * This agrees with glibc.
	 * @param d
	 * @return
	 */
	public final static long round(final double d)
	{
		final long l = ONE_HALF | sign(d); // equivalent to d < 0 ? -0.5 : 0.5;
		return (long)(d + Double.longBitsToDouble(l));
	}
	/**
	 * Returns a long containing the sign of d.
	 * Returns 1L<<63 if d < 0 and 0 otherwise.
	 * @param d
	 * @return
	 */
	public static final long sign(final double d)
	{
		return Double.doubleToRawLongBits(d) & Long.MIN_VALUE;
	}

	public final static long[] round(final double[] d)
	{
		final long[] ret = new long[d.length];
		for (int i = 0; i < d.length; i++)
			ret[i] = round(d[i]);
		return ret;
	}

	// TODO do we need more InPlace methods? they are significantly faster because of no allocation/copying overhead
	public final static void roundInPlace(final double[] d)
	{
		for (int i = 0; i < d.length; ++i)
		{
			d[i] = (double)round(d[i]);
		}
	}

	public final static double roundToCent(final double d)
	{
		return round(d*100.0) / 100.0;
	}

	// numerically stable calculation of mean
	public final static double mean(final double [] values)
	{
		return sum(values) / values.length;
	}

	public final static double meanFast(final double[] values)
	{
		return sumFast(values) / values.length;
	}

	public final static double meanSlow(final double[] values)
	{
		return sumSlow(values) / values.length;
	}

	public final static double max(double... values)
	{
		if (values.length==0) return Double.NaN;
		double ret = values[0];
		for (int i = 1; i < values.length; i++)
			if (values[i] > ret)
				ret = values[i];
		return ret;
	}

	public final static int max(int... values)
	{
		if (values.length==0) return Integer.MIN_VALUE;
		int ret = values[0];
		for (int i = 1; i < values.length; i++)
			if (values[i] > ret)
				ret = values[i];
		return ret;
	}

	public final static double min(double... values)
	{
		if (values.length==0) return Double.NaN;
		double ret = values[0];
		for (int i = 1; i < values.length; i++)
			if (values[i] < ret)
				ret = values[i];
		return ret;
	}

	// not for production code
	static final double multiplyAndSumSlow(final double[]values, final double scale)
	{
		BigDecimal ret = new BigDecimal(0);
		for (double x : values)
			ret = ret.add(new BigDecimal(x),MathContext.UNLIMITED);
		return ret.multiply(new BigDecimal(scale),MathContext.UNLIMITED).doubleValue();
	}
	// No side effects.
	// TODO implement without garbage
	static public final double multiplyAndSum(final double[]d, final double scale)
	{
		return sum(scale(d,scale));
	}

	// Note: this seems to be less accurate than multiplyAndSum
	static public final double multiplyAndSumFast(final double[]d,final double scale)
	{
		return scale*sum(d);
	}

	// for comparison purposes
	/*
	static public final double varianceApache(final double[] values)
	{
		return new org.apache.commons.math3.stat.descriptive.moment.Variance().evaluate(values);
	}
	*/

	// Alternative implementation of variance. Not sure which one is more precise.
	static public final double varianceTwo(final double [] values)
	{
		final long n = values.length;
		final long n1 = values.length - 1;
		final double sumOfSquares = sum(pow(values,2)) / n1;
		final double squareOfSum = Math.pow(sum(values),2) / (n*n1);
		return sumOfSquares - squareOfSum;
	}

	// if the mean is precalculated, no point in calculating it again!
	static public final double variance(final double [] values, final double mean)
	{
		final double [] squares = pow(shift(values,-mean),2);
		return sum(squares) / (squares.length - 1);
	}

	static public final double stdev(final double [] values, final double mean)
	{
		return Math.sqrt(variance(values,mean));
	}
	static public final double stdev(final double[] values)
	{
		final double mean = mean(values);
		return stdev(values,mean);
	}

	// Numerically stable calculation of variance
	// whose accuracy doesn't suffer for large n
	// IIRC this matches with varianceApache in tests
	// but it is much faster.
	static public final double variance(final double [] values)
	{
		return variance(values,mean(values));
	}

	// calculation of variance which is somewhat numerically stable
	public static final double variancePopulation(final double [] values)
	{
		final double mean = mean(values);
		final double [] centered = shift(values,-mean);
		final double [] squares = pow(centered,2);
		return sum(squares) / squares.length;
	}

	// theoretically same as R's 'diff' function
	final public static double [] diff(final double[] data)
	{
		double [] ret = new double[data.length - 1];
		for (int i = data.length - 1; i--!=0; )
			ret[i] = data[i+1] - data[i];
		return ret;
	}

	// This will not throw an exception if the arrays are not of equal length,
	// it will return an array whose length is the smaller of the two array lengths
	final public static double [] subtract(final double[] first, final double[] second)
	{
		final int len = Math.min(first.length,second.length);
		final double ret[] = new double[len];
		for (int i = 0; i < len; i++)
			ret[i] = first[i] - second[i];
		return ret;
	}

	/**
	 * Returns newly allocated array which is the result
	 * of termwise addition of the input arrays
	 * @param first
	 * @param second
	 * @throws ArrayIndexOutOfBoundsException if the arrays are of unequal length
	 * @return
	 * Side Effects: Allocation of the returned array
	 */
	final public static double [] add(final double[] first, final double[] second)
	{
		final int len = Math.min(first.length,second.length);
		final double ret[] = new double[len];
		for (int i = 0; i < len; i++)
			ret[i] = first[i]+second[i];
		return ret;
	}

	/**
	 * Returns newly allocated array which is the result
	 * of termwise multiplication of the input arrays
	 * @param first
	 * @param second
	 * @throws ArrayIndexOutOfBoundsException if the arrays are of unequal length
	 * @return
	 * Side Effects: Allocation of the returned array
	 */
	final public static double [] multiply(final double [] x, final double[] y)
	{
		final int len = x.length;
		if (len!=y.length)
			throw new ArrayIndexOutOfBoundsException("Tried to multiply two vectors of unequal length!");
		final double ret[] = new double[len];
		for (int i = 0; i < len; i++)
			ret[i] = x[i]*y[i];
		return ret;
	}

	/**
	 * Returns a newly allocated array whose elements
	 * are the termwise log of the input
	 * @param data
	 * @return
	 * Side Effects: Allocation of a new array
	 */
	final public static double [] log(final double[] data)
	{
		final double [] ret = new double[data.length];
		for (int i = data.length; i-- != 0; )
			ret[i] = Math.log(data[i]);
		return ret;
	}

	/**
	 * Calculates the mean to arbitrary precision. Internally uses
	 * BigDecimals so it suffers from both speed and memory.
	 * @param data
	 * @return
	 */
	final public static double meanArbitraryPrecision(final double data[])
	{ BigDecimal mean = new BigDecimal(0);
		for (double x : data)
			mean = mean.add(new BigDecimal(x),MathContext.UNLIMITED);
		mean = mean.divide(new BigDecimal(data.length),MathContext.UNLIMITED);
		return mean.doubleValue();
	}

	// funnily enough, even though this uses BigDecimals,
	// it has a bug in it and spits out wrong answers sometimes.
	// (I think the bug is in how it handles division and
	// repeating rational numbers
	private final static double varianceSlow(final double[] data)
	{
		BigDecimal mean = new BigDecimal(0);
		for (double x : data)
			mean = mean.add(new BigDecimal(x),MathContext.UNLIMITED);
		mean = mean.divide(new BigDecimal(data.length),MathContext.UNLIMITED);
		//mean = new BigDecimal(mean(data));
		BigDecimal ret = new BigDecimal(0);
		for (double x : data)
		{
			//BigDecimal summand = ret.add(new BigDecimal(x),MathContext.UNLIMITED);
			BigDecimal summand = new BigDecimal(x).subtract(mean,MathContext.UNLIMITED);
			ret = ret.add(summand.pow(2));
		}
		ret = ret.divide(new BigDecimal(data.length - 1),MathContext.DECIMAL128);
		return ret.doubleValue();
	}

	// numerically stable calculation of standard deviation
	public static final double stdevPopulation(final double [] values)
	{
		return Math.sqrt(variance(values));
	}

	/**
	 * Performs termwise Math.pow(double,double) on the elements.
	 * @param values
	 * @param exp
	 * @return Newly allocated array whose elements are
	 * termwise exponentiations of the input.
	 * Side Effects: Allocation of new array
	 */
	public static final double [] pow(final double[] values, final double exp)
	{
		final double[] ret = new double[values.length];
		for (int i = values.length; i--!=0; )
			ret[i] = Math.pow(values[i],exp);
		return ret;
	}

	/**
	 * Shifts all the elements of <code>values</code> by <code>constant</code>.
	 * @param values
	 * @param constant
	 * @return Newly allocated array whose values are values[i]+constant
	 * Side Effects: Allocation of new array
	 */
	public static final double [] shift(final double[] values, final double constant)
	{
		final double ret[] = new double[values.length];
		for (int i = values.length; i--!=0; )
			ret[i] = values[i] + constant;
		return ret;
	}

	// Returns a newly allocated array with all the values
	// of the original array multiplied by the scale.
	// TODO maybe rename this 'multiply'?
	public static final double [] scale(double[] values, double scale)
	{
		double [] ret = new double[values.length];
		for (int i = values.length; i--!=0;)
			ret[i] = values[i]*scale;
		return ret;
	}

	/** Numerically precise implementation of sum
	 * Optimized version of an implementation of Schewchuk's algorithm
	 * which keeps full precision by keeping O(n) space
	 * for the error, unlike Kahan's algorithm which keeps O(1) space.
	 * The tradeoff is that this is fully precise, but Kahan's algorithm
	 * is almost always precise anyways. It is about 12x slower
	 * than the naive implementation, but in turn about 10x faster than calculating
	 * the sum to full precision and then truncating.
	 * @param values
	 * @return Sum of the values
	 */
	public final static double sumSlow(double... values)
	{
		final double[] partials = new double[values.length];
		int size = 0;
		for (double x : values)
		{
			int i = 0;
			for (int j = 0; j < size; ++j) // size not necessarily == partials.length
			{
				double y = partials[j];

				if (abs(x) < abs(y))
				{
					final double tmp = x;
					x = y;
					y = tmp;
				}
				double hi = x + y;
				double lo = y - (hi - x);
				if (lo != 0.0)
					partials[i++] = lo;
				x = hi;
			}
			if (i < size)
			{
				partials[i] = x;
				Arrays.fill(partials,i+1,size,0);
			}
			else
			{
				partials[size++] = x;
			}
		}
		double sum = 0;
		for (double d : partials) sum += d;
		return sum;
	}

	// debugger function, not really needed
	private static long time = System.currentTimeMillis();
	static void logTime(String msg)
	{
		final long new_time = System.currentTimeMillis();
		System.out.println(msg+" "+(new_time-time));
		MathUtils.time = new_time;
	}

	/**
	 * returns newly allocated array of length len
	 * whose elements are doubles uniformaly distributed between 0.0 and 1.0
	 * generated by {@link MathUtils#random()}
	 * @param len
	 * @return
	 * Side Effects: none
	 */
	public static final double[] random(final int len)
	{
		final double ret[] = new double[len];
		for (int i = 0; i < len; i++)
			ret[i] = random();
		return ret;
	}

	static final private Random random = new Random();
	// not guaranteed to be thread safe
	// returns a double uniformly distributed between 0.0 and 1.0.
	// The implementation may change without warning, currently uses the Mersenne Twister.
	public static final double random()
	{
		return random.nextDouble();
	};

	// returns newly allocated array with same length as input
	// elements are max(input,0)
	final public static double[] positivePart(final double[] data)
	{
		final double [] ret = new double[data.length];
		for (int i = data.length; i--!=0;)
			ret[i] = max(data[i],0.0);
		return ret;
	}

	// returns newly allocated array
	// elements are min(input,-0)
	final public static double [] negativePart(final double[] data)
	{
		final double [] ret = new double[data.length];
		for (int i = data.length; 0!=i--;)
			ret[i] = min(data[i],-0.0);
		return ret;
	}

	// for comparison with apache commons. not for production code!
	/*
	static final double apacheSum(final double[] values)
	{
		return new org.apache.commons.math3.stat.descriptive.summary.Sum().evaluate(values);
	}
	*/

	/**
	 * Copies the sign of sign into magnitude.
	 * Faster than OpenJDK implementation.
	 * @param magnitude
	 * @param sign
	 * @return
	 */
	public static double copySign(final double magnitude, final double sign)
	{
		final long m = Double.doubleToLongBits(magnitude);
		final long s = Double.doubleToLongBits(sign);
		if (0<=(m^s)) 
			return -magnitude;
		return magnitude; // flip sign
	}
	/**
	 * Returns the absolute value of d (without branching).
	 * Faster than OpenJDK implementation.
	 * @param d
	 * @return
	 */
	public static final double abs(final double d)
	{
		return Double.longBitsToDouble(Long.MAX_VALUE & Double.doubleToRawLongBits(d));
	}
	/**
	 * Returns the absolute value of f (without branching).
	 * Faster than OpenJDK implementation.
	 * @param f
	 * @return
	 */
	public static final float abs(final float f)
	{
		return Float.intBitsToFloat(Integer.MAX_VALUE & Float.floatToRawIntBits(f));
	}
	/**
	 * Returns the absolute value of l (without branching).
	 * Faster than OpenJDK implementation.
	 * @param l
	 * @return
	 */
	public static final long abs(final long l)
	{
		final long sign = l>>>63;
		return (l^(~sign+1)) + sign;
	}
	/**
	 * Returns the absolute value of i (without branching).
	 * Faster than OpenJDK implementation.
	 * @param d
	 * @return
	 */
	public static final int abs(final int i)
	{
		final int sign = i>>>31;
		return (i^(~sign+1)) + sign;
	}

	/**
	 * Implementation of sum which is both more numerically
	 * stable _and faster_ than the naive implementation
	 * which is used in all standard numerical libraries I know of:
	 * Colt, OpenGamma, Apache Commons Math, EJML.
	 *
	 * Implementation uses variant of Kahan's algorithm keeping a running error
	 * along with the accumulator to try to cancel out the error at the end.
	 * This is much faster than Schewchuk's algorithm but not
	 * guaranteed to be perfectly precise
	 * In most cases, however, it is just as precise.
	 * Due to optimization it is about 30% faster
	 * even than the naive implementation on my machine.
	 * @param values
	 * @return
	 * Side Effects: none
	 */
	public static final double sum(final double... values)
	{
		double sum = 0;
		double err = 0;
		final int unroll = 6; // empirically it doesn't get much better than this
		final int len = values.length - values.length%unroll;

		// unroll the loop. due to IEEE 754 restrictions
		// the JIT shouldn't be allowed to unroll it dynamically, so it's
		// up to us to do it by hand ;)
		int i = 0;
		for (; i < len; i+=unroll)
		{
			final double val = values[i] + values[i+1]
				+ values[i+2] + values[i+3]
				+ values[i+4] + values[i+5];
			final double partial = val - err;
			final double hi = sum + val;
			err = (hi - sum) - partial;
			sum = hi;
		}
		for (; i < values.length; i++)
		{
			final double val = values[i];
			final double partial = val - err;
			final double hi = sum + val;
			err = (hi - sum) - partial;
			sum = hi;
		}
		return sum;
	}

	// Numerically naive, unoptimized version of sum. Intended for demonstrating superiority of
	// other methods only.
	public static final double sumNaive(final double ... values)
	{
		double ret = 0;
		for (final double x : values)
			ret += x;
		return ret;
	}

	/**
	 * Numerically naive implementation of sum which is faster than MathUtils.sum() and sumNaive()
	 * Generally exhibits rounding error which grows with the length of the sum
	 * Note that it may not agree with other implementations
	 * due to optimizations which change the order of iteration
	 * which can affect the rounding error.
	 * It is O(n) in the length of the array to be summed.
	 * It is faster than the naive, unoptimized implementation by 20-40%
	 * (dependent on the mood of the JIT) on my machine.
	 * @param values
	 * @return
	 * Side Effects: none
	 */
	public static final double sumFast(final double... values)
	{
		double ret = 0;
		// unroll the loop since the JIT shouldn't
		final int unroll = 4; // empirically unrolling more than 3 doesn't help much
		final int len = values.length - values.length%unroll;
		int i = 0;
		for (; i < len; i+=unroll)
			ret += values[i] + values[i+1] + values[i+2] + values[i+3];
		for (; i < values.length; i++)
			ret += values[i];
		return ret;
	}

	// Numerically precise implementation of sum.
	// Uses java.math.BigDecimal to internally keep an arbitrary
	// precision accumulator and then truncates at the end of the sum.
	// MathUtils.sumSlow(double[]), in addition to being about 10 times faster,
	// should (theoretically) return the same value as this method,
	// so this method shouldn't be used except as a sanity check.
	public static final double sumArbitraryPrecision(final double... values)
	{
		BigDecimal sum = new BigDecimal(0);
		for (int i = values.length; i-- != 0; )
			sum = sum.add(new BigDecimal(values[i]),MathContext.UNLIMITED);
		return sum.doubleValue();
	}

	static final boolean testSum (final double[] values)
	{
		return sum(values)==sumSlow(values);
	}

	/**
	 * Prints out the array in the format "[1.0, 2.0, 3.0]"
	 * @param d
	 */
	public final static void printArray(double[] d)
	{
		System.out.print("[");
		for (int i = 0; i < d.length; ++i)
		{
			if (0!=i) System.out.print(",");
			System.out.print(d[i]);
		}
		System.out.println("]");
	}

	static final void compare(final String msg, final double d1, final double d2)
	{
		final double diff = d1-d2;
		System.out.println("Diff "+msg+":"+diff);
		System.out.println("Exponent:"+Math.getExponent(diff));
		System.out.println("Mantissa:"+mantissa(diff));
		System.out.println("Precision:"+precision(diff));
	}

	/**
	 *  Returns the number of bits required to represent the mantissa of d
	 * by counting the number of trailing zeros in the mantissa.
	 * @param d
	 * @return
	 */
	public static final int precision(final double d)
	{
		final long l = Double.doubleToLongBits(d);
		return Math.max(0,(53 - Long.numberOfTrailingZeros(l)));
	}

	// Implementation which uses a BigDecimal for the accumulator
	// and so should have infinite precision.
	// Used for sanity checking.
	static final double linearCombinationArbitraryPrecision(final double []x, final double[]y)
	{
		final int len = Math.min(x.length,y.length);
		//final double [][] ret = new double[len][2];
		BigDecimal ret = new BigDecimal(0);
		for (int i = len; 0!=i--;)
		{
			BigDecimal product = new BigDecimal(x[i]).multiply(new BigDecimal(y[i]),MathContext.UNLIMITED);
			ret = ret.add(product,MathContext.UNLIMITED);
		}
		return ret.doubleValue();
	}

	/**
	 * Numerically precise dot product. Keeps a running error along with the
	 * accumulator. Equivalent to MathUtils.sum(MathUtils.multiply(x,y))
	 * but much faster and with O(1) memory overhead.
	 * O(n) with O(1) space.
	 * Even faster than the naive implementation ;).
	 * @param x Reference to first array
	 * @param y Reference to second array
	 * @param startOne Offset from beginning of first array
	 * @param startTwo Offset from beginning of second array
	 * @param len Number of terms to combine
	 * @throws ArrayIndexOutOfBoundsException if the input indices don't make sense
	 * (in particular, if startOne + len > x.length || startTwo + len > y.length)
	 * @return ddot(x,y)
	 * Side Effects: none
	 */
	public static final double linearCombination(final double[]x, final double[]y, final int startOne, final int startTwo, final int len)
	{
		//if (true) return MathUtils.sum(MathUtils.multiply(x,y));
		if (startOne + len > x.length || startTwo + len > y.length)
			throw new ArrayIndexOutOfBoundsException("Vector indices don't make sense!");
		final int unroll = 4; // don't blindly change this without changing the loop!
		// unroll was tuned to my machine. the optimal value is
		// probably architecture specific. one day java is give access to SIMD
		// instructions and then this can be optimized more.
		final int len_down = len - len%unroll;
		double sum = 0;
		double err = 0;
		int i = 0;
		int xPtr = startOne;
		int yPtr = startTwo;
		// hopefully the cpu will pipeline all of these things
		for (; i < len_down; i+= unroll,xPtr+=unroll,yPtr+=unroll)
		{
			// this line depends on unroll variable.
			final double prod = x[xPtr]*y[yPtr]
					+ x[xPtr+1]*y[yPtr+1]
					+ x[xPtr+2]*y[yPtr+2]
					+ x[xPtr+3]*y[yPtr+3];
			final double partial = prod - err;
			final double hi = sum + prod;
			err = (hi - sum) - partial;
			sum = hi;
		}
		for (; i < len; i++,xPtr++,yPtr++)
		{
			final double prod = x[xPtr]*y[yPtr];
			final double partial = prod - err;
			final double hi = sum + prod;
			err = (hi - sum) - partial;
			sum = hi;
		}
		return sum;
	}

	/**
	 * Calls {@link MathUtils#linearCombination(x,y,0,0,x.length)}
	 * @param x
	 * @param y
	 * @throws ArrayIndexOutOfBoundsException if x and y have unequal length
	 * @return
	 */
	public static final double linearCombination(final double[]x, final double[] y)
	{
		if (x.length!=y.length)
			throw new ArrayIndexOutOfBoundsException("Unequal length vectors!");
		return linearCombination(x,y,0,0,x.length);
	}

	/**
	 * Returns x*x, or x_1^2 + x_2^2 .., the Euclidean norm squared.
	 * @param x
	 * @return
	 */
	public static final double normSquared(final double[]x)
	{
		return linearCombination(x,x);
	}

	/**
	 * Returns sqrt(x*x), the Euclidean norm.
	 * @param x
	 * @return
	 */
	public static final double norm(final double[] x)
	{
		return Math.sqrt(normSquared(x));
	}

	// implementation of Strassen's algorithm.
	// not fast enough yet.
	private static final DoubleMatrix2D matrixMultiplyStrassen(final DoubleMatrix2D first, final DoubleMatrix2D second)
	{
		if (first.rows()!=first.cols())
			return null;
		if (second.rows()!= first.rows() || second.cols()!=second.rows())
			return null;
		if (!isPowerOfTwo(first.rows()))
			return null;
		final DoubleMatrix2D ret = new DoubleMatrix2D(first.rows(),first.rows());
		strassen(first, second, ret);
		return ret;
	}
	private static final int STRASSEN_IS_SLOWER = 256; // this should be tuned for different computers
	private static final void strassen(final DoubleMatrix2D first, final DoubleMatrix2D second, final DoubleMatrix2D ret)
	{
		final int n = first.rows();
		if (1==n) // 1x1
		{
			ret.set(0,0,first.get(0,0)*second.get(0,0));
			return;
		}
		// TODO case 2==n

		// see wikipedia for explanation of variable names
		final int m = n/2;
		final DoubleMatrix2D a11,a12,a21,a22,b11,b12,b21,b22, m1,m2,m3,m4,m5,m6,m7;
		a11 = first.submatrix(0,0,m,m);
		a12 = first.submatrix(0,m,m,m);
		a21 = first.submatrix(m,0,m,m);
		a22 = first.submatrix(m,m,m,m);
		b11 = second.submatrix(0,0,m,m);
		b12 = second.submatrix(0,m,m,m);
		b21 = second.submatrix(m,0,m,m);
		b22 = second.submatrix(m,m,m,m);
		m1 = new DoubleMatrix2D(m,m); m2 = new DoubleMatrix2D(m,m);
		m3 = new DoubleMatrix2D(m,m); m4 = new DoubleMatrix2D(m,m);
		m5 = new DoubleMatrix2D(m,m); m6 = new DoubleMatrix2D(m,m);
		m7 = new DoubleMatrix2D(m,m);
		if (m > STRASSEN_IS_SLOWER) // otherwise strassen is faster
		//if (false)
		//if (true)
		{
			DoubleMatrix2D slot1 = new DoubleMatrix2D(m,m);
			DoubleMatrix2D slot2 = new DoubleMatrix2D(m,m);
			// if we aren't on the first iteration anymore it is safe to overwrite the data in the arguments
			// (a11+a22)(b11+b22)
			strassen(DoubleMatrix2D.add(a11,a22,slot1), DoubleMatrix2D.add(b11,b22,slot2), m1);
			// (a21+a22)(b11)
			strassen(DoubleMatrix2D.add(a21,a22,slot1), b11, m2);
			// (a11)(b12-b22)
			strassen(a11,DoubleMatrix2D.subtract(b12,b22,slot2), m3);
			// (a22)(b21-b11)
			strassen(a22,DoubleMatrix2D.subtract(b21,b11,slot2), m4);
			// (a11+a12)(b22)
			strassen(DoubleMatrix2D.add(a11,a12,slot1),b22, m5);
			// (a21-a11)(b11+b12)
			strassen(DoubleMatrix2D.subtract(a21,a11,slot1),DoubleMatrix2D.add(b11,b12,slot2), m6);
			// (a12-a22)(b21+b22)
			strassen(DoubleMatrix2D.subtract(a12,a22,slot1),DoubleMatrix2D.add(b21,b22,slot2), m7);
		}
		else
		{
			DoubleMatrix2D slot1 = new DoubleMatrix2D(m,m);
			DoubleMatrix2D slot2 = new DoubleMatrix2D(m,m);
			// if we aren't on the first iteration anymore it is safe to overwrite the data in the arguments
			// (a11+a22)(b11+b22)
			matrixMultiplyFast(DoubleMatrix2D.add(a11,a22,slot1), DoubleMatrix2D.add(b11,b22,slot2), m1, false);
			// (a21+a22)(b11)
			matrixMultiplyFast(DoubleMatrix2D.add(a21,a22,slot1), b11, m2, false);
			// (a11)(b12-b22)
			matrixMultiplyFast(a11,DoubleMatrix2D.subtract(b12,b22,slot2), m3, false);
			// (a22)(b21-b11)
			matrixMultiplyFast(a22,DoubleMatrix2D.subtract(b21,b11,slot2), m4, false);
			// (a11+a12)(b22)
			matrixMultiplyFast(DoubleMatrix2D.add(a11,a12,slot1),b22, m5, false);
			// (a21-a11)(b11+b12)
			matrixMultiplyFast(DoubleMatrix2D.subtract(a21,a11,slot1),DoubleMatrix2D.add(b11,b12,slot2), m6, false);
			// (a12-a22)(b21+b22)
			matrixMultiplyFast(DoubleMatrix2D.subtract(a12,a22,slot1),DoubleMatrix2D.add(b21,b22,slot2), m7, false);
		}

		int idx = 0;
		int idx2 = ret.index(0,0);
		for (int i = 0; i < m; i++)
		{
			for (int j = 0; j < m; j++,idx++,idx2++) // try to help it pipeline the instructions
			{
				// fracking pointer arithmetic
				// TODO refactor to use DoubleMatrix2D.ptr()
				final int idx11 = idx2;
				final int idx12 = idx2 + m;
				final int idx21 = idx2 + m*n;
				final int idx22 = idx2 + m*n + m;
				ret.data()[idx11] = m1.data()[idx] + m4.data()[idx] - m5.data()[idx] + m7.data()[idx]; // C11
				ret.data()[idx12] = m3.data()[idx] + m5.data()[idx]; // C12
				ret.data()[idx21] = m2.data()[idx] + m4.data()[idx]; // C21
				ret.data()[idx22] = m1.data()[idx] - m2.data()[idx] + m3.data()[idx] + m6.data()[idx]; // C22
			}
			// since idx2 is going to be m more from where it was at the beginning of this iteration
			// and we want it to be n more.
			idx2+=m;
		}
		return ;
	}

	static final boolean isPowerOfTwo(final int i)
	{
		return (i & (~i + 1)) == i;
	}

	// TODO refactor to use DoubleMatrix2D.ptr()
	/**
	 * Returns the numerically precise result of multiplying matrix 1 by matrix 2.
	 * Faster than any other Java matrix multiplication I know of.
	 * @param first
	 * @param second
	 * @return
	 * Side Effects: Allocation of new matrix
	 */
	public static final DoubleMatrix2D matrixMultiply(final DoubleMatrix2D first, final DoubleMatrix2D second)
	{
		if (first.cols()!=second.rows())
			throw new ArrayIndexOutOfBoundsException("Trying to multiply matrices of different dimensions?!");
		final DoubleMatrix2D ret = new DoubleMatrix2D(first.rows(), second.cols());
		for (int i = 0; i < second.cols(); i++)
		{
			// extract the column beforehand to improve cache hits
			// TODO think about extracting on fly to save cost of read
			final double vector[] = second.extractColumn(i);
			int idx = ret.index(0,i);
			for (int j = 0; j < first.rows(); j++,idx+=ret.real_row_size())
			{
				final double val = linearCombination(first.data(), vector, first.index(j,0), 0, first.cols());
				ret.data()[idx] = val;
			}
		}
		return ret;
	}

	/**
	 * Fast matrix multiply, at the expense of numerical precision.
	 * @param first
	 * @param second
	 * @return
	 * Side Effects: Allocation of new matrix
	 */
	public static final DoubleMatrix2D matrixMultiplyFast(final DoubleMatrix2D first, final DoubleMatrix2D second)
	{
		final DoubleMatrix2D ret = new DoubleMatrix2D(first.rows(), second.cols());
		return matrixMultiplyFast(first,second,ret,true);
	}

	// TODO refactor to use DoubleMatrix2D.ptr()
	/**
	 * Fast matrix multiply, at the expense of numerical precision.
	 * This method will overwrite the contents of ret and return it.
	 * @param first
	 * @param second
	 * @param ret
	 * @param sanity_check Perform array index sanity checks and whatnot
	 * @return ret
	 * Side Effects: Overwriting of contents of ret.
	 */
	public static final DoubleMatrix2D matrixMultiplyFast(final DoubleMatrix2D first, final DoubleMatrix2D second, final DoubleMatrix2D ret, final boolean sanity_check)
	{
		if (sanity_check)
		{
			if (first.cols()!=second.rows())
				throw new ArrayIndexOutOfBoundsException("Trying to multiply matrices of different dimensions?!");
			if (ret.rows()!=first.rows()||ret.cols()!=second.cols())
				throw new IllegalArgumentException("Bad ret");
		}
		int first_idx = first.index(0,0);
		int ret_idx = ret.index(0,0);
		for (int i = 0; i < second.cols(); i++,ret_idx-=ret.rows()*ret.real_row_size()-1,first_idx-=first.rows()*first.real_row_size())
		{
			// extract the column beforehand to improve cache hits
			final double vector[] = second.extractColumn(i);
			// unroll the loop, since the JIT doesn't seem to want to do so.
			// in general i am guessing that if there is a cache miss on write
			// that is because the matrix is so big that the linear combination
			// will fill up the cache with other stuff anyways
			int j = 0;
			for (; j < first.rows()/2; j++,ret_idx+=ret.real_row_size()<<1,first_idx+=first.real_row_size()<<1)
			{
				// if (first_idx!=first.index(j,0)) throw new RuntimeException("You have bug!");
				// if (ret_idx!=ret.index(j,i)) throw new RuntimeException("You have bug!");
				ret.data()[ret_idx]
					= linearCombinationFast(first.data(), vector, first_idx, 0, first.cols());
				ret.data()[ret_idx + ret.real_row_size()]
					= linearCombinationFast(first.data(),vector,first_idx+first.real_row_size(),0,first.cols());
			}
			// extra stuff
			if (0!=(first.rows()&1))
			{
				ret.data()[ret_idx]
					= linearCombinationFast(first.data(), vector, first_idx, 0, first.cols());
			}
		}
		return ret;
	}

	/**
	 * Returns the i'th column of matrix as an array.
	 * @param matrix
	 * @param col
	 * @return
	 */
	public static final double[] extractColumn(final double[][] matrix, final int col)
	{
		final double[] ret = new double[matrix[0].length];
		for (int i = 0; i < matrix.length; i++)
			ret[i] = matrix[i][col];
		return ret;
	}

	// TODO think about parallelizing this
	public static final double[][] matrixMultiply(final double[][] first, final double[][] second)
	{
		// i,j,k
		final int firstRows = first.length;
		final int firstCols = first[0].length;
		final int secondRows = second.length;
		final int secondCols = second[0].length;
		if (firstCols!=secondRows)
			throw new ArrayIndexOutOfBoundsException("Trying to multiply matrices of different dimensions?!");
		final double ret[][] = new double[firstRows][secondCols];
		for (int i = 0; i < secondCols; i++)
		// iterate over columns so we can maintain cache locality!
		{
			// TODO think about extracting column on the fly
			final double[] vector = extractColumn(second,i);
			for (int k = 0; k < firstRows; k++)
			{
				ret[k][i] = linearCombination(first[k],vector);
			}
		}
		return ret;
	}

	// attempt to be faster than matrixMultiply. failure.
	public static final double[][] matrixMultiplyTwo(final double[][] first, final double[][] second)
	{
		final int firstRows = first.length;
		final int firstCols = first[0].length;
		final int secondRows = second.length;
		final int secondCols = second[0].length;
		if (firstCols!=secondRows)
			throw new ArrayIndexOutOfBoundsException("Trying to multiply matrices of different dimensions?!");
		final double ret[][] = new double[firstRows][secondCols];
		final double err[][] = new double[firstRows][secondCols];

		final int unroll = 3;
		final int len = secondCols - secondCols%unroll;

		for (int j = 0; j < firstRows; j++)
		{
			for (int i = 0; i < secondCols; i++)
			{
				Arrays.fill(ret[j],0);
				Arrays.fill(err[j],0);
				for (int k = 0; k < firstCols; k++)
				{
					final double prod = first[j][k] * second[k][i];
					final double sum = ret[j][i];
					final double hi = sum + prod;
					ret[j][i] = hi;
					err[j][i] += hi - sum - prod;
				}
			}
		}

		for (int i = 0; i < firstRows; i++)
			for (int j = 0; j < secondCols; j++)
				ret[i][j] -= err[i][j];
		return ret;
	}

	public static final double[][] matrixMultiplyFast(final double[][] first, final double[][] second)
	{
		final int firstRows = first.length;
		final int firstCols = first[0].length;
		final int secondRows = second.length;
		final int secondCols = second[0].length;
		if (firstCols!=secondRows)
			throw new ArrayIndexOutOfBoundsException("Trying to multiply matrices of different dimensions?!");
		final double ret[][] = new double[firstRows][secondCols];
		for (int i = 0; i < secondCols; i++)
		// iterate over columns so we can maintain cache locality!
		{
			final double[] vector = extractColumn(second,i);
			for (int k = 0; k < firstRows; k++)
			{
				ret[k][i] = linearCombinationFast(first[k],vector);
			}
		}
		return ret;
	}

	// Pre: matrix has m rows and n columns, and vector has n elements.
	// Note: no checking on the sizes of the inputs!
	// May throw ArrayIndexOutOfBoundsExceptions or other such
	// nasty things if you don't sanitize input!
	public static final double[] matrixMultiply(final double[][] matrix, final double[] vector)
	{
		// recall matrix is double[rows][cols], and matrix.length==rows
		final double ret[] = new double[matrix.length];
		for (int i = 0; i < matrix.length; i++)
			ret[i] = linearCombination(matrix[i],vector);
		return ret;
	}

	public static final double[] matrixMultiplyFast(final double[][] matrix, final double[] vector)
	{
		final double []ret = new double[matrix.length];
		for (int i = 0; i < matrix.length; ++i)
			ret[i] = linearCombinationFast(matrix[i],vector);
		return ret;
	}

	public static final double[] matrixMultiplyFast(final DoubleMatrix2D matrix, final double[] vector)
	{
		final double []ret = new double[matrix.rows()];
		for (int i = 0; i < matrix.rows(); i++)
		{
			ret[i] = linearCombinationFast(matrix.data(),vector,matrix.index(i,0),0,matrix.cols());
		}
		return ret;
	}

	public static final double[] matrixMultiply(final DoubleMatrix2D matrix, final double[] vector)
	{
		final double ret[] = new double[matrix.rows()];
		for (int i = 0; i < matrix.rows(); i++)
		{
			ret[i] = linearCombination(matrix.data(),vector,matrix.index(i,0),0,matrix.cols());
		}
		return ret;
	}

	// Numerically precise dot product. Returns MathUtils.sumSlow(MathUtils.multiply(x,y));
	// O(n) time and O(n) space.
	// TODO make if faster by not allocating new array in multiply().
	static final double linearCombinationSlow(final double[]x,final double[]y)
	{
		return sumSlow(multiply(x,y));
	}

	// Faster but lower precision than linearCombination.
	// Numerically naive implementation.
	public static final double linearCombinationFast(final double[]x, final double[]y)
	{
		if (x.length!=y.length)
			throw new ArrayIndexOutOfBoundsException("Dot product of vectors with different lengths!");
		return linearCombinationFast(x,y,0,0,x.length);
	}
	public static final double linearCombinationFast(final double []x, final double[]y, final int startOne, final int startTwo, final int len)
	{
		if (startOne + len > x.length || startTwo + len > y.length)
			throw new ArrayIndexOutOfBoundsException("Bad length!");
		double ret = 0;
		final int unroll = 4;
		final int len_down = len - len%unroll;
		int i = 0;
		int xPtr = startOne;
		int yPtr = startTwo;
		for (; i < len_down; i+=unroll,xPtr+=unroll,yPtr+=unroll)
		{
			ret+= x[xPtr]*y[yPtr]
				+ x[xPtr+1]*y[yPtr+1]
				+ x[xPtr+2]*y[yPtr+2]
				+ x[xPtr+3]*y[yPtr+3];
		}
		// get the terms at the end
		for (; i < len; i++,xPtr++,yPtr++)
		{
			ret += x[xPtr]*y[yPtr];
		}
		return ret;
	}

	static final double distance(final double[]p, final double[]q)
	{
		// TODO make this faster if needed, as it is it is going to loop like three times ;)
		return Math.sqrt(sum(pow(subtract(p,q),2)));
	}

	/**
	 * Returns a double with the exponent and sign set to zero.
	 * @param d
	 * @return
	 */
	static final double mantissa(final double d)
	{
		return abs(Math.scalb(d,-Math.getExponent(d)));
	}

	/**
	 * Returns new array which is reverse of the argument.
	 * @param data
	 * @return
	 * Side Effects: none
	 */
	public static final double[] reverse(final double[] data)
	{
		final double[] ret = new double[data.length];
		int center = data.length / 2;
		while (center--!=0)
		{
			final int left = center;
			final int right = data.length-center-1;
			ret[left] = data[right];
			ret[right] = data[left];
		}
		return ret;
	}

	/**
	 * Allocates new array which is sorted version of argument.
	 * O(n) in the allocation and then O(n log n) in the sort.
	 * Behavior should be identical to calling Arrays.sort(data.clone())
	 * @param data
	 * @return
	 * Side Effects: none
	 */
	public static final double[] sort(final double[] data)
	{
		final double ret[] = copyOf(data);
		Arrays.sort(ret);
		return ret;
	}

	/**
	 * Behavior is identical to calling data.clone() or Arrays.copyOf(data)
	 * But can be up to 30% faster if the JIT doesn't optimize those functions
	 * @param data
	 * @return
	 */
	public static final double[] copyOf(final double[] data)
	{
		return copyOfRange(data,0,data.length);
	}

	/**
	 * Behavior is identical to calling data.clone() or Arrays.copyOf(data)
	 * But can be up to 30% faster if the JIT doesn't optimize those functions
	 * @param data
	 * @return
	 */
	public static final double[] copyOfRange(final double[]data, final int start, final int end)
	{
		if (end > data.length || start < 0)
			throw new IllegalArgumentException("Bad array bounds!");
		final double ret[] = new double[end-start];
		for (int i = 0; i < (end-start)/3; i++)
		{
			final int x = i*3;
			ret[x] = data[x+start];
			ret[x+1] = data[x+1+start];
			ret[x+2] = data[x+2+start];
		}
		// Don't care about extra copies if data.length%3==0
		ret[ret.length-2] = data[end-2];
		ret[ret.length-1] = data[end-1];
		return ret;
	}

	/**
	 * Returns true if two arrays have same length
	 * and all the elements of the two arrays are equal
	 * @param x
	 * @param y
	 * @throws NullPointerException if either x or y are null.
	 * @return
	 */
	public static final boolean equals(final double[]x, final double[]y)
	{
		final int len = x.length;
		if (y.length!=len)
			return false;
		for (int i = len; i--!=0;)
			if (x[i]!=y[i]) // TODO deal with NaNs?
				return false;
		return true;
	}

	public static final double[][] covariance(final double[][] data)
	{
		final int len = data.length;
		final double []means = new double[len]; // precalculate the means
		final double[][] ret = new double[len][len];
		for (int i = 0; i < len; i++)
		{
			means[i] = mean(data[i]);
			for (int j = 0; j <= i; j++)
			{
				final double d = sum(multiply(shift(data[i],-means[i]),shift(data[j],-means[j]))) / (len);
				ret[i][j] = d;
				ret[j][i] = d;
			}
		}
		return ret;
	}

	final static boolean fuzzyEquals(final double d, final double e)
	{
		return (float)d==(float)e;
	}
	final static boolean fuzzyEquals(final double []x, final double []y)
	{
		if (x.length!=y.length)
			return false;
		for (int i = 0; i < x.length; i++)
			if (!fuzzyEquals(x[i],y[i]))
				return false;
		return true;
	}
	final static boolean fuzzyEquals(final DoubleMatrix2D a, final DoubleMatrix2D b)
	{
		if (a.cols()!=b.cols()||a.rows()!=b.rows())
			return false;
		for (int i = 0; i < a.rows(); i++)
			for (int j = 0; j < a.rows(); j++)
				if (!fuzzyEquals(a.get(i,j),b.get(i,j)))
					return false;
		return true;
	}

	/**
	 * Test that matrixMultiplyStrassen gives same results as matrixMultiply
	 * @param one
	 * @param two
	 */
	final static void testMatrixMultiply(final DoubleMatrix2D one, final DoubleMatrix2D two)
	{
		final DoubleMatrix2D resultOne = matrixMultiplyStrassen(one,two);
		final DoubleMatrix2D resultTwo = matrixMultiply(one,two);
		for (int i = 0; i < resultOne.data().length; i++)
			if ( (float)resultOne.data()[i] != (float)resultTwo.data()[i])
				throw new RuntimeException("i:"+i);
	}
	/**
	 * Tests matrixMultiply by checking
	 * that both left and right multiplication
	 * of a randomly generated two by two matrix with the identity
	 * returns the same matrix.
	 */
	final static void testMatrixMultiply()
	{
		final DoubleMatrix2D a = DoubleMatrix2D.identity(2);
		final DoubleMatrix2D b = DoubleMatrix2D.random(2,2);
		if (!fuzzyEquals(matrixMultiply(a,b),b))
			throw new RuntimeException("Failed matrixMultiply");
		if (!fuzzyEquals(matrixMultiply(b,a),b))
			throw new RuntimeException("Failed matrixMultiply");
		System.out.println("passed matrixMultiply");
	}
	/**
	 * Tests matrixMultiplyFast by checking
	 * that both left and right multiplication
	 * of a randomly generated two by two matrix with the identity
	 * returns the same matrix.
	 */
	final static void testMatrixMultiplyFast()
	{
		final DoubleMatrix2D a = DoubleMatrix2D.identity(2);
		final DoubleMatrix2D b = DoubleMatrix2D.random(2,2);
		if (!fuzzyEquals(matrixMultiplyFast(a,b),b))
			throw new RuntimeException("Failed matrixMultiplyFast");
		if (!fuzzyEquals(matrixMultiplyFast(b,a),b))
			throw new RuntimeException("Failed matrixMultiplyFast");
		System.out.println("passed matrixMultiplyFast");
	}
	/**
	 * Tests matrixMultiplyStrassen by checking
	 * that both left and right multiplication
	 * of a randomly generated two by two matrix with the identity
	 * returns the same matrix.
	 */
	final static void testMatrixMultiplyStrassen()
	{
		final DoubleMatrix2D a = DoubleMatrix2D.identity(2);
		final DoubleMatrix2D b = DoubleMatrix2D.random(2,2);
		if (!fuzzyEquals(matrixMultiplyStrassen(a,b),b))
			throw new RuntimeException("Failed matrixMultiplyStrassen");
		if (!fuzzyEquals(matrixMultiplyStrassen(b,a),b))
			throw new RuntimeException("Failed matrixMultiplyStrassen");
		System.out.println("passed matrixMultiplyStrassen");
	}
	/**
	 * Tests linearCombination and linearCombinationFast
	 * by checking that [1,1,1,1] * [1,1,1,1] == 4.
	 */
	final static void testLinearCombination()
	{
		final double[] one = new double[]{1,1,1,1};
		final double[] two = new double[]{1,1,1,1};
		if (4.0!=linearCombination(one,two))
			throw new RuntimeException("linearCombination Test Failed!");
		if (4.0!=linearCombinationFast(one,two))
			throw new RuntimeException("linearCombinationFast Test Failed!");
		System.out.println("passed linearCombination");
	}
	/**
	 * Test addition and subtraction of matrices
	 */
	final static void testAddSubtract()
	{
		DoubleMatrix2D test = DoubleMatrix2D.random(32,32);
		DoubleMatrix2D result = test.plus(test).minus(test);
		if (!fuzzyEquals(test,result))
			throw new RuntimeException("addSubtract Test Failed!");
		DoubleMatrix2D test1 = test.submatrix(0, 0, 16, 16);
		DoubleMatrix2D test2 = test.submatrix(16, 16, 16, 16);
		result = test1.plus(test2).minus(test2);
		if (!fuzzyEquals(test1,result))
			throw new RuntimeException("addSubtract Test Failed!");
		result = test2.plus(test1).minus(test1);
		if (!fuzzyEquals(test2,result))
			throw new RuntimeException("addSubtract Test Failed!");

		test = DoubleMatrix2D.identity(2);
		if (0.0!=test.submatrix(0,1,1,1).plus(test.submatrix(1,0,1,1)).get(0,0))
			throw new RuntimeException("addSubtractTest Failed!");
		if (0.0!=test.submatrix(0,1,1,1).minus(test.submatrix(1,0,1,1)).get(0,0))
			throw new RuntimeException("addSubtractTest Failed!");
		System.out.println("passed addSubtract");
	}
	public static final DoubleMatrix2D transpose(DoubleMatrix2D A){
		  DoubleMatrix2D AT = new DoubleMatrix2D(A.cols(),A.rows());
		  for(int i = 0; i < A.rows(); i++){
			  double temp[] = new double[A.rows()];
			  
			  //for(int k = 0; k < temp.length; k++){
			 
			    temp = A.extractColumn(i);
			//  }
			  for(int j = 0; j < temp.length; j++){
				  AT.set(temp[j]);
				  AT.incrementColumn();
			  }
		
		  }
		  AT.resetPtr();
		  return AT; 
	}

	public static void main(String[] args)
	{
		/*
		logTime("start");
		final int len = 1000*1000*100;
		final double[] x = random(len);
		final double[] y = random(len);
		for (int i = 0; i < 5; i++)
			linearCombinationFast(x,y);
		logTime("warmed up");
		linearCombinationFast(x,y);
		logTime("done");
		if (true) return;
		logTime("start");
		testAddSubtract();
		testLinearCombination();
		testMatrixMultiply();
		testMatrixMultiplyFast();
		testMatrixMultiplyStrassen();
		//final int len = 1024;
		final DoubleMatrix2D one = DoubleMatrix2D.random(len,len);
		final DoubleMatrix2D two = DoubleMatrix2D.random(len,len);
//		final DoubleMatrix2D one = DoubleMatrix2D.diagonal(new double[]{1,2,1,1});
//		final DoubleMatrix2D two = DoubleMatrix2D.diagonal(new double[]{2,1,1,1});
		DoubleMatrix2D foo = null;
		testMatrixMultiply(one,two); System.out.println("passed");
		logTime("one");
		for (int i = 0; i < 5; i++)
			foo = matrixMultiplyStrassen(one,two);
		logTime("Warmed up");
		for (int i = 0; i < 1000; i++)
			foo = matrixMultiplyStrassen(two,one);
		logTime("done");

		if (true) return;
		*/
		logTime("start");
		final int len = 1000*1000*10;
		final double[] data = shift(random(len),-.5);
		logTime("random");
		double fast,medium,slow;
		fast = medium = slow = 0;

		for (int i = 0; i < 100; i++)
		{
			fast = sumFast(data);
			slow = sum(data);
			medium = sumNaive(data);
		}
		compare("FS",fast,slow);
		compare("MS",medium,slow);
		compare("FM",fast,medium);

		logTime("warmup");
		for (int i = 0; i < 100; i++)
			sum(data);
		logTime("slow");
		for (int i = 0; i < 100; i++)
			sumFast(data);
		logTime("fast");
		for (int i = 0; i < 100; i++)
			sumNaive(data);
		logTime("mine");
	}
}

