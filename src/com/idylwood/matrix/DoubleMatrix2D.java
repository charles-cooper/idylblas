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

package com.idylwood.matrix;

import com.idylwood.utils.MathUtils;

public class DoubleMatrix2D
{
	private final double[] data;
	private final int rows;
	private final int cols;
	public DoubleMatrix2D(double[][] data)
	{
		this(data.length, data[0].length);
		for (int i = 0; i < rows; i++)
		{
			if (cols != data[i].length)
				throw new ArrayIndexOutOfBoundsException("Illegal matrix: uneven number of columns.");
			System.arraycopy(data[i], 0, this.data, i*cols, cols);
		}
	}
	public DoubleMatrix2D(final int rows, final int cols)
	{
		this.rows = rows;
		this.cols = cols;
		this.data = new double[rows*cols];
	}
	// square matrix construction. does not copy the array.
	private DoubleMatrix2D(final double data[], final int rows, final int cols)
	{
		this.rows = rows; this.cols = cols;
		if (this.rows*this.cols != data.length)
			throw new RuntimeException("Matrix must be square!");
		this.data = data;
	}

	public int rows() { return this.rows; }
	public double[] data() { return this.data; }
	public int cols() { return this.cols; }
	public static DoubleMatrix2D random(final int rows, final int cols)
	{
		final double[] data = MathUtils.random(rows*cols);
		return new DoubleMatrix2D(data,rows,cols);
	}
	public int index(final int row, final int col)
	{
		return row * cols + col;
	}
	public double get(final int row, final int col)
	{
		return data[index(row,col)];
	}
	public void set(final int row, final int col, final double val)
	{
		data[index(row,col)] = val;
	}
	public void increment(final int row, final int col, final double incr)
	{
		data[index(row,col)] += incr;
	}
	public double[] extractColumn(final int col)
	{
		final double[] ret = new double[rows];
		for (int i = 0; i < rows; i++)
		{
			ret[i] = get(i,col);
		}
		return ret;
	}
	public DoubleMatrix2D plus(final DoubleMatrix2D other)
	{
		return add(this,other);
	}
	public DoubleMatrix2D minus(final DoubleMatrix2D other)
	{
		return subtract(this,other);
	}
	public static DoubleMatrix2D add(final DoubleMatrix2D one, final DoubleMatrix2D two)
	{
		if (one.cols()!=two.cols() || one.rows()!=two.rows())
			throw new RuntimeException("uh oh");
		final double[] data = MathUtils.add(one.data(),two.data());
		final DoubleMatrix2D ret = new DoubleMatrix2D(data,one.rows(),one.cols());
		return ret;
	}
	public static DoubleMatrix2D subtract(final DoubleMatrix2D minuend, final DoubleMatrix2D subtrahend)
	{
		if (minuend.cols()!=subtrahend.cols() || minuend.rows()!=subtrahend.rows())
			throw new RuntimeException("uh oh");
		final double[] data = MathUtils.subtract(minuend.data(),subtrahend.data());
		final DoubleMatrix2D ret = new DoubleMatrix2D(data,minuend.rows(),minuend.cols());
		return ret;
	}
	public double[] extractRow(int row)
	{
		return MathUtils.copyOfRange(data,row*cols,(row+1)*cols);
	}
}

