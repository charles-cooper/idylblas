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
	DoubleMatrix2D(final double data[], final int rows, final int cols)
	{
		this.rows = rows; this.cols = cols;
//		if (this.rows*this.cols != data.length)
//			throw new RuntimeException("Matrix must be square!");
		this.data = data;
	}

	public int rows() { return this.rows; }
	public double[] data() { return this.data; }
	public int cols() { return this.cols; }
	public static final DoubleMatrix2D random(final int rows, final int cols)
	{
		final double[] data = MathUtils.random(rows*cols);
		return new DoubleMatrix2D(data,rows,cols);
	}
	public static final DoubleMatrix2D diagonal(final double []data)
	{
		final double[] input = new double[data.length*data.length];
		for (int i = 0; i < data.length; i++)
			input[i*data.length + i] = data[i];
		return new DoubleMatrix2D(input,data.length,data.length);
	}
	public DoubleMatrix2D submatrix(final int row_offset, final int col_offset, final int rows, final int cols)
	{
		return new SubDoubleMatrix2D(this, row_offset, col_offset, rows, cols);
	}
	/** Calculates the index of a row and column in the underlying data structure.
	 * For example for an ordinary matrix laid out in memory row by row,
	 * it will return row * this.cols() + col.
	 * @param row
	 * @param col
	 * @return The index of the row and column.
	 */
	public int index(final int row, final int col)
	{
		return row * cols + col;
	}
	/**
	 * Returns the entry at row, col.
	 * @param row
	 * @param col
	 * @return
	 */
	public final double get(final int row, final int col)
	{
		return data[index(row,col)];
	}
	public void set(final int row, final int col, final double val)
	{
		data[index(row,col)] = val;
	}
	/**
	 * Copies the other matrix into the block matrix starting from row,col.
	 * @param row
	 * @param col
	 * @param other
	 */
	public void set(final int row, final int col, final DoubleMatrix2D other)
	{
		for (int i = 0; i < other.rows();i++)
		{
			System.arraycopy(other.data(), i*other.cols(), this.data(), col + (row+i)*this.cols(), other.cols());
		}
	}
	/**
	 * Increments the value at <code>row</code>,<code>col</code> by <code>incr</code>.
	 * @param row
	 * @param col
	 * @param incr
	 */
	public void increment(final int row, final int col, final double incr)
	{
		data[index(row,col)] += incr;
	}
	public DoubleMatrix2D plus(final DoubleMatrix2D other)
	{
		return add(this,other);
	}
	public DoubleMatrix2D plus(final DoubleMatrix2D other, final int this_idx, final int other_idx, final int num_rows, final int num_cols)
	{
		final DoubleMatrix2D ret = new DoubleMatrix2D(num_rows,num_cols);
		for (int i = 0; i < num_rows; i++)
		{
//			int this_idx = this.index(this_i+i, this_j);
//			int other_idx = other.index(other_i+i, other_j);
			int ret_idx = ret.index(i, 0);
			for (int j = 0; j < num_cols; j++)
			{
				final double x = this.data[this_idx+i*this.cols+j];
				final double y = other.data[other_idx+i*other.cols+j];
				ret.data[ret_idx+j] = x+y;
			}
		}
		return ret;
	}
	public DoubleMatrix2D minus(final DoubleMatrix2D other, final int this_i, final int this_j, final int other_i, final int other_j, final int num_rows, final int num_cols)
	{
		final DoubleMatrix2D ret = new DoubleMatrix2D(num_rows,num_cols);
		for (int i = 0; i < num_rows; i++)
		{
			int this_idx = this.index(this_i+i, this_j);
			int other_idx = other.index(other_i+i, other_j);
			int ret_idx = ret.index(i, 0);
			for (int j = 0; j < num_cols; j++)
			{
				final double x = this.data[this_idx+j];
				final double y = other.data[other_idx+j];
				ret.data[ret_idx+j] = x-y;
			}
		}
		return ret;
	}

	public DoubleMatrix2D minus(final DoubleMatrix2D other)
	{
		return subtract(this,other);
	}
	public static DoubleMatrix2D add(final DoubleMatrix2D one, final DoubleMatrix2D two)
	{
		if (one.cols()!=two.cols())
			throw new RuntimeException("uh oh");
		if (one.rows()!=two.rows())
			throw new RuntimeException("uh oh");
		final DoubleMatrix2D ret = new DoubleMatrix2D(one.rows(),one.cols());
		int idx = 0;
		for (int i = 0; i < one.rows();i++)
			for (int j = 0; j < one.cols();j++,idx++)
				ret.data()[idx]= one.get(i,j)+two.get(i,j);
		return ret;
	}
	public static DoubleMatrix2D subtract(final DoubleMatrix2D minuend, final DoubleMatrix2D subtrahend)
	{
		if (minuend.cols()!=subtrahend.cols() || minuend.rows()!=subtrahend.rows())
			throw new RuntimeException("uh oh");
		final DoubleMatrix2D ret = new DoubleMatrix2D(minuend.rows(),minuend.cols());
		int idx = 0;
		for (int i = 0; i < minuend.rows();i++)
			for (int j = 0; j < minuend.cols();j++,idx++)
				ret.data()[idx] = minuend.get(i,j) - subtrahend.get(i,j);
		return ret;
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
	public double[] extractRow(int row)
	{
		final double[] ret = new double[cols];
		for (int i = 0; i < cols; i++)
		{
			ret[i] = get(row,i);
		}
		return ret;
		//return MathUtils.copyOfRange(data,row*cols,(row+1)*cols);
	}
}

