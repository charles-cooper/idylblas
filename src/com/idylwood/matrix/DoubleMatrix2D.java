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

/**
 * Class which represents a two dimensional matrix.
 * Currently the data is stored internally in a one dimensional array
 * in (C-like) row major format, but this is due to change.
 * It also provides an iterator-like interface through ptr()
 * and related methods.
 * @author charles
 *
 */
// TODO experiment with block recursive storage.
public class DoubleMatrix2D
{
	private final double[] data;
	private final int rows;
	private final int cols;
	private final int real_rows;
	private final int real_cols;
	private final int row_offset;
	private final int col_offset;
	private int ptr;
	public DoubleMatrix2D(final double[][] data)
	{
		this(data.length,data[0].length);
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
		this.real_rows = rows;
		this.real_cols = cols;
		this.row_offset = 0;
		this.col_offset = 0;
		this.data = new double[rows*cols];
	}
	// square matrix construction. does not copy the array.
	DoubleMatrix2D(final double data[],final int real_rows, final int real_cols, final int rows,final int cols,final int row_offset, final int col_offset)
	{
		this.rows = rows; this.cols = cols;
		this.real_rows = real_rows; this.real_cols = real_cols;
		this.row_offset = row_offset; this.col_offset = col_offset;
		/*
		 * TODO put this check back in, it should not be failing!
		if (this.real_rows*this.real_cols != data.length)
			throw new RuntimeException("Matrix must be square!");
			*/
		this.data = data;
	}
	DoubleMatrix2D(final double data[], final int rows, final int cols)
	{
		this(data,rows,cols,rows,cols,0,0);
	}

	public final int rows() { return this.rows; }
	public final double[] data() { return this.data; }
	public final int cols() { return this.cols; }
	public final int real_row_size() { return this.real_cols; }

	public static final DoubleMatrix2D identity(final int size)
	{
		final double data[] = new double[size];
		java.util.Arrays.fill(data,1);
		return diagonal(data);
	}
	/**
	 * Returns a newly allocated matrix whose entries are randomly generated
	 * numbers in [0,1).
	 * @param rows
	 * @param cols
	 * @return
	 */
	public static final DoubleMatrix2D random(final int rows, final int cols)
	{
		final double[] data = MathUtils.random(rows*cols);
		return new DoubleMatrix2D(data,rows,cols);
	}
	/**
	 * Returns a newly allocated diagonal matrix with entries given by <code>data</code>
	 * @param data
	 * @return
	 */
	public static final DoubleMatrix2D diagonal(final double []data)
	{
		final double[] input = new double[data.length*data.length];
		for (int i = 0; i < data.length; i++)
			input[i*data.length + i] = data[i];
		return new DoubleMatrix2D(input,data.length,data.length);
	}
	/**
	 * Returns newly allocated matrix with <code>rows</code> rows and <code>cols</code> columns
	 * where all the entries are equal to <code>val</code>
	 * @param rows
	 * @param cols
	 * @param val
	 * @return
	 */
	public static final DoubleMatrix2D constant(final int rows, final int cols, final double val)
	{
		final DoubleMatrix2D ret = new DoubleMatrix2D(rows,cols);
		java.util.Arrays.fill(ret.data, val);
		return ret;
	}
	/**
	 * Abstraction which returns a submatrix. THIS IS BACKED BY THE SAME MEMORY.
	 * ANY MODIFICATIONS TO THIS SUBMATRIX WILL MODIFY THE PARENT MATRIX.
	 * If you want a copy, then call DoubleMatrix2D.submatrix().copy().
	 * @param row_offset
	 * @param col_offset
	 * @param rows
	 * @param cols
	 * @return
	 */
	public final DoubleMatrix2D submatrix(final int row_offset, final int col_offset, final int rows, final int cols)
	{
		return new DoubleMatrix2D(this.data,this.real_rows,this.real_cols,rows,cols,this.row_offset+row_offset,this.col_offset+col_offset);
	}
	/**
	 * Returns a (newly allocated) copy of the matrix.
	 * @return
	 */
	public final DoubleMatrix2D copy()
	{
		final DoubleMatrix2D ret = new DoubleMatrix2D(this.rows,this.cols);
		int idx = 0;
		int idx2 = this.index(0,0);
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j < cols; j++,idx++,idx2++)
				ret.data[idx] = this.data[idx2];

			// idx2 has been incremented by cols, real_cols - cols left
			idx2+=this.real_cols - this.cols;
		}
		return ret;
	}
	/** Calculates the index of a row and column in the underlying data structure.
	 * For example for an ordinary matrix laid out in memory row by row,
	 * it will return row * this.cols() + col.
	 * @param row
	 * @param col
	 * @return The index of the row and column.
	 */
	public final int index(final int row, final int col)
	{
		return (row+row_offset)*real_cols + col + col_offset;
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
	/**
	 * Sets the entry at row, col with val.
	 * @param row
	 * @param col
	 * @param val
	 */
	public final void set(final int row, final int col, final double val)
	{
		data[index(row,col)] = val;
	}
	/**
	 * Increments the value at <code>row</code>,<code>col</code> by <code>incr</code>.
	 * @param row
	 * @param col
	 * @param incr
	 */
	public final void increment(final int row, final int col, final double incr)
	{
		data[index(row,col)] += incr;
	}
	/**
	 * Returns an integer which is the pointer to the current index.
	 * @return
	 */
	public final int ptr()
	{
		return ptr;
	}
	public final void incrementRow()
	{
		ptr += this.real_cols;
	}
	public final void incrementColumn()
	{
		++ptr;
	}
	/**
	 * Equivalent to calling ptr(); incrementRow()
	 * @return
	 */
	public final int getPtrAndIncrementRow()
	{
		int oldPtr = ptr;
		ptr += this.real_cols;
		return oldPtr;
	}
	/**
	 * Equivalent to calling incrementRow(); ptr()
	 * @return
	 */
	public final int incrementRowAndGetPtr()
	{
		ptr += this.real_cols;
		return ptr;
	}
	public final int getPtrAndIncrementColumn()
	{
		return ptr++;
	}
	public final int incrementColumnAndGetPtr()
	{
		return ++ptr;
	}
	/**
	 * Returns newly allocated matrix which is the result of summing this with other
	 * @param other
	 * @return
	 */
	public final DoubleMatrix2D plus(final DoubleMatrix2D other)
	{
		return add(this,other);
	}
	public final DoubleMatrix2D minus(final DoubleMatrix2D other)
	{
		return subtract(this,other);
	}
	public final static DoubleMatrix2D add(final DoubleMatrix2D one, final DoubleMatrix2D two)
	{
		final DoubleMatrix2D ret = new DoubleMatrix2D(one.rows(),one.cols());
		return add(one,two,ret);
	}
	public final static DoubleMatrix2D add(final DoubleMatrix2D one, final DoubleMatrix2D two, final DoubleMatrix2D ret)
	{
		if (one.cols()!=two.cols())
			throw new RuntimeException("uh oh");
		if (one.rows()!=two.rows())
			throw new RuntimeException("uh oh");
		if (one.rows()!=ret.rows()||one.cols()!=ret.cols())
			throw new RuntimeException("uh oh");
		int idx = ret.index(0,0);
		int idx1 = one.index(0,0);
		int idx2 = two.index(0,0);
		for (int i = 0; i < one.rows(); i++)
		{
			for (int j = 0; j < one.cols(); j++,idx++,idx1++,idx2++)
			{
				ret.data[idx] = one.data[idx1] + two.data[idx2];
			}
			idx1 += one.real_cols - one.cols;
			idx2 += two.real_cols - two.cols;
		}
		return ret;
	}
	public static final DoubleMatrix2D subtract(final DoubleMatrix2D minuend, final DoubleMatrix2D subtrahend)
	{
		final DoubleMatrix2D ret = new DoubleMatrix2D(minuend.rows(),minuend.cols());
		return subtract(minuend,subtrahend,ret);
	}
	public static final DoubleMatrix2D subtract(final DoubleMatrix2D minuend, final DoubleMatrix2D subtrahend, final DoubleMatrix2D ret)
	{
		if (minuend.cols()!=subtrahend.cols() || minuend.rows()!=subtrahend.rows())
			throw new RuntimeException("uh oh");
		if (minuend.rows()!=ret.rows()||minuend.cols()!=ret.cols())
			throw new RuntimeException("uh oh");
		int idx = ret.index(0,0);
		int idx1 = minuend.index(0,0);
		int idx2 = subtrahend.index(0,0);
		for (int i = 0; i < minuend.rows(); i++)
		{
			for (int j = 0; j < minuend.cols(); j++,idx++,idx1++,idx2++)
			{
				ret.data[idx] = minuend.data[idx1] - subtrahend.data[idx2];
			}
			idx1 += minuend.real_cols - minuend.cols;
			idx2 += subtrahend.real_cols - subtrahend.cols;
		}
		return ret;
	}
	public final double[] extractColumn(final int col)
	{
		final double[] ret = new double[rows];
		int idx = this.index(0,col);
		for (int i = 0; i < rows; i++,idx+=real_cols)
		{
			ret[i] = data[idx];
		}
		return ret;
	}
	public final double[] extractRow(final int row)
	{
		final double[] ret = new double[cols];
		int idx = this.index(row,0);
		for (int i = 0; i < cols; i++,idx++)
		{
			ret[i] = data[idx];
		}
		return ret;
	}
}

