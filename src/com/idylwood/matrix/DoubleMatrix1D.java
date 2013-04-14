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

// Light wrapper for double vector. Mostly for readable comparison with other linear algebra libraries
public class DoubleMatrix1D
{
	private final double[] data;
	public final int length;
	// Copy constructor
	public DoubleMatrix1D(DoubleMatrix1D other)
	{
		this(other.data,true);
	}
	public DoubleMatrix1D(final double[] data)
	{
		this(data,true);
	}
	public DoubleMatrix1D(final double[] data, final boolean copy)
	{
		if (copy)
			this.data = com.idylwood.utils.MathUtils.copyOf(data);
		else
			this.data = data;

		this.length = data.length;
	}
	// no side effects
	public double get(final int idx)
	{
		return data[idx];
	}
	// can modify array
	public void set(final int idx, final double val)
	{
		data[idx] = val;
	}
}

