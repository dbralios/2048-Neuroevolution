#ifndef _MATRIX_
#define _MATRIX_

#include "matrixContainer.h"

template <class Type, class A = std::allocator<Type>>
class matrix : public matrixCont<Type, A>
{
public:

	using size		= typename matrixCont<Type, A>::size;
	using reference = typename matrixCont<Type, A>::reference;
	using matrixCont<Type, A>::getColumns;
	using matrixCont<Type, A>::getRows;
	using matrixCont<Type, A>::at;
	using matrixCont<Type, A>::columns;
	using matrixCont<Type, A>::rows;
	using matrixCont<Type, A>::data;

	matrix() {}
	matrix(size r, size c) : matrixCont<Type, A>::matrixCont(r, c) {}
	matrix(std::initializer_list<std::initializer_list<Type>> list) : matrixCont<Type, A>::matrixCont(list) {};
	matrix(const matrix& obj) : matrixCont<Type, A>::matrixCont(obj) {};
	~matrix() 
	{

	};

	matrix operator+(const matrix& obj);
	matrix operator-(const matrix& obj);
	matrix operator*(const matrix& obj);
	matrix operator*(const Type scalar);
	matrix operator+=(const matrix& obj);
	matrix operator-=(const matrix& obj);
	matrix operator*=(const matrix& obj);
	matrix operator*=(const Type scalar);

	matrix trans();
	matrix inv();
	Type trace();
	Type det();

	void swapRows(size first, size second);
	void swapCols(size first, size second);
	void multiplyRowBy(size row, Type multiplier);
	void multiplyColBy(size col, Type multiplier);
	void addRows(size first, size second, Type scalar);
	void addCols(size first, size second, Type scalar);

	matrix getPrincipalSubmatrix(size row, size column);

	class rowIterator;
	class colIterator;
	
	class iterator
	{
	private:
		matrix* parent;
		size index;
	
		void iteratorSwap(iterator& first, iterator& second);
	public:
		iterator(size ind, matrix* p) : parent(p), index(ind) {}
		iterator(const iterator& it) : parent(it.parent), index(it.index) {};
		~iterator() {};

		iterator& operator=(iterator);
		bool operator==(const iterator&) const;
		bool operator!=(const iterator&) const;
		bool operator<(const iterator&) const;
		bool operator>(const iterator&) const;
		bool operator<=(const iterator&) const;
		bool operator>=(const iterator&) const;

		iterator& operator++();
		iterator operator++(int);
		iterator& operator--();
		iterator operator--(int);
		reference operator*() const;
	};

	iterator begin();
	iterator end();
private:
};

template <class Type, class A = std::allocator<Type>>
matrix<Type, A> matrix<Type, A>::operator+=(const matrix& obj)
{
	size contSize = rows * columns;

	for (size i = 0; i < contSize; i++)
		data[i] += obj.data[i];

	return *this;
}

template <class Type, class A = std::allocator<Type>>
matrix<Type, A> matrix<Type, A>::operator+(const matrix& obj)
{
	matrix sum(rows, columns);

	size contSize = rows * columns;

	for (size i = 0; i < contSize; i++)
		sum.data[i] = (data[i] + obj.data[i]);

	return sum;
}

template <class Type, class A = std::allocator<Type>>
matrix<Type, A> matrix<Type, A>::operator-=(const matrix& obj)
{
	size contSize = rows * columns;

	for (size i = 0; i < contSize; i++)
		data[i] -= obj.data[i];

	return *this;
}

template <class Type, class A = std::allocator<Type>>
matrix<Type, A> matrix<Type, A>::operator-(const matrix& obj)
{
	matrix difference(rows, columns);

	size contSize = rows * columns;

	for (size i = 0; i < contSize; i++)
		difference.data[i] = (data[i] - obj.data[i]);

	return difference;
}

template <class Type, class A = std::allocator<Type>>
matrix<Type, A> matrix<Type, A>::operator*=(const matrix& obj)
{
	*this = *this * obj;
	return *this;
}

template <class Type, class A = std::allocator<Type>>
matrix<Type, A> matrix<Type, A>::operator*(const matrix& obj)
{
	matrix product(rows, obj.columns);

	size contSize = rows * columns;

	for (size rowInd = 0; rowInd < product.getRows(); rowInd++)
	{
		for (size colInd = 0; colInd < product.getColumns(); colInd++)
		{
			product(rowInd, colInd) = 0;
			for (size sumInd = 0; sumInd < getColumns(); sumInd++)
				product(rowInd, colInd) += at(rowInd, sumInd)*obj(sumInd, colInd);
		}
	}

	return product;
}

template <class Type, class A = std::allocator<Type>>
matrix<Type, A> matrix<Type, A>::operator*(const Type scalar)
{
	matrix product(rows, columns);

	size contSize = rows * columns;

	for (size i = 0; i < contSize; i++)
		product.data[i] = data[i] * scalar;

	return product;
}

template <class Type, class A = std::allocator<Type>>
matrix<Type, A> matrix<Type, A>::operator*=(const Type scalar)
{
	size contSize = rows * columns;

	for (size i = 0; i < contSize; i++)
		data[i] *= scalar;

	return *this;
}

template <class Type, class A = std::allocator<Type>>
matrix<Type, A> matrix<Type, A>::trans()
{
	matrix trans(columns, rows);

	for (size rowInd = 0; rowInd < trans.getRows(); rowInd++)
		for (size colInd = 0; colInd < trans.getColumns(); colInd++)
			trans(rowInd, colInd) = at(colInd, rowInd);
		
	return trans;
}

template <class Type, class A = std::allocator<Type>>
Type matrix<Type, A>::trace()
{
	Type trace = 0;

	for (size rowInd = 0; rowInd < this.getRows(); rowInd++)
		trace += this.at(rowInd, rowInd);

	return trace;
}

template <class Type, class A = std::allocator<Type>>
Type matrix<Type, A>::det()
{
	size N = rows;

	if (N == 1) // a
		return at(0, 0);

	if (N == 2) // ad - bc
		return (at(0, 0)*at(1, 1) - at(0, 1)*at(1, 0));
	
	Type det = 0;
	for (size i = 0; i < N; i++)
	{
		matrix minor(N - 1, N - 1);
		minor = getPrincipalSubmatrix(0,i); // Get submatrix by deleting row 0 and column i
		(i % 2 == 0) ? det += at(0, i)*minor.det() : det += (-1)*at(0, i)*minor.det(); // Find det recursively
	}
	return det;
}

template <class Type, class A = std::allocator<Type>>
matrix<Type, A> matrix<Type, A>::getPrincipalSubmatrix(size row, size column)
{
	matrix subMatrix(rows - 1, columns - 1);
	size subMatrixSize = (rows - 1) * (columns - 1);

	size oldInd = -1;
	size oldCurRow, oldCurCol;
	for (size i = 0; i < subMatrixSize; i++)
	{
		oldInd++;

		oldCurRow = (oldInd - (oldInd % columns)) / columns;
		if (oldCurRow == row) 
			oldInd += columns; // Skip row

		oldCurCol = oldInd % columns;
		if (oldCurCol == column) 
			oldInd++; // Skip column element

		subMatrix.data[i] = data[oldInd];
	}
	return subMatrix;
}

template <class Type, class A = std::allocator<Type>>
void matrix<Type, A>::swapRows(size first, size second)
{
	Type * tempRow;
	A tempAlloc;
	size rowSize = columns;

	tempRow = tempAlloc.allocate(rowSize);

	for (size i = 0; i < rowSize; i++)
	{
		tempRow[i] = at(first, i);
		at(first, i) = at(second, i);
		at(second, i) = tempRow[i];
	}

	tempAlloc.destroy(tempRow);
	tempAlloc.deallocate(tempRow, rowSize);
}

template <class Type, class A = std::allocator<Type>>
void matrix<Type, A>::swapCols(size first, size second)
{
	Type * tempCol;
	A tempAlloc;
	size colSize = rows;

	tempCol = tempAlloc.allocate(colSize);

	for (size i = 0; i < colSize; i++)
	{
		tempCol[i] = at(i, first);
		at(i, first) = at(i, second);
		at(i, second) = tempCol[i];
	}

	tempAlloc.destroy(tempCol);
	tempAlloc.deallocate(tempCol, colSize);
}

template <class Type, class A = std::allocator<Type>>
void matrix<Type, A>::multiplyRowBy(size row, Type multiplier)
{
	size rowSize = columns;

	for (size i = 0; i < rowSize; i++)
		at(row, i) *= multiplier;
}

template <class Type, class A = std::allocator<Type>>
void matrix<Type, A>::multiplyColBy(size col, Type multiplier)
{
	size colSize = rows;

	for (size i = 0; i < colSize; i++)
		at(i, col) *= multiplier;
}

template <class Type, class A = std::allocator<Type>>
void matrix<Type, A>::addRows(size first, size second, Type scalar)
{
	size rowSize = columns;

	for (size i = 0; i < rowSize; i++)
		at(first, i) += scalar * at(second, i);
}

template <class Type, class A = std::allocator<Type>>
void matrix<Type, A>::addCols(size first, size second, Type scalar)
{
	size colSize = rows;

	for (size i = 0; i < colSize; i++)
		at(i, first) += scalar * at(i, second);
}

template <class Type, class A = std::allocator<Type>>
typename matrix<Type, A>::iterator matrix<Type, A>::begin()
{
	iterator begin(0, this);
	return begin;
}

template <class Type, class A = std::allocator<Type>>
typename matrix<Type, A>::iterator matrix<Type, A>::end()
{
	size contSize = this->rows * this->columns;
	iterator end(contSize, this);
	return end;
}

template <class Type, class A = std::allocator<Type>>
void matrix<Type, A>::iterator::iteratorSwap(iterator& first, iterator& second)
{
	std::swap(first.index, second.index);
	std::swap(first.parent, second.parent);
}

template <class Type, class A = std::allocator<Type>>
typename matrix<Type, A>::iterator& matrix<Type, A>::iterator::operator=(iterator it)
{
	iteratorSwap(*this, it);
	return *this;
}

template <class Type, class A = std::allocator<Type>>
bool matrix<Type, A>::iterator::operator==(const iterator &it) const
{
	if (parent != it.parent)
		return false;
	if (index != it.index)
		return false;
	return true;
}

template <class Type, class A = std::allocator<Type>>
bool matrix<Type, A>::iterator::operator!=(const iterator &it) const
{
	return (!(*this == it));
}

template <class Type, class A = std::allocator<Type>>
bool matrix<Type, A>::iterator::operator<(const iterator &it) const
{
	return (index < it.index);
}

template <class Type, class A = std::allocator<Type>>
bool matrix<Type, A>::iterator::operator>(const iterator &it) const
{
	return (index > it.index);
}

template <class Type, class A = std::allocator<Type>>
bool matrix<Type, A>::iterator::operator<=(const iterator &it) const
{
	return (index <= it.index);
}

template <class Type, class A = std::allocator<Type>>
bool matrix<Type, A>::iterator::operator>=(const iterator &it) const
{
	return (index >= it.index);
}

template <class Type, class A = std::allocator<Type>>
typename matrix<Type, A>::iterator& matrix<Type, A>::iterator::operator++()
{
	index++;
	return *this;
}

template <class Type, class A = std::allocator<Type>>
typename matrix<Type, A>::iterator matrix<Type, A>::iterator::operator++(int)
{
	iterator temp = *this;
	++*this;
	return temp;
}

template <class Type, class A = std::allocator<Type>>
typename matrix<Type, A>::iterator& matrix<Type, A>::iterator::operator--()
{
	index--;
	return *this;
}

template <class Type, class A = std::allocator<Type>>
typename matrix<Type, A>::iterator matrix<Type, A>::iterator::operator--(int)
{
	iterator temp = *this;
	--*this;
	return temp;
}

template <class Type, class A = std::allocator<Type>>
typename matrix<Type, A>::reference matrix<Type, A>::iterator::operator*() const
{
	size row, col;
	col = index % parent->getColumns();
	row = (index - col) / parent->getColumns();

	if (row >= parent->rows || col >= parent->columns)
		throw index;
	return parent->at(row, col);
}
#endif