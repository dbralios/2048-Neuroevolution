#ifndef _MATRIX_CONT_
#define _MATRIX_CONT_

#define _SILENCE_STDEXT_ALLOCATORS_DEPRECATION_WARNING
#include <allocators>

template <class Type, class A = std::allocator<Type>>
class matrixCont
{
public:
	typedef A allocator_type;
	typedef typename A::size_type		capacity;
	typedef typename A::size_type		size;
	typedef typename A::reference		reference;
	typedef typename A::const_reference	constReference;

	matrixCont();
	matrixCont(size r, size c);
	matrixCont(std::initializer_list<std::initializer_list<Type>> list);
	matrixCont(const matrixCont& obj);
	~matrixCont();

	matrixCont& operator=(const matrixCont& obj);
	bool operator==(const matrixCont& obj) const;
	bool operator!=(const matrixCont& obj) const;
	reference operator()(size row, size column);
	constReference operator()(size row, size column) const;

	void resize(size rowsNew, size columnsNew); 
	void resizeNoCopy(size rowsNew, size columnsNew);
	void addRow(std::initializer_list<Type> list);
	void addColumn(std::initializer_list<Type> list);

	reference at(size row, size column)
	{
		if (row >= rows)
			throw row;
		if (column >= columns)
			throw column;
		return data[getIndex(row, column)];
	}
	constReference at(size row, size column) const
	{
		if (row >= rows)
			throw row;
		if (column >= columns)
			throw column;
		return data[getIndex(row, column)];
	}

	size getRows()
	{
		return rows;
	}
	size getColumns()
	{
		return columns;
	}

	void containerSwap(matrixCont& first, matrixCont& second);

protected:
	Type *	data;
	A	allocator;

	size rows;
	size columns;

	size getIndex(size row, size column) const
	{
		return (row*columns + column);
	}
};

template <class Type, class A = std::allocator<Type>>
void matrixCont<Type, A>::containerSwap(matrixCont& first, matrixCont& second)
{
	std::swap(first.rows, second.rows);
	std::swap(first.columns, second.columns);

	std::swap(first.data, second.data);
	std::swap(first.allocator, second.allocator);
}

template <class Type, class A = std::allocator<Type>>
matrixCont<Type, A>::matrixCont() : rows(0), columns(0)
{
	data = allocator.allocate(0);
}

template <class Type, class A = std::allocator<Type>> // default constructors
matrixCont<Type, A>::matrixCont(size r, size c) : rows(r), columns(c)
{
	data = allocator.allocate(rows*columns);
}

template <class Type, class A = std::allocator<Type>>
matrixCont<Type, A>::matrixCont(std::initializer_list<std::initializer_list<Type>> list) : rows(list.size()), columns((*list.begin()).size())
{
	data = allocator.allocate(rows*columns);

	for (size i = 0; i < rows; ++i)
		for (size j = 0; j < columns; ++j)
			at(i, j) = *((*(list.begin() + i)).begin() + j);
}

template <class Type, class A = std::allocator<Type>> // copy constructor
matrixCont<Type, A>::matrixCont(const matrixCont& obj) : rows(obj.rows), columns(obj.columns)
{
	size contSize = obj.rows*obj.columns;

	data = allocator.allocate(contSize);

	for (size i = 0; i < contSize; ++i)
		data[i] = obj.data[i];
}

template <class Type, class A = std::allocator<Type>> // destructor
matrixCont<Type, A>::~matrixCont()
{
	allocator.destroy(data);
	allocator.deallocate(data, rows*columns);
}

template <class Type, class A = std::allocator<Type>>
matrixCont<Type, A>& matrixCont<Type, A>::operator=(const matrixCont& obj)
{
	allocator.destroy(data);
	allocator.deallocate(data, rows*columns);

	size contSize = obj.rows*obj.columns;
	data = allocator.allocate(contSize);
	rows = obj.rows;
	columns = obj.columns;

	for (size i = 0; i < contSize; ++i)
		data[i] = obj.data[i];
	return *this;
}

template <class Type, class A = std::allocator<Type>>
bool matrixCont<Type, A>::operator==(const matrixCont& obj) const
{
	if (rows != obj.rows || columns != obj.columns)
		return false;

	size contSize = rows*columns;

	for (size i = 0; i < contSize; ++i)
		if (data[i] != obj.data[i])
			return false;

	return true;
}

template <class Type, class A = std::allocator<Type>>
bool matrixCont<Type, A>::operator!=(const matrixCont& obj) const
{
	return (!(*this == obj));
}

template <class Type, class A = std::allocator<Type>>
typename matrixCont<Type, A>::reference matrixCont<Type, A>::operator()(size row, size column)
{
	if (row >= rows)
		throw row;
	if (column >= columns)
		throw column;
	return data[getIndex(row, column)];
}

template <class Type, class A = std::allocator<Type>>
typename matrixCont<Type, A>::constReference matrixCont<Type, A>::operator()(size row, size column) const
{
	if (row >= rows)
		throw row;
	if (column >= columns)
		throw column;
	return data[getIndex(row, column)];
}

template <class Type, class A = std::allocator<Type>>
void matrixCont<Type, A>::resize(size rowsNew, size columnsNew)
{
	size newSize = rowsNew*columnsNew;
	size oldSize = rows*columns;

	if (newSize < 0)
		throw newSize;

	Type * temp;
	A tempAlloc;

	temp = tempAlloc.allocate(oldSize);
	for (size i = 0; i < oldSize; ++i)
		temp[i] = data[i];

	allocator.destroy(data);
	allocator.deallocate(data, oldSize);

	data = allocator.allocate(newSize);

	size indNew;
	for (size indOld = 0; indOld < oldSize; ++indOld)
	{
		indNew = (indOld / columns)*columnsNew + indOld % columns;
		data[indNew] = temp[indOld];
	}

	tempAlloc.destroy(temp);
	tempAlloc.deallocate(temp, oldSize);

	columns = columnsNew;
	rows = rowsNew;
}

template <class Type, class A = std::allocator<Type>>
void matrixCont<Type, A>::resizeNoCopy(size rowsNew, size columnsNew)
{
	size newSize = rowsNew*columnsNew;
	size oldSize = rows*columns;

	allocator.destroy(data);
	allocator.deallocate(data, oldSize);
	data = allocator.allocate(newSize);

	columns = columnsNew;
	rows = rowsNew;
}

template <class Type, class A = std::allocator<Type>>
void matrixCont<Type, A>::addRow(std::initializer_list<Type> list)
{
	if (columns > 0)
		resize(rows + 1, columns);
	else
		resize(rows + 1, list.size());

	for (size i = 0; i < columns; ++i)
		at(rows - 1, i) = *(list.begin() + i);
}

template <class Type, class A = std::allocator<Type>>
void matrixCont<Type, A>::addColumn(std::initializer_list<Type> list)
{
	if (rows > 0)
		resize(rows, columns + 1);
	else
		resize(list.size(), columns + 1);

	for (size i = 0; i < rows; ++i)
		at(i, columns - 1) = *(list.begin() + i);
}
#endif