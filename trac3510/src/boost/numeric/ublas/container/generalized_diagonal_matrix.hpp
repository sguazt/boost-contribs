/**
 *  \file generalized_diagonal_matrix.hpp
 *
 *  \brief Generalized diagonal matrix.
 *
 *  Copyright (c) 2009, Marco Guazzone
 *
 *  Distributed under the Boost Software License, Version 1.0. (See
 *  accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 *
 *  \author Marco Guazzone, marco.guazzone@gmail.com
 */

#ifndef BOOST_NUMERIC_UBLAS_CONTAINER_DIAGONAL_MATRIX_HPP
#define BOOST_NUMERIC_UBLAS_CONTAINER_DIAGONAL_MATRIX_HPP

#include <algorithm>
//#include <boost/numeric/ublas/detail/matrix_assign.hpp>
#include <boost/numeric/ublas/detail/temporary.hpp>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/storage.hpp>


namespace boost { namespace numeric { namespace ublas {

/**
 * \brief Generalized diagonal matrix.
 * \tparam ValueT The type of matrix values.
 * \tparam LayoutT The matrix layout type.
 *  Default to \c row_major.
 * \tparam ArrayT The type of the array storing the matrix.
 *  Default to \c unbounded_array<ValueT>
 *
 * This class can be used to represent:
 * - <em>diagonal matrices</em>: a square matrix \f$A\f$ of order \f$n\f$ is a
 *   <em>diagonal matrix</em> if \f$a_{ij}=0\f$ for \f$i \ne j\f$,
 *   \f$i,j=0,\ldots,n\f$;
 * - <em>rectangular diagonal matrices</em>: an \f$m \times n\f$ matrix \f$A\f$
 *   is a <em>rectangular diagonal matrix</em> if \f$a_{ij}=0\f$ for
 *   \f$i \ne j\f$,\f$i=1,\ldots,m\f$ and \f$j=1,\ldots,n\f$. Thus the only
 *   entries in \f$A\f$ that may be non-zero are the \f$d_{ii}\f$,
 *   \f$i=1,\ldots,\min\{m,n\}\f$;
 * - <em>sub-diagonal matrices</em>: a square matrix \f$A\f$ of order
 *   \f$n\f$ is a <em>k-th sub-diagonal matrix</em> if \f$a_{ij}=0\f$ for
 *   \f$i \ne (j+k)\f$, \f$i=1,\ldots,n\f$ and \f$j=1,\ldots,n\f$;
 * - <em>super-diagonal matrices</em>: a square matrix \f$A\f$ of order \f$n\f$
 *   is a <em>k-th super-diagonal matrix</em> if \f$a_{ij}=0\f$ for
 *   \f$(i+k) \ne j\f$, \f$i=1,\ldots,n\f$ and \f$j=1,\ldots,n\f$;
 * - <em>rectangular sub-diagonal matrices</em>: a \f$m \times n\f$ matrix
 *   \f$A\f$ is a <em>k-th rectangular sub-diagonal matrix</em> if
 *   \f$a_{ij}=0\f$ for \f$i \ne (j+k)\f$, \f$i=1,\ldots,m\f$ and
 *   \f$j=1,\ldots,n\f$;
 * - <em>rectangular super-diagonal matrices</em>: a \f$m \times n\f$ matrix
 *   \f$A\f$ is a <em>k-th rectangular super-diagonal matrix</em> if
 *   \f$a_{ij}=0\f$ for \f$(i+k) \ne j\f$, \f$i=1,\ldots,m\f$ and
 *   \f$j=1,\ldots,n\f$;
 * .
 *
 * \author Marco Guazzone, marco.guazzone@gmail.com
 */
template <typename ValueT, typename LayoutT = row_major, typename ArrayT = unbounded_array<ValueT> >
class generalized_diagonal_matrix: public matrix_container<generalized_diagonal_matrix<ValueT, LayoutT, ArrayT> >
{

	private: typedef ValueT *pointer;
	private: typedef LayoutT layout_type;
	private: typedef generalized_diagonal_matrix<ValueT, LayoutT, ArrayT> self_type;
	public: typedef typename ArrayT::size_type size_type;
	public: typedef typename ArrayT::difference_type difference_type;
	public: typedef ValueT value_type;
	public: typedef const ValueT &const_reference;
	public: typedef ValueT &reference;
	public: typedef ArrayT array_type;
	public: typedef const matrix_reference<const self_type> const_closure_type;
	public: typedef matrix_reference<self_type> closure_type;
	public: typedef vector<ValueT, ArrayT> vector_temporary_type;
	public: typedef matrix<ValueT, LayoutT, ArrayT> matrix_temporary_type;  // general sub-matrix
	public: typedef packed_tag storage_category;
	public: typedef typename LayoutT::orientation_category orientation_category;
	private: typedef const value_type const_value_type;
	// Iterator types
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
	public: typedef indexed_iterator1<self_type, packed_random_access_iterator_tag> iterator1;
	public: typedef indexed_iterator2<self_type, packed_random_access_iterator_tag> iterator2;
	public: typedef indexed_const_iterator1<self_type, packed_random_access_iterator_tag> const_iterator1;
	public: typedef indexed_const_iterator2<self_type, packed_random_access_iterator_tag> const_iterator2;
#else
	public: class const_iterator1;
	public: class iterator1;
	public: class const_iterator2;
	public: class iterator2;
#endif
	public: typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
	public: typedef reverse_iterator_base1<iterator1> reverse_iterator1;
	public: typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
	public: typedef reverse_iterator_base2<iterator2> reverse_iterator2;


#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
	public: using matrix_container<self_type>::operator();
#endif


	//@{ Construction and destruction

	public: BOOST_UBLAS_INLINE
		generalized_diagonal_matrix ()
		: matrix_container<self_type>(),
		  size1_(0),
		  size2_(0),
		  k_(0),
		  r_(0),
		  c_(0),
		  data_(0)
	{
		// Empty
	}


	/**
	 * \brief Create a diagonal matrix of order \a size whose non-zero elements
	 *  are on diagonal \a k.
	 */
	public: BOOST_UBLAS_INLINE
		generalized_diagonal_matrix(size_type size, difference_type k=0)
		: matrix_container<self_type>(),
		  size1_(size),
		  size2_(size),
		  k_(k),
		  r_(k < 0 ? -k : 0),
		  c_(k > 0 ?  k : 0),
		  data_(size - (k >= 0 ? k : -k))
	{
		// preconditions
		BOOST_UBLAS_CHECK(size_type(k >= 0 ? k : -k) < size, bad_size());
	}


	/**
	 * \brief Create a rectangular diagonal matrix of size \a size1 by \a size2
	 * whose non-zero elements are on diagonal \a k.
	 */
	public: BOOST_UBLAS_INLINE
		generalized_diagonal_matrix(size_type size1, size_type size2, difference_type k)
		: matrix_container<self_type>(),
		  size1_(size1),
		  size2_(size2),
		  k_(k),
		  r_(k < 0 ? -k : 0),
		  c_(k > 0 ?  k : 0),
		  data_(std::min(size1 - r_, size2 - c_))
	{
		// preconditions
		BOOST_UBLAS_CHECK(r_ < size1_, bad_size());
		BOOST_UBLAS_CHECK(c_ < size2_, bad_size());
	}


	public: BOOST_UBLAS_INLINE
		generalized_diagonal_matrix(size_type size, difference_type k, array_type const& data)
		: matrix_container<self_type>(),
		  size1_(size),
		  size2_(size),
		  k_(k),
		  r_(k < 0 ? -k : 0),
		  c_(k > 0 ?  k : 0),
		  data_(data)
	{
		// preconditions
		BOOST_UBLAS_CHECK(r_ < size1_, bad_size());
		BOOST_UBLAS_CHECK(c_ < size2_, bad_size());

		size_type real_size = size - (r_ + c_);

		if (data_.size() > real_size)
		{
			data_.resize(real_size, value_type(0));
		}
	}


	public: BOOST_UBLAS_INLINE
		generalized_diagonal_matrix(size_type size1, size_type size2, difference_type k, array_type const& data)
		: matrix_container<self_type>(),
		  size1_(size1),
		  size2_(size2),
		  k_(k),
		  r_(k < 0 ? -k : 0),
		  c_(k > 0 ?  k : 0),
		  data_(data)
	{
		// preconditions
		BOOST_UBLAS_CHECK(r_ < size1_, bad_size());
		BOOST_UBLAS_CHECK(c_ < size2_, bad_size());

		size_type real_size = std::min(size1 - r_, size2 - c_);

		if (data_.size() > real_size)
		{
			data_.resize(real_size, value_type(0));
		}
	}


	public: BOOST_UBLAS_INLINE
		generalized_diagonal_matrix(generalized_diagonal_matrix const& m)
		: matrix_container<self_type>(),
		  size1_(m.size1_),
		  size2_(m.size2_),
		  k_(m.k_),
		  r_(m.r_),
		  c_(m.c_),
		  data_(m.data_)
	{
		// Empty
		// preconditions
		BOOST_UBLAS_CHECK(r_ < size1_, bad_size());
		BOOST_UBLAS_CHECK(c_ < size2_, bad_size());
	}


	public: template <typename ExprT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix(matrix_expression<ExprT> const& me, difference_type k=0)
		: matrix_container<self_type>(),
		  //size_(std::min(me().size1(), me().size2())),
		  size1_(me().size1()),
		  size2_(me().size2()),
		  k_(k),
		  r_(k < 0 ? -k : 0),
		  c_(k > 0 ?  k : 0),
		  //data_(size_ - (k >= 0 ? k : -k))
		  data_(std::min(me().size1() - r_, me().size2() - c_))
	{
		// preconditions
		BOOST_UBLAS_CHECK(r_ < size1_, bad_size());
		BOOST_UBLAS_CHECK(c_ < size2_, bad_size());

		matrix_assign<scalar_assign>(*this, me);
	}


	public: template <typename ExprT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix(vector_expression<ExprT> const& ve, difference_type k=0)
		: matrix_container<self_type>(),
		  size1_(ve().size() + (k < 0 ? -k : 0)),
		  size2_(ve().size() + (k >= 0 ? k : 0)),
		  k_(k),
		  r_(k < 0 ? -k : 0),
		  c_(k > 0 ?  k : 0),
		  data_(ve().size())
	{
		// preconditions
		BOOST_UBLAS_CHECK(r_ < size1_, bad_size());
		BOOST_UBLAS_CHECK(c_ < size2_, bad_size());

		typedef typename ExprT::size_type ve_size_type;

		ve_size_type ve_size = ve().size();

		for (
			ve_size_type i = 0;
			i < ve_size;
			++i
		) {
			(*this)(i+r_, i+c_) = ve()(i);
		}
	}


	//@} Construction/Destruction

	//@{ Accessors


	public: BOOST_UBLAS_INLINE
		size_type size1() const
	{
		return size1_;
	}


	public: BOOST_UBLAS_INLINE
		size_type size2() const
	{
		return size2_;
	}


	public: BOOST_UBLAS_INLINE
		difference_type offset() const
	{
		return k_;
	}


	// Storage accessors
	public: BOOST_UBLAS_INLINE
		array_type const& data() const
	{
		return data_;
	}


	public: BOOST_UBLAS_INLINE
		array_type& data()
	{
		return data_;
	}


	//@} Accessors

	//@{ Resizing


	public: BOOST_UBLAS_INLINE
		void resize(size_type size, difference_type k=0, bool preserve=true)
	{
		if (preserve)
		{
			self_type temporary(size, k);
			detail::matrix_resize_preserve<layout_type>(*this, temporary);
		}
		else
		{
			size1_ = size;
			size2_ = size;
			k_ = k;
			r_ = k < 0 ? -k : 0;
			c_ = k > 0 ?  k : 0;
			data().resize(size - (k > 0 ? k : -k));
		}
	}


	public: BOOST_UBLAS_INLINE
		void resize(size_type size1, size_type size2, difference_type k=0, bool preserve=true)
	{
		if (preserve)
		{
			self_type temporary(size1, size2, k);
			detail::matrix_resize_preserve<layout_type>(*this, temporary);
		}
		else
		{
			size1_ = size1;
			size2_ = size2;
			k_ = k;
			r_ = k < 0 ? -k : 0;
			c_ = k > 0 ?  k : 0;
			data().resize(std::min(size1 - r_, size2 - c_));
		}
	}


	public: BOOST_UBLAS_INLINE
		void resize_packed_preserve(size_type size, difference_type k=0)
	{
		size1_ = size;
		size2_ = size;
		k_ = k;
		r_ = k < 0 ? -k : 0;
		c_ = k > 0 ?  k : 0;
		data().resize(size - (k > 0 ? k : -k), value_type());
	}


	public: BOOST_UBLAS_INLINE
		void resize_packed_preserve(size_type size1, size_type size2, difference_type k=0)
	{
		size1_ = size1;
		size2_ = size2;
		k_ = k;
		r_ = k < 0 ? -k : 0;
		c_ = k > 0 ?  k : 0;
		data().resize(std::min(size1 - r_, size2 - c_), value_type());
	}


	//@} Resizing

	//@{ Element access


	public: BOOST_UBLAS_INLINE
		const_reference operator()(size_type i, size_type j) const
	{
		BOOST_UBLAS_CHECK(i < size1_, bad_index());
		BOOST_UBLAS_CHECK(j < size2_, bad_index());

		const difference_type dr(i-r_);
		const difference_type dc(j-c_);

		if (dr == dc)
		{
			if (k_ > 0)
			{
				return data()[layout_type::element(i, size1_, 0, 1)];
			}
			else
			{
				return data()[layout_type::element(0, 1, j, size2_)];
			}
		}

		return zero_;
	}


	public: BOOST_UBLAS_INLINE
		reference at_element(size_type i, size_type j)
	{
		BOOST_UBLAS_CHECK(i < size1_, bad_index());
		BOOST_UBLAS_CHECK(j < size2_, bad_index());

		if (k_ > 0)
		{
			return data()[layout_type::element(i, size1_, 0, 1)];
		}
		else
		{
			return data()[layout_type::element(0, 1, j, size2_)];
		}
	}


	public: BOOST_UBLAS_INLINE
		reference operator()(size_type i, size_type j)
	{
		BOOST_UBLAS_CHECK(i < size1_, bad_index());
		BOOST_UBLAS_CHECK(j < size2_, bad_index());

		const difference_type dr(i-r_);
		const difference_type dc(j-c_);

		if (dr != dc)
		{
			bad_index().raise();
		}

		if (k_ > 0)
		{
			return data()[layout_type::element(i, size1_, 0, 1)];
		}
		else
		{
			return data()[layout_type::element(0, 1, j, size2_)];
		}
	}


	//@} Element access

	//@{ Element assignment


	public: BOOST_UBLAS_INLINE
		reference insert_element(size_type i, size_type j, const_reference t)
	{
		return (operator()(i, j) = t);
	}


	public: BOOST_UBLAS_INLINE
		void erase_element (size_type i, size_type j)
	{
		operator()(i, j) = value_type/*zero*/();
	}


	//@} Element assignment

	//@{ Zeroing


	public: BOOST_UBLAS_INLINE
		void clear()
	{
		std::fill(data().begin(), data().end(), value_type/*zero*/());
	}


	//@} Zeroing

	//@{ Assignment


	public: BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& operator=(generalized_diagonal_matrix const& m)
	{
		size1_ = m.size1_;
		size2_ = m.size2_;
		k_ = m.k_;
		r_ = m.r_;
		c_ = m.c_;
		data() = m.data();

		return *this;
	}


	public: BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& assign_temporary(generalized_diagonal_matrix& m)
	{
		swap(m);

		return *this;
	}


	public: template <typename ExprT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& operator=(matrix_expression<ExprT> const& me)
	{
		self_type temporary(me, k_);

		return assign_temporary(temporary);
	}


	public: template <typename ExprT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& assign(matrix_expression<ExprT> const& me)
	{
		matrix_assign<scalar_assign>(*this, me);

		return *this;
	}


	public: template <typename ExprT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& operator=(vector_expression<ExprT> const& ve)
	{
		self_type temporary(ve, k_);

		return assign_temporary(temporary);
	}


	public: template <typename ExprT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& assign(vector_expression<ExprT> const& ve)
	{
		return operator=(ve);
	}


	public: template <typename ExprT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& operator+=(matrix_expression<ExprT> const& me)
	{
		self_type temporary(*this + me, k_);

		return assign_temporary(temporary);
	}


	public: template <typename ExprT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& plus_assign(matrix_expression<ExprT> const& me)
	{
		matrix_assign<scalar_plus_assign>(*this, me);

		return *this;
	}


	public: template <typename ExprT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& operator-=(matrix_expression<ExprT> const& me)
	{
		self_type temporary(*this - me, k_);

		return assign_temporary(temporary);
	}


	public: template <typename ExprT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& minus_assign(matrix_expression<ExprT> const& me)
	{
		matrix_assign<scalar_minus_assign>(*this, me);

		return *this;
	}


	public: template <typename ScalarT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& operator*=(ScalarT const& se)
	{
		matrix_assign_scalar<scalar_multiplies_assign>(*this, se);

		return *this;
	}


	public: template <typename ScalarT>
		BOOST_UBLAS_INLINE
		generalized_diagonal_matrix& operator/=(ScalarT const& se)
	{
		matrix_assign_scalar<scalar_divides_assign>(*this, se);

		return *this;
	}


	//@} Assignment

	//@{ Swapping


	public: BOOST_UBLAS_INLINE
		void swap(generalized_diagonal_matrix& m)
	{
		if (this != &m)
		{
			std::swap(size1_, m.size1_);
			std::swap(size2_, m.size2_);
			std::swap(k_, m.k_);
			std::swap(r_, m.r_);
			std::swap(c_, m.c_);
			data().swap(m.data());
		}
	}


	public: BOOST_UBLAS_INLINE
		friend void swap(generalized_diagonal_matrix& m1, generalized_diagonal_matrix& m2)
	{
		m1.swap(m2);
	}


	//@} Swapping

	//@{ Element lookup

	public: BOOST_UBLAS_INLINE
		const_iterator1 find1(int /*rank*/, size_type i, size_type j) const
	{
		if (i < r_)
		{
			i = r_;
		}
		if (i > (r_ + std::min(size1_-r_, size2_-c_)))
		{
			i = r_ + std::min(size1_-r_, size2_-c_);
		}
//		if (j < c_)
//		{
//			j = c_;
//		}
//		if (j > (c_ + std::min(size1_-r_, size2_-c_)))
//		{
//			j = c_ + std::min(size1_-r_, size2_-c_);
//		}

		return const_iterator1(*this, i, j);
	}


	public: BOOST_UBLAS_INLINE
		iterator1 find1(int /*rank*/, size_type i, size_type j)
	{
		if (i < r_)
		{
			i = r_;
		}
		if (i > (r_ + std::min(size1_-r_, size2_-c_)))
		{
			i = r_ + std::min(size1_-r_, size2_-c_);
		}
//		if (j < c_)
//		{
//			j = c_;
//		}
//		if (j > (c_ + std::min(size1_-r_, size2_-c_)))
//		{
//			j = c_ + std::min(size1_-r_, size2_-c_);
//		}

		return iterator1(*this, i, j);
	}


	public: BOOST_UBLAS_INLINE
		const_iterator2 find2(int /*rank*/, size_type i, size_type j) const
	{
//		if (i < r_)
//		{
//			i = r_;
//		}
//		if (i > (r_ + std::min(size1_-r_, size2_-c_)))
//		{
//			i = r_ + std::min(size1_-r_, size2_-c_);
//		}
		if (j < c_)
		{
			j = c_;
		}
		if (j > (c_ + std::min(size1_-r_, size2_-c_)))
		{
			j = c_ + std::min(size1_-r_, size2_-c_);
		}

		return const_iterator2(*this, i, j);
	}


	public: BOOST_UBLAS_INLINE
		iterator2 find2(int /*rank*/, size_type i, size_type j)
	{
//		if (i < r_)
//		{
//			i = r_;
//		}
//		if (i > (r_ + std::min(size1_-r_, size2_-c_)))
//		{
//			i = r_ + std::min(size1_-r_, size2_-c_);
//		}
		if (j < c_)
		{
			j = c_;
		}
		if (j > (c_ + std::min(size1_-r_, size2_-c_)))
		{
			j = c_ + std::min(size1_-r_, size2_-c_);
		}

		return iterator2(*this, i, j);
	}


	//@} Element lookup

	//@{ Forward Iterators

	public: BOOST_UBLAS_INLINE
		const_iterator1 begin1() const
	{
		return find1(0, r_, c_);
	}


	public: BOOST_UBLAS_INLINE
		const_iterator1 end1() const
	{
		return find1(0, r_ + std::min(size1_-r_, size2_-c_), c_);
	}


	public: BOOST_UBLAS_INLINE
		iterator1 begin1()
	{
		return find1(0, r_, c_);
	}


	public: BOOST_UBLAS_INLINE
		iterator1 end1()
	{
		return find1(0, r_ + std::min(size1_-r_, size2_-c_), c_);
	}


	public: BOOST_UBLAS_INLINE
		const_iterator2 begin2() const
	{
		return find2(0, r_, c_);
	}


	public: BOOST_UBLAS_INLINE
		const_iterator2 end2() const
	{
		return find2(0, r_, c_ + std::min(size1_-r_, size2_-c_));
	}


	public: BOOST_UBLAS_INLINE
		iterator2 begin2()
	{
		return find2(0, r_, c_);
	}


	public: BOOST_UBLAS_INLINE
		iterator2 end2()
	{
		return find2(0, r_, c_ + std::min(size1_-r_, size2_-c_));
	}


	//@} Forward Iterators

	//@{ Reverse iterators

	public: BOOST_UBLAS_INLINE
		const_reverse_iterator1 rbegin1() const
	{
		return const_reverse_iterator1(end1());
	}


	public: BOOST_UBLAS_INLINE
		const_reverse_iterator1 rend1() const
	{
		return const_reverse_iterator1(begin1());
	}


	public: BOOST_UBLAS_INLINE
		reverse_iterator1 rbegin1()
	{
		return reverse_iterator1(end1());
	}


	public: BOOST_UBLAS_INLINE
		reverse_iterator1 rend1()
	{
		return reverse_iterator1(begin1());
	}


	public: BOOST_UBLAS_INLINE
		const_reverse_iterator2 rbegin2() const
	{
		return const_reverse_iterator2(end2());
	}


	public: BOOST_UBLAS_INLINE
		const_reverse_iterator2 rend2() const
	{
		return const_reverse_iterator2(begin2());
	}


	public: BOOST_UBLAS_INLINE
		reverse_iterator2 rbegin2()
	{
		return reverse_iterator2(end2());
	}


	public: BOOST_UBLAS_INLINE
		reverse_iterator2 rend2()
	{
		return reverse_iterator2(begin2());
	}


	//@} Reverse iterators

	//@{ Iterators types

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
	public: class const_iterator1: 	public container_const_reference<generalized_diagonal_matrix>,
									public random_access_iterator_base<packed_random_access_iterator_tag, const_iterator1, value_type>
	{
		public: typedef typename generalized_diagonal_matrix::value_type value_type;
		public: typedef typename generalized_diagonal_matrix::difference_type difference_type;
		public: typedef typename generalized_diagonal_matrix::const_reference reference;
		public: typedef const typename generalized_diagonal_matrix::pointer pointer;

		public: typedef const_iterator2 dual_iterator_type;
		public: typedef const_reverse_iterator2 dual_reverse_iterator_type;

		// Construction and destruction
		public: BOOST_UBLAS_INLINE
			const_iterator1()
			: container_const_reference<self_type>(),
			  it1_(),
			  it2_()
		{
		}


		public: BOOST_UBLAS_INLINE
			const_iterator1(self_type const& m, size_type it1, size_type it2)
			: container_const_reference<self_type>(m),
			  it1_(it1),
			  it2_(it2)
		{
		}


		public: BOOST_UBLAS_INLINE
			const_iterator1(const iterator1 &it)
			: container_const_reference<self_type>(it()),
			  it1_(it.it1_),
			  it2_(it.it2_)
		{
		}


		// Arithmetic
		public: BOOST_UBLAS_INLINE
			const_iterator1& operator++()
		{
			++it1_;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			const_iterator1& operator--()
		{
			--it1_;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			const_iterator1& operator+=(difference_type n)
		{
			it1_ += n;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			const_iterator1& operator-=(difference_type n)
		{
			it1_ -= n;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			difference_type operator-(const_iterator1 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it2_ == it.it2_, external_logic());
			return it1_ - it.it1_;
		}


		// Dereference
		public: BOOST_UBLAS_INLINE
			const_reference operator*() const
		{
			return (*this)()(it1_, it2_);
		}


		public: BOOST_UBLAS_INLINE
			const_reference operator[](difference_type n) const
		{
			return *(*this + n);
		}


#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			const_iterator2 begin() const
		{
			difference_type k = (*this)().offset();
			//return (*this)().find2(1, it1_, it1_ + (k > 0 ? k : 0));
			return (*this)().find2(1, it1_, it1_ + k);
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			const_iterator2 end() const
		{
			difference_type k = (*this)().offset();
			//return (*this)().find2(1, it1_, it1_ + (k > 0 ? k : 0) + 1);
			return (*this)().find2(1, it1_, it1_ + k + 1);
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			const_reverse_iterator2 rbegin() const
		{
			return const_reverse_iterator2(end());
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			const_reverse_iterator2 rend() const
		{
			return const_reverse_iterator2(begin());
		}
#endif


		// Indices
		public: BOOST_UBLAS_INLINE
			size_type index1() const
		{
			return it1_;
		}


		public: BOOST_UBLAS_INLINE
			size_type index2() const
		{
			return it2_;
		}


		// Assignment
		public: BOOST_UBLAS_INLINE
			const_iterator1& operator=(const_iterator1 const& it)
		{
			container_const_reference<self_type>::assign (&it ());
			it1_ = it.it1_;
			it2_ = it.it2_;
			return *this;
		}


		// Comparison
		public: BOOST_UBLAS_INLINE
			bool operator==(const_iterator1 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it2_ == it.it2_, external_logic());
			return it1_ == it.it1_;
		}


		public: BOOST_UBLAS_INLINE
			bool operator<(const_iterator1 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it2_ == it.it2_, external_logic());
			return it1_ < it.it1_;
		}


		private: size_type it1_;
		private: size_type it2_;
	};


	public: class iterator1: 	public container_reference<generalized_diagonal_matrix>,
								public random_access_iterator_base<packed_random_access_iterator_tag, iterator1, value_type>
	{
		public: typedef typename generalized_diagonal_matrix::value_type value_type;
		public: typedef typename generalized_diagonal_matrix::difference_type difference_type;
		public: typedef typename generalized_diagonal_matrix::reference reference;
		public: typedef typename generalized_diagonal_matrix::pointer pointer;

		public: typedef iterator2 dual_iterator_type;
		public: typedef reverse_iterator2 dual_reverse_iterator_type;


		// Construction and destruction
		public: BOOST_UBLAS_INLINE
			iterator1()
			: container_reference<self_type>(),
			  it1_(),
			  it2_()
		{
		}


		public: BOOST_UBLAS_INLINE
			iterator1(self_type& m, size_type it1, size_type it2)
			: container_reference<self_type>(m),
			  it1_(it1),
			  it2_(it2)
		{
		}


		// Arithmetic
		public: BOOST_UBLAS_INLINE
			iterator1& operator++()
		{
			++it1_;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			iterator1& operator--()
		{
			--it1_;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			iterator1& operator+=(difference_type n)
		{
			it1_ += n;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			iterator1& operator-=(difference_type n)
		{
			it1_ -= n;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			difference_type operator-(iterator1 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it2_ == it.it2_, external_logic());
			return it1_ - it.it1_;
		}


		// Dereference
		public: BOOST_UBLAS_INLINE
			reference operator*() const
		{
			return (*this)().at_element(it1_, it2_);
		}


		public: BOOST_UBLAS_INLINE
			reference operator[](difference_type n) const
		{
			return *(*this + n);
		}


#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			iterator2 begin() const
		{
			difference_type k = (*this)().offset();
			//return (*this)().find2(1, it1_, it1_ + (k > 0 ? k : 0));
			return (*this)().find2(1, it1_, it1_ + k);
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			iterator2 end() const
		{
			difference_type k = (*this)().offset();
			//return (*this)().find2(1, it1_, it1_ + (k > 0 ? k : 0) + 1);
			return (*this)().find2(1, it1_, it1_ + k + 1);
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			reverse_iterator2 rbegin() const
		{
			return reverse_iterator2(end());
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			reverse_iterator2 rend() const
		{
			return reverse_iterator2(begin());
		}
#endif


		// Indices
		public: BOOST_UBLAS_INLINE
			size_type index1() const
		{
			return it1_;
		}


		public: BOOST_UBLAS_INLINE
			size_type index2() const
		{
			return it2_;
		}


		// Assignment
		public: BOOST_UBLAS_INLINE
			iterator1& operator=(iterator1 const& it)
		{
			container_reference<self_type>::assign(&it());
			it1_ = it.it1_;
			it2_ = it.it2_;
			return *this;
		}


		// Comparison
		public: BOOST_UBLAS_INLINE
			bool operator==(iterator1 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it2_ == it.it2_, external_logic());
			return it1_ == it.it1_;
		}


		public: BOOST_UBLAS_INLINE
			bool operator<(iterator1 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it2_ == it.it2_, external_logic());
			return it1_ < it.it1_;
		}


		private: size_type it1_;
		private: size_type it2_;
		private: friend class const_iterator1;
	};


	public: class const_iterator2: 	public container_const_reference<generalized_diagonal_matrix>,
									public random_access_iterator_base<packed_random_access_iterator_tag, const_iterator2, value_type>
	{
		public: typedef typename generalized_diagonal_matrix::value_type value_type;
		public: typedef typename generalized_diagonal_matrix::difference_type difference_type;
		public: typedef typename generalized_diagonal_matrix::const_reference reference;
		public: typedef const typename generalized_diagonal_matrix::pointer pointer;

		public: typedef const_iterator1 dual_iterator_type;
		public: typedef const_reverse_iterator1 dual_reverse_iterator_type;

		// Construction and destruction
		public: BOOST_UBLAS_INLINE
			const_iterator2()
			: container_const_reference<self_type>(),
			  it1_(),
			  it2_()
		{
		}


		public: BOOST_UBLAS_INLINE
			const_iterator2(self_type const& m, size_type it1, size_type it2)
			: container_const_reference<self_type>(m),
			  it1_(it1),
			  it2_(it2)
		{
		}


		public: BOOST_UBLAS_INLINE
			const_iterator2(iterator2 const& it)
			: container_const_reference<self_type>(it()),
			  it1_(it.it1_),
			  it2_(it.it2_)
		{
		}


		// Arithmetic
		public: BOOST_UBLAS_INLINE
			const_iterator2& operator++()
		{
			++it2_;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			const_iterator2& operator--()
		{
			--it2_;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			const_iterator2& operator+=(difference_type n)
		{
			it2_ += n;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			const_iterator2& operator-=(difference_type n)
		{
			it2_ -= n;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			difference_type operator-(const_iterator2 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it1_ == it.it1_, external_logic());
			return it2_ - it.it2_;
		}


		// Dereference
		public: BOOST_UBLAS_INLINE
			const_reference operator*() const
		{
			return (*this)()(it1_, it2_);
		}


		public: BOOST_UBLAS_INLINE
			const_reference operator[](difference_type n) const
		{
			return *(*this + n);
		}


#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			const_iterator1 begin() const
		{
			difference_type k = (*this)().offset();
			//return (*this)().find1(1, it2_ - (k < 0 ? k : 0), it2_);
			return (*this)().find1(1, it2_ - k, it2_);
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			const_iterator1 end() const
		{
			difference_type k = (*this)().offset();
			//return (*this)().find1(1, it2_ - (k < 0 ? k : 0) + 1, it2_);
			return (*this)().find1(1, it2_ - k + 1, it2_);
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			const_reverse_iterator1 rbegin() const
		{
			return const_reverse_iterator1(end());
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			const_reverse_iterator1 rend() const
		{
			return const_reverse_iterator1(begin());
		}
#endif


		// Indices
		public: BOOST_UBLAS_INLINE
			size_type index1() const
		{
			return it1_;
		}


		public: BOOST_UBLAS_INLINE
			size_type index2() const
		{
			return it2_;
		}


		// Assignment
		public: BOOST_UBLAS_INLINE
			const_iterator2& operator=(const_iterator2 const& it)
		{
			container_const_reference<self_type>::assign(&it());
			it1_ = it.it1_;
			it2_ = it.it2_;

			return *this;
		}


		// Comparison
		public: BOOST_UBLAS_INLINE
			bool operator==(const_iterator2 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it1_ == it.it1_, external_logic());
			return it2_ == it.it2_;
		}


		public: BOOST_UBLAS_INLINE
			bool operator<(const_iterator2 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it1_ == it.it1_, external_logic());
			return it2_ < it.it2_;
		}


		private: size_type it1_;
		private: size_type it2_;
	};


	public: class iterator2: 	public container_reference<generalized_diagonal_matrix>,
								public random_access_iterator_base<packed_random_access_iterator_tag, iterator2, value_type>
	{
		public: typedef typename generalized_diagonal_matrix::value_type value_type;
		public: typedef typename generalized_diagonal_matrix::difference_type difference_type;
		public: typedef typename generalized_diagonal_matrix::reference reference;
		public: typedef typename generalized_diagonal_matrix::pointer pointer;

		public: typedef iterator1 dual_iterator_type;
		public: typedef reverse_iterator1 dual_reverse_iterator_type;

		// Construction and destruction
		public: BOOST_UBLAS_INLINE
			iterator2()
			: container_reference<self_type>(),
			  it1_(),
			  it2_()
		{
		}


		public: BOOST_UBLAS_INLINE
			iterator2(self_type& m, size_type it1, size_type it2)
			: container_reference<self_type>(m),
			  it1_(it1),
			  it2_(it2)
		{
		}


		// Arithmetic
		public: BOOST_UBLAS_INLINE
			iterator2& operator++()
		{
			++it2_;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			iterator2& operator--()
		{
			--it2_;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			iterator2& operator+=(difference_type n)
		{
			it2_ += n;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			iterator2& operator-=(difference_type n)
		{
			it2_ -= n;
			return *this;
		}


		public: BOOST_UBLAS_INLINE
			difference_type operator-(iterator2 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it1_ == it.it1_, external_logic());
			return it2_ - it.it2_;
		}


		// Dereference
		public: BOOST_UBLAS_INLINE
			reference operator*() const
		{
			return (*this)().at_element(it1_, it2_);
		}


		public: BOOST_UBLAS_INLINE
			reference operator[](difference_type n) const
		{
			return *(*this + n);
		}


#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			iterator1 begin() const
		{
			difference_type k = (*this)().offset();
			//return (*this)().find1(1, it2_ - (k < 0 ? k : 0), it2_);
			return (*this)().find1(1, it2_ - k, it2_);
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			iterator1 end() const
		{
			difference_type k = (*this)().offset();
			//return (*this)().find1(1, it2_ - (k < 0 ? k : 0) + 1, it2_);
			return (*this)().find1(1, it2_ - k + 1, it2_);
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			reverse_iterator1 rbegin() const
		{
			return reverse_iterator1(end());
		}


		public: BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
			typename self_type::
#endif
			reverse_iterator1 rend() const
		{
			return reverse_iterator1(begin());
		}
#endif


		// Indices
		public: BOOST_UBLAS_INLINE
			size_type index1 () const
		{
			return it1_;
		}


		public: BOOST_UBLAS_INLINE
			size_type index2() const
		{
			return it2_;
		}

		// Assignment
		public: BOOST_UBLAS_INLINE
			iterator2& operator=(iterator2 const& it)
		{
			container_reference<self_type>::assign (&it ());
			it1_ = it.it1_;
			it2_ = it.it2_;
			return *this;
		}


		// Comparison
		public: BOOST_UBLAS_INLINE
			bool operator==(iterator2 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this)() == &it(), external_logic());
			BOOST_UBLAS_CHECK(it1_ == it.it1_, external_logic());
			return it2_ == it.it2_;
		}


		public: BOOST_UBLAS_INLINE
			bool operator<(iterator2 const& it) const
		{
			BOOST_UBLAS_CHECK(&(*this) () == &it (), external_logic());
			BOOST_UBLAS_CHECK(it1_ == it.it1_, external_logic());
			return it2_ < it.it2_;
		}


		private: size_type it1_;
		private: size_type it2_;
		private: friend class const_iterator2;
	};
#endif // BOOST_UBLAS_USE_INDEXED_ITERATOR


	//@} Iterators types


	private: size_type size1_;
	private: size_type size2_;
	private: difference_type k_;
	private: size_type r_;
	private: size_type c_;
	private: array_type data_;
	private: static const_value_type zero_;
};


template<class ValueT, class LayoutT, class ArrayT>
typename generalized_diagonal_matrix<ValueT,LayoutT,ArrayT>::const_value_type generalized_diagonal_matrix<ValueT,LayoutT,ArrayT>::zero_ = generalized_diagonal_matrix<ValueT,LayoutT,ArrayT>::value_type/*zero*/();


/*TODO: the code belowe is taken from banded.hpp
// Generalized diagonal matrix adaptor class
template <typename MatrixT>
class generalized_diagonal_adaptor: public matrix_expression< generalized_diagonal_adaptor<MatrixT> >
{
	private: typedef generalized_diagonal_adaptor<MatrixT> self_type;
#ifdef BOOST_UBLAS_ENABLE_PROXY_SHORTCUTS
	public: using matrix_expression<self_type>::operator ();
#endif
	public: typedef const M const_matrix_type;
	public: typedef M matrix_type;
	public: typedef typename M::size_type size_type;
	public: typedef typename M::difference_type difference_type;
	public: typedef typename M::value_type value_type;
	public: typedef typename M::const_reference const_reference;
	public: typedef typename boost::mpl::if_<boost::is_const<M>,
									  typename M::const_reference,
									  typename M::reference>::type reference;
	public: typedef typename boost::mpl::if_<boost::is_const<M>,
									  typename M::const_closure_type,
									  typename M::closure_type>::type matrix_closure_type;
	public: typedef const self_type const_closure_type;
	public: typedef self_type closure_type;
	// Replaced by _temporary_traits to avoid type requirements on M
	//typedef typename M::vector_temporary_type vector_temporary_type;
	//typedef typename M::matrix_temporary_type matrix_temporary_type;
	public: typedef typename storage_restrict_traits<typename M::storage_category,
											 packed_proxy_tag>::storage_category storage_category;
	public: typedef typename M::orientation_category orientation_category;

	// Construction and destruction
	public: BOOST_UBLAS_INLINE
		diagonal_adaptor (matrix_type &data, size_type lower = 0, size_type upper = 0)
		: matrix_expression<self_type> (),
		  data_ (data), lower_ (lower), upper_ (upper) {}


	public: BOOST_UBLAS_INLINE
		diagonal_adaptor (const diagonal_adaptor &m)
		: matrix_expression<self_type> (),
		  data_ (m.data_), lower_ (m.lower_), upper_ (m.upper_) {}


	// Accessors
	public: BOOST_UBLAS_INLINE
	size_type size1 () const {
		return data_.size1 ();
	}
	public: BOOST_UBLAS_INLINE
	size_type size2 () const {
		return data_.size2 ();
	}
	BOOST_UBLAS_INLINE
	difference_type offset() const {
		return k_;
	}


	// Storage accessors
	BOOST_UBLAS_INLINE
	const matrix_closure_type &data () const {
		return data_;
	}
	BOOST_UBLAS_INLINE
	matrix_closure_type &data () {
		return data_;
	}

	// Element access
#ifndef BOOST_UBLAS_PROXY_CONST_MEMBER
	BOOST_UBLAS_INLINE
	const_reference operator () (size_type i, size_type j) const {
		BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
		BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
		//TODO
		return zero_;
	}
	BOOST_UBLAS_INLINE
	reference operator () (size_type i, size_type j) {
		BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
		BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
		//TODO
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
		bad_index ().raise ();
#endif
		return const_cast<reference>(zero_);
	}
#else
	BOOST_UBLAS_INLINE
	reference operator () (size_type i, size_type j) const {
		BOOST_UBLAS_CHECK (i < size1 (), bad_index ());
		BOOST_UBLAS_CHECK (j < size2 (), bad_index ());
#ifndef BOOST_UBLAS_REFERENCE_CONST_MEMBER
		bad_index ().raise ();
#endif
		return const_cast<reference>(zero_);
	}
#endif

	// Assignment
	BOOST_UBLAS_INLINE
	diagonal_adaptor &operator = (const diagonal_adaptor &m) {
		matrix_assign<scalar_assign> (*this, m);
		return *this;
	}
	BOOST_UBLAS_INLINE
	diagonal_adaptor &assign_temporary (diagonal_adaptor &m) {
		*this = m;
		return *this;
	}
	template<class AE>
	BOOST_UBLAS_INLINE
	diagonal_adaptor &operator = (const matrix_expression<AE> &ae) {
		matrix_assign<scalar_assign> (*this, matrix<value_type> (ae));
		return *this;
	}
	template<class AE>
	BOOST_UBLAS_INLINE
	diagonal_adaptor &assign (const matrix_expression<AE> &ae) {
		matrix_assign<scalar_assign> (*this, ae);
		return *this;
	}
	template<class AE>
	BOOST_UBLAS_INLINE
	diagonal_adaptor& operator += (const matrix_expression<AE> &ae) {
		matrix_assign<scalar_assign> (*this, matrix<value_type> (*this + ae));
		return *this;
	}
	template<class AE>
	BOOST_UBLAS_INLINE
	diagonal_adaptor &plus_assign (const matrix_expression<AE> &ae) {
		matrix_assign<scalar_plus_assign> (*this, ae);
		return *this;
	}
	template<class AE>
	BOOST_UBLAS_INLINE
	diagonal_adaptor& operator -= (const matrix_expression<AE> &ae) {
		matrix_assign<scalar_assign> (*this, matrix<value_type> (*this - ae));
		return *this;
	}
	template<class AE>
	BOOST_UBLAS_INLINE
	diagonal_adaptor &minus_assign (const matrix_expression<AE> &ae) {
		matrix_assign<scalar_minus_assign> (*this, ae);
		return *this;
	}
	template<class AT>
	BOOST_UBLAS_INLINE
	diagonal_adaptor& operator *= (const AT &at) {
		matrix_assign_scalar<scalar_multiplies_assign> (*this, at);
		return *this;
	}
	template<class AT>
	BOOST_UBLAS_INLINE
	diagonal_adaptor& operator /= (const AT &at) {
		matrix_assign_scalar<scalar_divides_assign> (*this, at);
		return *this;
	}

	// Closure comparison
	BOOST_UBLAS_INLINE
	bool same_closure (const diagonal_adaptor &ba) const {
		return (*this).data ().same_closure (ba.data ());
	}

	// Swapping
	BOOST_UBLAS_INLINE
	void swap (diagonal_adaptor &m) {
		if (this != &m) {
			BOOST_UBLAS_CHECK (lower_ == m.lower_, bad_size ());
			BOOST_UBLAS_CHECK (upper_ == m.upper_, bad_size ());
			matrix_swap<scalar_swap> (*this, m);
		}
	}
	BOOST_UBLAS_INLINE
	friend void swap (diagonal_adaptor &m1, diagonal_adaptor &m2) {
		m1.swap (m2);
	}

	// Iterator types
private:
	// Use the matrix iterator
	typedef typename M::const_iterator1 const_subiterator1_type;
	typedef typename boost::mpl::if_<boost::is_const<M>,
									  typename M::const_iterator1,
									  typename M::iterator1>::type subiterator1_type;
	typedef typename M::const_iterator2 const_subiterator2_type;
	typedef typename boost::mpl::if_<boost::is_const<M>,
									  typename M::const_iterator2,
									  typename M::iterator2>::type subiterator2_type;

public:
#ifdef BOOST_UBLAS_USE_INDEXED_ITERATOR
	typedef indexed_iterator1<self_type, packed_random_access_iterator_tag> iterator1;
	typedef indexed_iterator2<self_type, packed_random_access_iterator_tag> iterator2;
	typedef indexed_const_iterator1<self_type, packed_random_access_iterator_tag> const_iterator1;
	typedef indexed_const_iterator2<self_type, packed_random_access_iterator_tag> const_iterator2;
#else
	class const_iterator1;
	class iterator1;
	class const_iterator2;
	class iterator2;
#endif
	typedef reverse_iterator_base1<const_iterator1> const_reverse_iterator1;
	typedef reverse_iterator_base1<iterator1> reverse_iterator1;
	typedef reverse_iterator_base2<const_iterator2> const_reverse_iterator2;
	typedef reverse_iterator_base2<iterator2> reverse_iterator2;

	// Element lookup
	BOOST_UBLAS_INLINE
	const_iterator1 find1 (int rank, size_type i, size_type j) const {
		if (rank == 1) {
			size_type lower_i = (std::max) (difference_type (j - upper_), difference_type (0));
			i = (std::max) (i, lower_i);
			size_type upper_i = (std::min) (j + 1 + lower_, size1 ());
			i = (std::min) (i, upper_i);
		}
		return const_iterator1 (*this, data ().find1 (rank, i, j));
	}
	BOOST_UBLAS_INLINE
	iterator1 find1 (int rank, size_type i, size_type j) {
		if (rank == 1) {
			size_type lower_i = (std::max) (difference_type (j - upper_), difference_type (0));
			i = (std::max) (i, lower_i);
			size_type upper_i = (std::min) (j + 1 + lower_, size1 ());
			i = (std::min) (i, upper_i);
		}
		return iterator1 (*this, data ().find1 (rank, i, j));
	}
	BOOST_UBLAS_INLINE
	const_iterator2 find2 (int rank, size_type i, size_type j) const {
		if (rank == 1) {
			size_type lower_j = (std::max) (difference_type (i - lower_), difference_type (0));
			j = (std::max) (j, lower_j);
			size_type upper_j = (std::min) (i + 1 + upper_, size2 ());
			j = (std::min) (j, upper_j);
		}
		return const_iterator2 (*this, data ().find2 (rank, i, j));
	}
	BOOST_UBLAS_INLINE
	iterator2 find2 (int rank, size_type i, size_type j) {
		if (rank == 1) {
			size_type lower_j = (std::max) (difference_type (i - lower_), difference_type (0));
			j = (std::max) (j, lower_j);
			size_type upper_j = (std::min) (i + 1 + upper_, size2 ());
			j = (std::min) (j, upper_j);
		}
		return iterator2 (*this, data ().find2 (rank, i, j));
	}

	// Iterators simply are indices.

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
	class const_iterator1:
		public container_const_reference<diagonal_adaptor>,
		public random_access_iterator_base<typename iterator_restrict_traits<
											   typename const_subiterator1_type::iterator_category, packed_random_access_iterator_tag>::iterator_category,
										   const_iterator1, value_type> {
	public:
		typedef typename const_subiterator1_type::value_type value_type;
		typedef typename const_subiterator1_type::difference_type difference_type;
		typedef typename const_subiterator1_type::reference reference;
		typedef typename const_subiterator1_type::pointer pointer;

		typedef const_iterator2 dual_iterator_type;
		typedef const_reverse_iterator2 dual_reverse_iterator_type;

		// Construction and destruction
		BOOST_UBLAS_INLINE
		const_iterator1 ():
			container_const_reference<self_type> (), it1_ () {}
		BOOST_UBLAS_INLINE
		const_iterator1 (const self_type &m, const const_subiterator1_type &it1):
			container_const_reference<self_type> (m), it1_ (it1) {}
		BOOST_UBLAS_INLINE
		const_iterator1 (const iterator1 &it):
			container_const_reference<self_type> (it ()), it1_ (it.it1_) {}

		// Arithmetic
		BOOST_UBLAS_INLINE
		const_iterator1 &operator ++ () {
			++ it1_;
			return *this;
		}
		BOOST_UBLAS_INLINE
		const_iterator1 &operator -- () {
			-- it1_;
			return *this;
		}
		BOOST_UBLAS_INLINE
		const_iterator1 &operator += (difference_type n) {
			it1_ += n;
			return *this;
		}
		BOOST_UBLAS_INLINE
		const_iterator1 &operator -= (difference_type n) {
			it1_ -= n;
			return *this;
		}
		BOOST_UBLAS_INLINE
		difference_type operator - (const const_iterator1 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it1_ - it.it1_;
		}

		// Dereference
		BOOST_UBLAS_INLINE
		const_reference operator * () const {
			size_type i = index1 ();
			size_type j = index2 ();
			BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
			BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
			size_type k = j;
			size_type l = (*this) ().upper () + i - j;
			if (k < (*this) ().size2 () &&
				l < (*this) ().lower () + 1 + (*this) ().upper ())
				return *it1_;
			return (*this) () (i, j);
		}
		BOOST_UBLAS_INLINE
		const_reference operator [] (difference_type n) const {
			return *(*this + n);
		}

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		const_iterator2 begin () const {
			return (*this) ().find2 (1, index1 (), 0);
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		const_iterator2 end () const {
			return (*this) ().find2 (1, index1 (), (*this) ().size2 ());
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		const_reverse_iterator2 rbegin () const {
			return const_reverse_iterator2 (end ());
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		const_reverse_iterator2 rend () const {
			return const_reverse_iterator2 (begin ());
		}
#endif

		// Indices
		BOOST_UBLAS_INLINE
		size_type index1 () const {
			return it1_.index1 ();
		}
		BOOST_UBLAS_INLINE
		size_type index2 () const {
			return it1_.index2 ();
		}

		// Assignment
		BOOST_UBLAS_INLINE
		const_iterator1 &operator = (const const_iterator1 &it) {
			container_const_reference<self_type>::assign (&it ());
			it1_ = it.it1_;
			return *this;
		}

		// Comparison
		BOOST_UBLAS_INLINE
		bool operator == (const const_iterator1 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it1_ == it.it1_;
		}
		BOOST_UBLAS_INLINE
		bool operator < (const const_iterator1 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it1_ < it.it1_;
		}

	private:
		const_subiterator1_type it1_;
	};
#endif

	BOOST_UBLAS_INLINE
	const_iterator1 begin1 () const {
		return find1 (0, 0, 0);
	}
	BOOST_UBLAS_INLINE
	const_iterator1 end1 () const {
		return find1 (0, size1 (), 0);
	}

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
	class iterator1:
		public container_reference<diagonal_adaptor>,
		public random_access_iterator_base<typename iterator_restrict_traits<
											   typename subiterator1_type::iterator_category, packed_random_access_iterator_tag>::iterator_category,
										   iterator1, value_type> {
	public:
		typedef typename subiterator1_type::value_type value_type;
		typedef typename subiterator1_type::difference_type difference_type;
		typedef typename subiterator1_type::reference reference;
		typedef typename subiterator1_type::pointer pointer;

		typedef iterator2 dual_iterator_type;
		typedef reverse_iterator2 dual_reverse_iterator_type;

		// Construction and destruction
		BOOST_UBLAS_INLINE
		iterator1 ():
			container_reference<self_type> (), it1_ () {}
		BOOST_UBLAS_INLINE
		iterator1 (self_type &m, const subiterator1_type &it1):
			container_reference<self_type> (m), it1_ (it1) {}

		// Arithmetic
		BOOST_UBLAS_INLINE
		iterator1 &operator ++ () {
			++ it1_;
			return *this;
		}
		BOOST_UBLAS_INLINE
		iterator1 &operator -- () {
			-- it1_;
			return *this;
		}
		BOOST_UBLAS_INLINE
		iterator1 &operator += (difference_type n) {
			it1_ += n;
			return *this;
		}
		BOOST_UBLAS_INLINE
		iterator1 &operator -= (difference_type n) {
			it1_ -= n;
			return *this;
		}
		BOOST_UBLAS_INLINE
		difference_type operator - (const iterator1 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it1_ - it.it1_;
		}

		// Dereference
		BOOST_UBLAS_INLINE
		reference operator * () const {
			size_type i = index1 ();
			size_type j = index2 ();
			BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
			BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
			size_type k = j;
			size_type l = (*this) ().upper () + i - j;
			if (k < (*this) ().size2 () &&
				l < (*this) ().lower () + 1 + (*this) ().upper ())
				return *it1_;
			return (*this) () (i, j);
		}
		BOOST_UBLAS_INLINE
		reference operator [] (difference_type n) const {
			return *(*this + n);
		}

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		iterator2 begin () const {
			return (*this) ().find2 (1, index1 (), 0);
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		iterator2 end () const {
			return (*this) ().find2 (1, index1 (), (*this) ().size2 ());
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		reverse_iterator2 rbegin () const {
			return reverse_iterator2 (end ());
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		reverse_iterator2 rend () const {
			return reverse_iterator2 (begin ());
		}
#endif

		// Indices
		BOOST_UBLAS_INLINE
		size_type index1 () const {
			return it1_.index1 ();
		}
		BOOST_UBLAS_INLINE
		size_type index2 () const {
			return it1_.index2 ();
		}

		// Assignment
		BOOST_UBLAS_INLINE
		iterator1 &operator = (const iterator1 &it) {
			container_reference<self_type>::assign (&it ());
			it1_ = it.it1_;
			return *this;
		}

		// Comparison
		BOOST_UBLAS_INLINE
		bool operator == (const iterator1 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it1_ == it.it1_;
		}
		BOOST_UBLAS_INLINE
		bool operator < (const iterator1 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it1_ < it.it1_;
		}

	private:
		subiterator1_type it1_;

		friend class const_iterator1;
	};
#endif

	BOOST_UBLAS_INLINE
	iterator1 begin1 () {
		return find1 (0, 0, 0);
	}
	BOOST_UBLAS_INLINE
	iterator1 end1 () {
		return find1 (0, size1 (), 0);
	}

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
	class const_iterator2:
		public container_const_reference<diagonal_adaptor>,
		public random_access_iterator_base<packed_random_access_iterator_tag,
										   const_iterator2, value_type> {
	public:
		typedef typename iterator_restrict_traits<typename const_subiterator2_type::iterator_category,
												  packed_random_access_iterator_tag>::iterator_category iterator_category;
		typedef typename const_subiterator2_type::value_type value_type;
		typedef typename const_subiterator2_type::difference_type difference_type;
		typedef typename const_subiterator2_type::reference reference;
		typedef typename const_subiterator2_type::pointer pointer;

		typedef const_iterator1 dual_iterator_type;
		typedef const_reverse_iterator1 dual_reverse_iterator_type;

		// Construction and destruction
		BOOST_UBLAS_INLINE
		const_iterator2 ():
			container_const_reference<self_type> (), it2_ () {}
		BOOST_UBLAS_INLINE
		const_iterator2 (const self_type &m, const const_subiterator2_type &it2):
			container_const_reference<self_type> (m), it2_ (it2) {}
		BOOST_UBLAS_INLINE
		const_iterator2 (const iterator2 &it):
			container_const_reference<self_type> (it ()), it2_ (it.it2_) {}

		// Arithmetic
		BOOST_UBLAS_INLINE
		const_iterator2 &operator ++ () {
			++ it2_;
			return *this;
		}
		BOOST_UBLAS_INLINE
		const_iterator2 &operator -- () {
			-- it2_;
			return *this;
		}
		BOOST_UBLAS_INLINE
		const_iterator2 &operator += (difference_type n) {
			it2_ += n;
			return *this;
		}
		BOOST_UBLAS_INLINE
		const_iterator2 &operator -= (difference_type n) {
			it2_ -= n;
			return *this;
		}
		BOOST_UBLAS_INLINE
		difference_type operator - (const const_iterator2 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it2_ - it.it2_;
		}

		// Dereference
		BOOST_UBLAS_INLINE
		const_reference operator * () const {
			size_type i = index1 ();
			size_type j = index2 ();
			BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
			BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
			size_type k = j;
			size_type l = (*this) ().upper () + i - j;
			if (k < (*this) ().size2 () &&
				l < (*this) ().lower () + 1 + (*this) ().upper ())
				return *it2_;
			return (*this) () (i, j);
		}
		BOOST_UBLAS_INLINE
		const_reference operator [] (difference_type n) const {
			return *(*this + n);
		}

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		const_iterator1 begin () const {
			return (*this) ().find1 (1, 0, index2 ());
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		const_iterator1 end () const {
			return (*this) ().find1 (1, (*this) ().size1 (), index2 ());
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		const_reverse_iterator1 rbegin () const {
			return const_reverse_iterator1 (end ());
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		const_reverse_iterator1 rend () const {
			return const_reverse_iterator1 (begin ());
		}
#endif

		// Indices
		BOOST_UBLAS_INLINE
		size_type index1 () const {
			return it2_.index1 ();
		}
		BOOST_UBLAS_INLINE
		size_type index2 () const {
			return it2_.index2 ();
		}

		// Assignment
		BOOST_UBLAS_INLINE
		const_iterator2 &operator = (const const_iterator2 &it) {
			container_const_reference<self_type>::assign (&it ());
			it2_ = it.it2_;
			return *this;
		}

		// Comparison
		BOOST_UBLAS_INLINE
		bool operator == (const const_iterator2 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it2_ == it.it2_;
		}
		BOOST_UBLAS_INLINE
		bool operator < (const const_iterator2 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it2_ < it.it2_;
		}

	private:
		const_subiterator2_type it2_;
	};
#endif

	BOOST_UBLAS_INLINE
	const_iterator2 begin2 () const {
		return find2 (0, 0, 0);
	}
	BOOST_UBLAS_INLINE
	const_iterator2 end2 () const {
		return find2 (0, 0, size2 ());
	}

#ifndef BOOST_UBLAS_USE_INDEXED_ITERATOR
	class iterator2:
		public container_reference<diagonal_adaptor>,
		public random_access_iterator_base<typename iterator_restrict_traits<
											   typename subiterator2_type::iterator_category, packed_random_access_iterator_tag>::iterator_category,
										   iterator2, value_type> {
	public:
		typedef typename subiterator2_type::value_type value_type;
		typedef typename subiterator2_type::difference_type difference_type;
		typedef typename subiterator2_type::reference reference;
		typedef typename subiterator2_type::pointer pointer;

		typedef iterator1 dual_iterator_type;
		typedef reverse_iterator1 dual_reverse_iterator_type;

		// Construction and destruction
		BOOST_UBLAS_INLINE
		iterator2 ():
			container_reference<self_type> (), it2_ () {}
		BOOST_UBLAS_INLINE
		iterator2 (self_type &m, const subiterator2_type &it2):
			container_reference<self_type> (m), it2_ (it2) {}

		// Arithmetic
		BOOST_UBLAS_INLINE
		iterator2 &operator ++ () {
			++ it2_;
			return *this;
		}
		BOOST_UBLAS_INLINE
		iterator2 &operator -- () {
			-- it2_;
			return *this;
		}
		BOOST_UBLAS_INLINE
		iterator2 &operator += (difference_type n) {
			it2_ += n;
			return *this;
		}
		BOOST_UBLAS_INLINE
		iterator2 &operator -= (difference_type n) {
			it2_ -= n;
			return *this;
		}
		BOOST_UBLAS_INLINE
		difference_type operator - (const iterator2 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it2_ - it.it2_;
		}

		// Dereference
		BOOST_UBLAS_INLINE
		reference operator * () const {
			size_type i = index1 ();
			size_type j = index2 ();
			BOOST_UBLAS_CHECK (i < (*this) ().size1 (), bad_index ());
			BOOST_UBLAS_CHECK (j < (*this) ().size2 (), bad_index ());
			size_type k = j;
			size_type l = (*this) ().upper () + i - j;
			if (k < (*this) ().size2 () &&
				l < (*this) ().lower () + 1 + (*this) ().upper ())
				return *it2_;
			return (*this) () (i, j);
		}
		BOOST_UBLAS_INLINE
		reference operator [] (difference_type n) const {
			return *(*this + n);
		}

#ifndef BOOST_UBLAS_NO_NESTED_CLASS_RELATION
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		iterator1 begin () const {
			return (*this) ().find1 (1, 0, index2 ());
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		iterator1 end () const {
			return (*this) ().find1 (1, (*this) ().size1 (), index2 ());
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		reverse_iterator1 rbegin () const {
			return reverse_iterator1 (end ());
		}
		BOOST_UBLAS_INLINE
#ifdef BOOST_UBLAS_MSVC_NESTED_CLASS_RELATION
		typename self_type::
#endif
		reverse_iterator1 rend () const {
			return reverse_iterator1 (begin ());
		}
#endif

		// Indices
		BOOST_UBLAS_INLINE
		size_type index1 () const {
			return it2_.index1 ();
		}
		BOOST_UBLAS_INLINE
		size_type index2 () const {
			return it2_.index2 ();
		}

		// Assignment
		BOOST_UBLAS_INLINE
		iterator2 &operator = (const iterator2 &it) {
			container_reference<self_type>::assign (&it ());
			it2_ = it.it2_;
			return *this;
		}

		// Comparison
		BOOST_UBLAS_INLINE
		bool operator == (const iterator2 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it2_ == it.it2_;
		}
		BOOST_UBLAS_INLINE
		bool operator < (const iterator2 &it) const {
			BOOST_UBLAS_CHECK (&(*this) () == &it (), external_logic ());
			return it2_ < it.it2_;
		}

	private:
		subiterator2_type it2_;

		friend class const_iterator2;
	};
#endif

	BOOST_UBLAS_INLINE
	iterator2 begin2 () {
		return find2 (0, 0, 0);
	}
	BOOST_UBLAS_INLINE
	iterator2 end2 () {
		return find2 (0, 0, size2 ());
	}

	// Reverse iterators

	BOOST_UBLAS_INLINE
	const_reverse_iterator1 rbegin1 () const {
		return const_reverse_iterator1 (end1 ());
	}
	BOOST_UBLAS_INLINE
	const_reverse_iterator1 rend1 () const {
		return const_reverse_iterator1 (begin1 ());
	}

	BOOST_UBLAS_INLINE
	reverse_iterator1 rbegin1 () {
		return reverse_iterator1 (end1 ());
	}
	BOOST_UBLAS_INLINE
	reverse_iterator1 rend1 () {
		return reverse_iterator1 (begin1 ());
	}

	BOOST_UBLAS_INLINE
	const_reverse_iterator2 rbegin2 () const {
		return const_reverse_iterator2 (end2 ());
	}
	BOOST_UBLAS_INLINE
	const_reverse_iterator2 rend2 () const {
		return const_reverse_iterator2 (begin2 ());
	}

	BOOST_UBLAS_INLINE
	reverse_iterator2 rbegin2 () {
		return reverse_iterator2 (end2 ());
	}
	BOOST_UBLAS_INLINE
	reverse_iterator2 rend2 () {
		return reverse_iterator2 (begin2 ());
	}

	private: matrix_closure_type data_;
	private: size_type k_;
	public: typedef const value_type const_value_type;
	public: static const_value_type zero_;
};


// Specialization for temporary_traits
template <typename M>
struct vector_temporary_traits< generalized_diagonal_adaptor<M> >
: vector_temporary_traits< M > {} ;

template <typename M>
struct vector_temporary_traits< const generalized_diagonal_adaptor<M> >
: vector_temporary_traits< M > {} ;

template <typename M>
struct matrix_temporary_traits< generalized_diagonal_adaptor<M> >
: matrix_temporary_traits< M > {} ;

template <typename M>
struct matrix_temporary_traits< const generalized_diagonal_adaptor<M> >
: matrix_temporary_traits< M > {} ;
*/

}}} // Namespace boost::numeric::ublas


#endif // BOOST_NUMERIC_UBLAS_CONTAINER_DIAGONAL_MATRIX_HPP
