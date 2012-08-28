/**
 * -*- c++ -*-
 *
 * \file size.hpp
 *
 * \brief The \c size operation.
 *
 * Copyright (c) 2009-2010, Marco Guazzone
 *
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * \author Marco Guazzone, marco.guazzone@gmail.com
 */

#ifndef BOOST_NUMERIC_UBLAS_OPERATION_SIZE_HPP
#define BOOST_NUMERIC_UBLAS_OPERATION_SIZE_HPP


#include <boost/numeric/ublas/detail/config.hpp>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/numeric/ublas/fwd.hpp>
#include <boost/numeric/ublas/tags.hpp>
#include <boost/numeric/ublas/traits.hpp>
#include <cstddef>


namespace boost { namespace numeric { namespace ublas {


//@{ Declarations

/**
 * \brief Return the number of columns.
 * \tparam VectorExprT A type which models the matrix expression concept.
 * \param m A matrix expression.
 * \return The number of columns.
 */
template <typename VectorExprT>
typename vector_traits<VectorExprT>::size_type size(vector_expression<VectorExprT> const& ve);

/**
 * \brief Return the size of the given dimension for the given expression.
 * \tparam Dim The dimension number (starting from 1).
 * \tparam VectorExprT An expression type.
 * \param e An expression.
 * \return The number of columns.
 * \return The size associated to the dimension \a Dim.
 */
template <std::size_t Dim, typename VectorExprT>
typename vector_traits<VectorExprT>::size_type size(vector_expression<VectorExprT> const& ve);

/**
 * \brief Return the size of the given dimension for the given expression.
 * \tparam Dim The dimension number (starting from 1).
 * \tparam MatrixExprT An expression type.
 * \param e An expression.
 * \return The number of columns.
 * \return The size associated to the dimension \a Dim.
 */
template <std::size_t Dim, typename MatrixExprT>
typename matrix_traits<MatrixExprT>::size_type size(matrix_expression<MatrixExprT> const& me);

/**
 * \brief Return the size of the given dimension tag for the given expression.
 * \tparam TagT The dimension tag type (e.g., tag::major).
 * \tparam MatrixExprT An expression type.
 * \param e An expression.
 * \return The size associated to the dimension tag \a TagT.
 */
template <typename TagT, typename MatrixExprT>
typename matrix_traits<MatrixExprT>::size_type size_by_tag(matrix_expression<MatrixExprT> const& me);
//typename matrix_traits<MatrixExprT>::size_type size(matrix_expression<MatrixExprT> const& me); //FIXME: don't work

//@} Declarations


//@{ Definitions

namespace detail {

/**
 * \brief Auxiliary class for computing the size of the given dimension for
 *  a container of the given category.
 * \tparam Dim The dimension number (starting from 1).
 * \tparam CategoryT The category type (e.g., vector_tag).
 */
template <std::size_t Dim, typename CategoryT>
struct size_by_dim_impl;


/**
 * \brief Auxiliary class for computing the size of the given dimension for
 *  a container of the given category and with the given orientation.
 * \tparam Dim The dimension number (starting from 1).
 * \tparam CategoryT The category type (e.g., vector_tag).
 * \tparam OrientationT The orientation category type (e.g., row_major_tag).
 */
template <typename TagT, typename CategoryT, typename OrientationT>
struct size_by_tag_impl;


/// \brief Specialization of \c size_by_dim_impl for computing the size of a
///  vector
template <>
struct size_by_dim_impl<1, vector_tag>
{
	/**
	 * \brief Compute the size of the given vector.
	 * \tparam ExprT A vector expression type.
	 * \pre ExprT must be a model of VectorExpression.
	 */
	template <typename ExprT>
	BOOST_UBLAS_INLINE
	static typename vector_traits<ExprT>::size_type apply(vector_expression<ExprT> const& ve)
	{
		return ve().size();
	}
};


/// \brief Specialization of \c size_by_dim_impl for computing the number of
///  rows of a matrix
template <>
struct size_by_dim_impl<1, matrix_tag>
{
	/**
	 * \brief Compute the number of rows of the given matrix.
	 * \tparam ExprT A matrix expression type.
	 * \pre ExprT must be a model of MatrixExpression.
	 */
	template <typename ExprT>
	BOOST_UBLAS_INLINE
	static typename matrix_traits<ExprT>::size_type apply(matrix_expression<ExprT> const& me)
	{
		return me().size1();
	}
};


/// \brief Specialization of \c size_by_dim_impl for computing the number of
///  columns of a matrix
template <>
struct size_by_dim_impl<2, matrix_tag>
{
	/**
	 * \brief Compute the number of columns of the given matrix.
	 * \tparam ExprT A matrix expression type.
	 * \pre ExprT must be a model of MatrixExpression.
	 */
	template <typename ExprT>
	BOOST_UBLAS_INLINE
	static typename matrix_traits<ExprT>::size_type apply(matrix_expression<ExprT> const& me)
	{
		return me().size2();
	}
};


/// \brief Specialization of \c size_by_tag_impl for computing the size of the
///  major dimension of a row-major oriented matrix.
template <>
struct size_by_tag_impl<tag::major, matrix_tag, row_major_tag>
{
	/**
	 * \brief Compute the number of rows of the given matrix.
	 * \tparam ExprT A matrix expression type.
	 * \pre ExprT must be a model of MatrixExpression.
	 */
	template <typename ExprT>
	BOOST_UBLAS_INLINE
	static typename matrix_traits<ExprT>::size_type apply(matrix_expression<ExprT> const& me)
	{
		return me().size1();
	}
};


/// \brief Specialization of \c size_by_tag_impl for computing the size of the
///  minor dimension of a row-major oriented matrix.
template <>
struct size_by_tag_impl<tag::minor, matrix_tag, row_major_tag>
{
	/**
	 * \brief Compute the number of columns of the given matrix.
	 * \tparam ExprT A matrix expression type.
	 * \pre ExprT must be a model of MatrixExpression.
	 */
	template <typename ExprT>
	BOOST_UBLAS_INLINE
	static typename matrix_traits<ExprT>::size_type apply(matrix_expression<ExprT> const& me)
	{
		return me().size2();
	}
};


/// \brief Specialization of \c size_by_tag_impl for computing the size of the
///  leading dimension of a row-major oriented matrix.
template <>
struct size_by_tag_impl<tag::leading, matrix_tag, row_major_tag>
{
	/**
	 * \brief Compute the number of columns of the given matrix.
	 * \tparam ExprT A matrix expression type.
	 * \pre ExprT must be a model of MatrixExpression.
	 */
	template <typename ExprT>
	BOOST_UBLAS_INLINE
	static typename matrix_traits<ExprT>::size_type apply(matrix_expression<ExprT> const& me)
	{
		return me().size2();
	}
};


/// \brief Specialization of \c size_by_tag_impl for computing the size of the
///  major dimension of a column-major oriented matrix.
template <>
struct size_by_tag_impl<tag::major, matrix_tag, column_major_tag>
{
	/**
	 * \brief Compute the number of columns of the given matrix.
	 * \tparam ExprT A matrix expression type.
	 * \pre ExprT must be a model of MatrixExpression.
	 */
	template <typename ExprT>
	BOOST_UBLAS_INLINE
	static typename matrix_traits<ExprT>::size_type apply(matrix_expression<ExprT> const& me)
	{
		return me().size2();
	}
};


/// \brief Specialization of \c size_by_tag_impl for computing the size of the
///  minor dimension of a column-major oriented matrix.
template <>
struct size_by_tag_impl<tag::minor, matrix_tag, column_major_tag>
{
	/**
	 * \brief Compute the number of rows of the given matrix.
	 * \tparam ExprT A matrix expression type.
	 * \pre ExprT must be a model of MatrixExpression.
	 */
	template <typename ExprT>
	BOOST_UBLAS_INLINE
	static typename matrix_traits<ExprT>::size_type apply(matrix_expression<ExprT> const& me)
	{
		return me().size1();
	}
};


/// \brief Specialization of \c size_by_tag_impl for computing the size of the
///  leading dimension of a column-major oriented matrix.
template <>
struct size_by_tag_impl<tag::leading, matrix_tag, column_major_tag>
{
	/**
	 * \brief Compute the number of rows of the given matrix.
	 * \tparam ExprT A matrix expression type.
	 * \pre ExprT must be a model of MatrixExpression.
	 */
	template <typename ExprT>
	BOOST_UBLAS_INLINE
	static typename matrix_traits<ExprT>::size_type apply(matrix_expression<ExprT> const& me)
	{
		return me().size1();
	}
};


/// \brief Specialization of \c size_by_tag_impl for computing the size of the
///  given dimension of a unknown oriented expression.
template <typename TagT, typename CategoryT>
struct size_by_tag_impl<TagT, CategoryT, unknown_orientation_tag>: size_by_tag_impl<TagT, CategoryT, row_major_tag>
{
	// Empty
};


} // Namespace detail


template <typename VectorExprT>
BOOST_UBLAS_INLINE
typename vector_traits<VectorExprT>::size_type size(vector_expression<VectorExprT> const& ve)
{
	return ve().size();
}


template <std::size_t Dim, typename VectorExprT>
BOOST_UBLAS_INLINE
typename vector_traits<VectorExprT>::size_type size(vector_expression<VectorExprT> const& ve)
{
	return detail::size_by_dim_impl<Dim, vector_tag>::template apply(ve);
}


template <std::size_t Dim, typename MatrixExprT>
BOOST_UBLAS_INLINE
typename matrix_traits<MatrixExprT>::size_type size(matrix_expression<MatrixExprT> const& me)
{
	return detail::size_by_dim_impl<Dim, matrix_tag>::template apply(me);
}


//[FIXME]: don't work
///**
// * \brief Return the size of the given dimension tag for the given expression.
// * \tparam TagT The dimension tag type (e.g., tag::major).
// * \tparam MatrixExprT An expression type.
// * \param e An expression.
// * \return The size associated to the dimension tag \a TagT.
// */
//template <typename TagT, typename MatrixExprT>
//BOOST_UBLAS_INLINE
//typename matrix_traits<MatrixExprT>::size_type size(matrix_expression<MatrixExprT> const& me)
//{
//	return detail::size_by_tag_impl<TagT, matrix_tag, typename matrix_traits<MatrixExprT>::orientation_category>::template apply(me);
//}
//[/FIXME]


template <typename TagT, typename MatrixExprT>
BOOST_UBLAS_INLINE
typename matrix_traits<MatrixExprT>::size_type size_by_tag(matrix_expression<MatrixExprT> const& me)
{
	return detail::size_by_tag_impl<TagT, matrix_tag, typename matrix_traits<MatrixExprT>::orientation_category>::template apply(me);
}

//@} Definitions

}}} // Namespace boost::numeric::ublas


#endif // BOOST_NUMERIC_UBLAS_OPERATION_SIZE_HPP
