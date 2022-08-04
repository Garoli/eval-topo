/**
*  This program is free software: you can redistribute it and/or modify
*  it under the terms of the GNU Lesser General Public License as
*  published by the Free Software Foundation, either version 3 of the
*  License, or  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
**/

#pragma once
/**
 * @file AffineTransformation.h
 * @author Sarah Brood (\c sarahbrood\@ecole.ensicaen.fr)
 * @author Heithem Dridi (\c heithemdridi\@ecole.ensicaen.fr)
 *
 * @date 2022/01
 *
 * This file is part of the project : Topological evaluation of image transformations
 */

#if defined(AffineTransformation_RECURSES)
#error Recursive header files inclusion detected in AffineTransformation.h
#else
#define AffineTransformation_RECURSES

#if !defined AffineTransformation_h
#define AffineTransformation_h

#include <iostream>
#include <cmath>
#include <climits>
#include <utility>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/kernel/domains/CDomain.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/math/linalg/SimpleMatrix.h>

namespace DGtal{
    namespace functors {
        template<typename TSpace>

        /**
        * Description of template class 'AffineTransformation' <p>
        * \brief Aim: Transform pixel
        * from a point to another point within a chosen transformation
        */
        class AffineTransformation{
        public:
            typedef typename TSpace::RealPoint RealPoint;
            typedef typename TSpace::Point Point;
        public:

            /**
             * Compute the new coordinates of the given the point with backward transformation
             * @param input point to transform
             */
            virtual RealPoint backward(const Point &input) const = 0;

            /**
             * Compute the new coordinates of the given the point with forward transformation
             * @param input point to transform
             */
            virtual RealPoint forward(const Point &input) const = 0;

            /**
             * Compute the new coordinates of the given the point with backward transformation
             * /!\ Backward transformation /!\
             * @param input point to transform
             */
            RealPoint operator()(const Point &input) const{
                return backward(input);
            }
        };

        /**
        * Description of class 'AffineTransformation2D' <p>
        * \brief Aim: Transform pixel, 2D implementation
        * from a point to another point within a chosen transformation
        *
        * @code

        typedef AffineTransformation2D<Space> AffineTransformation;
        AffineTransformation transformation = AffineTransformation::rotation(RealPoint(1, 0), M_PI_4);

        * @endcode
        *
        * @tparam TSpace the type of the 2d space, should be TSpace::dimension == 2
        */
        template<typename TSpace>
        class AffineTransformation2D : public AffineTransformation<TSpace>{
            BOOST_CONCEPT_ASSERT((concepts::CSpace<TSpace>));
            BOOST_STATIC_ASSERT((TSpace::dimension == 2));
        public:
            typedef typename AffineTransformation<TSpace>::RealPoint RealPoint;
            typedef typename AffineTransformation<TSpace>::Point Point;
            typedef SimpleMatrix<double, TSpace::dimension, TSpace::dimension> LinearMatrix;
            typedef SimpleMatrix<double, 1, TSpace::dimension> TranslateMatrix;
            typedef AffineTransformation2D<TSpace> Transformation;
        public:
            /**
             * Constructor.
             * @param linear part of the transformation matrix :  should be 2x2
             * @param translate  part of the transfromation : shoulb be 2x1
             */
            AffineTransformation2D(const LinearMatrix &linear, const TranslateMatrix &translate) :
                    _linear(linear), _inverse(linear.inverse()), _translate(translate){}

            /**
             * Constructor.
             * @param trans transformation to use
             */
            AffineTransformation2D(const Transformation &trans) :
                    _linear(trans._linear), _inverse(trans._inverse), _translate(trans._translate){}

            /**
             * Rotation operator
             * @param origin of the rotation
             * @param angle of the rotation
             * @return Affine transformation corresponding to rotation
             */
            static Transformation rotation(const RealPoint &origin, const double &angle){
                double cos = std::cos(angle);
                double sin = std::sin(angle);

                return Transformation(LinearMatrix({cos, -sin, sin, cos}), TranslateMatrix({0, 0}));
            }
            /**
             * Backward operator, apply the transformation backward
             * @param input
             * @return Point with transformed coordinates
             */
            RealPoint backward(const Point &input) const{
                RealPoint p;

                double a = input[0] - _translate(0, 0);
                double b = input[1] - _translate(0, 1);

                p[0] = _inverse(0, 0) * a + _inverse(0, 1) * b;
                p[1] = _inverse(1, 0) * a + _inverse(1, 1) * b;
                return p;
            }
            /**
              * Backward operator, apply the transformation forward
              * @param input
              * @return Point with transformed coordinates
              */
            RealPoint forward(const Point &input) const{
                RealPoint p;

                p[0] = _linear(0, 0) * input[0] + _linear(0, 1) * input[1] + _translate(0, 0);
                p[1] = _linear(1, 0) * input[0] + _linear(1, 1) * input[1] + _translate(0, 1);
                return p;
            }

        protected:
            LinearMatrix _linear;
            LinearMatrix _inverse;
            TranslateMatrix _translate;
        };

        /**
        * Description of class 'AffineTransformation3D' <p>
        * \brief Aim: Transform pixel, 3D implementation
        * from a point to another point within a chosen transformation
        *
        * @code

        typedef AffineTransformation3D<Space> AffineTransformation;
        AffineTransformation transformation = AffineTransformation::rotation(RealPoint(0, 0, 0), RealPoint(1, 0, 0), M_PI_4);

        * @endcode
        *
        * @tparam TSpace the type of the 3d space, should be TSpace::dimension == 3
        */
        template<typename TSpace>
        class AffineTransformation3D : public AffineTransformation<TSpace>{
            BOOST_CONCEPT_ASSERT((concepts::CSpace<TSpace>));
            BOOST_STATIC_ASSERT((TSpace::dimension == 3));
        public:
            typedef typename AffineTransformation<TSpace>::RealPoint RealPoint;
            typedef typename AffineTransformation<TSpace>::Point Point;
            typedef SimpleMatrix<double, TSpace::dimension, TSpace::dimension> LinearMatrix;
            typedef SimpleMatrix<double, 1, TSpace::dimension> TranslateMatrix;
            typedef AffineTransformation3D<TSpace> Transformation;
        public:
            /**
            * Constructor.
            * @param linear part of the transformation matrix :  should be 3x3
            * @param translate  part of the transformation : should be 3x1
            * @param origin origin of the transformation on the 3D object
            */
            AffineTransformation3D(const LinearMatrix &linear, const TranslateMatrix &translate,
                                   const TranslateMatrix &origin) :
                    _linear(linear), _inverse(linear.inverse()), _translate(translate), _origin(origin){}

            /**
             * Constructor.
             * @param trans transformation to use
             */
            AffineTransformation3D(const Transformation &trans) :
                    _linear(trans._linear), _inverse(trans._inverse), _translate(trans._translate),
                    _origin(trans._origin){}
            /**
            * Rotation operator
            * @param origin of the rotation
            * @param axis of rotation
            * @param angle of the rotation
            * @return Affine transformation corresponding to rotation
            */
            static Transformation rotation(const RealPoint &origin, const RealPoint &axis, const double &angle){
                double cos = std::cos(angle);
                double sin = std::sin(angle);
                double cosm = (1 - cos);
                RealPoint axNorm = axis.getNormalized();
                double ap = axNorm[0] * axNorm[1] * cosm + axNorm[2] * sin;
                double am = axNorm[0] * axNorm[1] * cosm - axNorm[2] * sin;
                double bp = axNorm[0] * axNorm[2] * cosm + axNorm[1] * sin;
                double bm = axNorm[0] * axNorm[2] * cosm - axNorm[1] * sin;
                double cp = axNorm[1] * axNorm[2] * cosm + axNorm[0] * sin;
                double cm = axNorm[1] * axNorm[2] * cosm - axNorm[0] * sin;
                LinearMatrix linear({(axNorm[0] * axNorm[0] * cosm + cos), am, bp,
                                     ap, (axNorm[1] * axNorm[1] * cosm + cos), cm,
                                     bm, cp, (axNorm[2] * axNorm[2] * cosm + cos)});

                for(int i = 0; i < 3; i++){
                    for(int j = 0; j < 3; j++){
                        linear(i, j) = std::round(linear(i, j) * 100) / 100;
                    }
                }

                return Transformation(linear, TranslateMatrix({0, 0, 0}),
                                      TranslateMatrix({origin[0], origin[1], origin[2]}));
            }
            /**
            * Backward operator, apply the transformation backward
            * @param input
            * @return Point with transformed coordinates
            */
            RealPoint backward(const Point &input) const{
                RealPoint p;
                double a = input[0] - _translate(0, 0) - _origin(0, 0);
                double b = input[1] - _translate(0, 1) - _origin(0, 1);
                double c = input[2] - _translate(0, 2) - _origin(0, 2);

                p[0] = _inverse(0, 0) * a + _inverse(0, 1) * b + _inverse(0, 2) * c + _origin(0, 0);
                p[1] = _inverse(1, 0) * a + _inverse(1, 1) * b + _inverse(1, 2) * c + _origin(0, 1);
                p[2] = _inverse(2, 0) * a + _inverse(2, 1) * b + _inverse(2, 2) * c + _origin(0, 2);
                return p;
            }
            /**
             * Backward operator, apply the transformation forward
             * @param input
             * @return Point with transformed coordinates
             */
            RealPoint forward(const Point &input) const{
                RealPoint p;

                p[0] = _linear(0, 0) * input[0] +
                       _linear(0, 1) * input[1] +
                       _linear(0, 2) * input[2] + _translate(0, 0);
                p[1] = _linear(1, 0) * input[0] +
                       _linear(1, 1) * input[1] +
                       _linear(1, 2) * input[2] + _translate(0,1);
                p[1] = _linear(2, 0) * input[0] +
                       _linear(2, 1) * input[1] +
                       _linear(2, 2) * input[2] + _translate(0,2);
                return p;
            }

        protected:
            LinearMatrix _linear;
            LinearMatrix _inverse;
            TranslateMatrix _translate;
            TranslateMatrix _origin;
        };


        /**
         * Description of template functor like class 'DomainRigidTransformation2D' <p>
         * \brief Aim: implements bounds of transformed domain.
         *
         * @tparam TDomain domain of the domain transfrer
         */
        template<typename TDomain>
        class DomainTransformer{
        public:
            typedef std::pair<typename TDomain::Space::Point, typename TDomain::Space::Point> Bounds;
            typedef typename TDomain::Space::Point Point;
            typedef typename TDomain::Space::RealPoint RealPoint;
        public:
            virtual Bounds operator()(const TDomain &input) const = 0;
        };


        /**
         * Description of template functor like class 'DomainTransformer2D' <p>
         * \brief Aim: implements bounds of transformed domain.
         *
         * @tparam TDomain a 2 dimensional domain.
         * @tparam TFunctor a functor which represent two dimensional rigid transformation.
         */
        template<typename TDomain, typename TFunctor>
        class DomainTransformer2D : public DomainTransformer<TDomain>{
            BOOST_CONCEPT_ASSERT((concepts::CDomain<TDomain>));
            BOOST_STATIC_ASSERT((TDomain::dimension == 2));
        public:
            typedef typename DomainTransformer<TDomain>::Point Point;
            typedef typename DomainTransformer<TDomain>::RealPoint RealPoint;
            typedef typename DomainTransformer<TDomain>::Bounds Bounds;
        public:
            explicit DomainTransformer2D(const TFunctor &functor)
                    : _functor(functor){}
            /**
            * Operator
            *
            * @return bounds of the transformed domain.
            */
            Bounds operator()(const TDomain &input) const{
                Point points[4];
                VectorRounding<RealPoint, Point> round;
                points[0] = round(_functor.forward(input.lowerBound()));
                points[1] = round(_functor.forward(input.upperBound()));
                points[2] = round(_functor.forward(Point(input.upperBound()[0], input.lowerBound()[1])));
                points[3] = round(_functor.forward(Point(input.lowerBound()[0], input.upperBound()[1])));

                Point t_min(INT_MAX, INT_MAX), t_max(INT_MIN, INT_MIN);
                for(auto &point: points){
                    if(point[0] < t_min[0]) t_min[0] = point[0];
                    if(point[1] < t_min[1]) t_min[1] = point[1];

                    if(point[0] > t_max[0]) t_max[0] = point[0];
                    if(point[1] > t_max[1]) t_max[1] = point[1];
                }

                Bounds bounds;
                bounds.first = t_min;
                bounds.second = t_max;
                return bounds;
            }

        protected:
            const TFunctor &_functor;
        };
        /**
         * Description of template functor like class 'DomainTransformer3D' <p>
         * \brief Aim: implements bounds of transformed domain.
         *
         * @tparam TDomain a 3 dimensional domain.
         * @tparam TFunctor a functor which represent two dimensional rigid transformation.
         */
        template<typename TDomain, typename TFunctor>
        class DomainTransformer3D : public DomainTransformer<TDomain>{
            BOOST_CONCEPT_ASSERT((concepts::CDomain<TDomain>));
            BOOST_STATIC_ASSERT((TDomain::dimension == 3));
        public:
            typedef typename DomainTransformer<TDomain>::Point Point;
            typedef typename DomainTransformer<TDomain>::RealPoint RealPoint;
            typedef typename DomainTransformer<TDomain>::Bounds Bounds;
        public:
            explicit DomainTransformer3D(const TFunctor &functor)
                    : _functor(functor){}

            Bounds operator()(const TDomain &input) const{
                Point points[8];
                VectorRounding<RealPoint, Point> round;
                points[0] = round(_functor(input.lowerBound()));
                points[1] = round(_functor(input.upperBound()));
                points[2] = round(_functor(Point(input.upperBound()[0], input.lowerBound()[1], input.lowerBound()[2])));
                points[3] = round(_functor(Point(input.lowerBound()[0], input.upperBound()[1], input.upperBound()[2])));
                points[4] = round(_functor(Point(input.upperBound()[0], input.lowerBound()[1], input.upperBound()[2])));
                points[5] = round(_functor(Point(input.lowerBound()[0], input.upperBound()[1], input.lowerBound()[2])));
                points[6] = round(_functor(Point(input.lowerBound()[0], input.lowerBound()[1], input.upperBound()[2])));
                points[7] = round(_functor(Point(input.upperBound()[0], input.upperBound()[1], input.lowerBound()[2])));

                Point t_min(INT_MAX, INT_MAX, INT_MAX), t_max(INT_MIN, INT_MIN, INT_MIN);
                for(Point &point: points){
                    if(point[0] < t_min[0]) t_min[0] = point[0];
                    if(point[1] < t_min[1]) t_min[1] = point[1];
                    if(point[2] < t_min[2]) t_min[2] = point[2];

                    if(point[0] > t_max[0]) t_max[0] = point[0];
                    if(point[1] > t_max[1]) t_max[1] = point[1];
                    if(point[2] > t_max[2]) t_max[2] = point[2];
                }

                Bounds bounds;
                bounds.first = t_min;
                bounds.second = t_max;
                return bounds;
            }

        protected:
            const TFunctor &_functor;
        };
    }
}
#endif

#undef AffineTransformation_RECURSES
#endif

