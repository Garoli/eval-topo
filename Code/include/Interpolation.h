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
* @file Interpolation.h
* @author
*
*
* @date
*
* .
*/

#if defined(Interpolation_RECURSES)
#error Recursive header files inclusion detected in Interpolation.h
#else // defined(Interpolation_RECURSES)
/** Prevents recursive inclusion of headers. */
#define Interpolation_RECURSES

#if !defined Interpolation_h
/** Prevents repeated inclusion of headers. */
#define Interpolation_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include <cmath>
#include <algorithm>
#include <climits>
#include <utility>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/kernel/domains/CDomain.h>
#include <DGtal/kernel/CSpace.h>
//////////////////////////////////////////////////////////////////////////////

namespace DGtal::functors{
    template<typename TSpace, typename TImage>
    class Interpolation{
    public:
        typedef typename TSpace::RealPoint RealPoint;
        typedef typename TSpace::Point Point;
        typedef typename TImage::Value Value;

        explicit Interpolation(TImage *image, int width = 0) : _width(width), _image(image){
            Point adj;
            for(int i = -_width; i <= width; i++){
                adj[0] = i;
                for(int j = -_width; j <= width; j++){
                    adj[1] = j;
                    if(TSpace::dimension == 3){
                        for(int k = -_width; k <= width; k++){
                            adj[2] = k;
                            _adj.push_back(std::move(adj));
                        }
                    }else{
                        _adj.push_back(std::move(adj));
                    }
                }
            }
        }

        virtual double kernel(const double &t) const = 0;

        Value operator()(const RealPoint &pos) const{
            VectorRounding<RealPoint, Point> round;
            Point u = round(pos);
            Point v;
            double value = 0;
            double K;

            for(Point p : _adj){
                v = u + p;
                K = 1;
                for(int i = 0; i < TSpace::dimension; i++){
                    K *= kernel(pos[i] - v[i]);
                }
                if(_image->domain().isInside(v)){
                    value += _image->operator()(v) * K;
                }
            }

            value = std::min(value, (double) std::numeric_limits<Value>::max());

            return (Value) value;
        }

    protected:
        int _width;
        TImage *_image;
        std::vector<Point> _adj;
    };

    template<typename TSpace, typename TImage>
    class NearestNeighbour : public Interpolation<TSpace, TImage>{
    public:
        explicit NearestNeighbour(TImage *image, int width = 0) : Interpolation<TSpace, TImage>(image, width){}

        double kernel(const double &t) const{
            if((-0.5 <= t) && (t < 0.5)){
                return 1;
            }
            return 0;
        }

    };

    template<typename TSpace, typename TImage>
    class Bilinear : public Interpolation<TSpace, TImage>{
    public:
        explicit Bilinear(TImage *image, int width = 1) : Interpolation<TSpace, TImage>(image, width){}

        double kernel(const double &t) const{
            return std::max(0.0, 1 - std::abs(t));
        }

    };

    template<typename TSpace, typename TImage>
    class Bicubic : public Interpolation<TSpace, TImage>{
    public:
        explicit Bicubic(TImage *image, int width = 2, double alpha = -0.5) :
                Interpolation<TSpace, TImage>(image, width), _alpha(alpha){}

        double kernel(const double &t) const{
            double t_abs = std::abs(t);

            if(t_abs <= 1){
                return (_alpha + 2) * std::pow(t_abs, 3)
                       - (_alpha + 3) * std::pow(t_abs, 2)
                       + 1;
            }

            if((1 < t_abs) && (t_abs < 2)){
                return _alpha * std::pow(t_abs, 3)
                       - 5 * _alpha * std::pow(t_abs, 2)
                       + 8 * _alpha * t_abs
                       - 4 * _alpha;
            }

            return 0;
        }

        double _alpha;
    };

}// namespace DGtal

#endif // !defined Interpolation_h

#undef Interpolation_RECURSES
#endif // else defined(Interpolation_RECURSES)

