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

/**
 * @file Topology.h
 * @author Sarah Brood (\c sarahbrood\@ecole.ensicaen.fr)
 * @author Heithem Dridi (\c heithemdridi\@ecole.ensicaen.fr)
 *
 * @date 2022/01
 *
 * This file is part of the project : Topological evaluation of image transformations
 */
#pragma once

#if defined(TopologyProject_RECURSES)
#error Recursive header files inclusion detected in TopologyProject.h
#else
#define TopologyProject_RECURSES

#if !defined TopologyProject_h
#define TopologyProject_h

#include <climits>
#include <cmath>
#include <iostream>
#include <iterator>
#include <utility>
#include <DGtal/base/Common.h>
#include <DGtal/geometry/curves/GridCurve.h>
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/images/imagesSetsUtils/SetFromImage.h>
#include <DGtal/io/boards/Board2D.h>
#include <DGtal/io/Color.h>
#include <DGtal/kernel/BasicPointFunctors.h>
#include <DGtal/kernel/CSpace.h>
#include <DGtal/kernel/domains/CDomain.h>
#include <DGtal/topology/helpers/Surfaces.h>
#include <DGtal/topology/SurfelAdjacency.h>
#include <DGtal/topology/CubicalComplex.h>
#include "DGtal/images/ConstImageAdapter.h"

namespace DGtal{

/**
* Description of template class 'ConnectedComponent' <p>
* \brief Aim: compute and returns the connected components of an image
* @code

typedef ImageTransformation<Space, Image, AffineTransformation, DomainTransformer, Interpolation> ImageTransformation;
ImageTransformation imageTransformation(transformation, domainTransformer, interpolation);
* @endcode
*
 * @tparam TImage
 * @tparam TDigitalSet Space digital set
 * @tparam TDigitalTopology Topology adjacency
 * @tparam TObject Object adjacency
 *
*/
    template<typename TImage, typename TDigitalSet, typename TDigitalTopology, typename TObject>
    struct ConnectedComponent{
        const TDigitalTopology &_dt;

        /**
         * Constructor.
         * @param dt adjacency
         */
        ConnectedComponent(TDigitalTopology dt) : _dt(dt){};

        /**
        * operator,
        * @param input adjacency returns the connected components of the image
        * @return Point with transformed coordinates
        */
        int operator()(const TImage &image, std::back_insert_iterator<std::vector<TObject>> inserter, int min, int max){
            TDigitalSet set(image.domain());
            SetFromImage<TDigitalSet>::template append<TImage>(set, image, min, max);
            TObject components(_dt, set);

            return components.writeComponents(inserter);
        }
    };


    /**
    * Description of template class 'Topology' <p>
    * \brief Aim: compute topological characteristic of an image
    *
    * @tparam TImage type of the image
    * @tparam TDigitalSet digital set of dimensions
    * @tparam TDigitalTopology adjacency
    * @tparam TObject type of the objects (connected components)
    */
    template<typename TImage, typename TDigitalSet, typename TDigitalTopology, typename TObject>
    class Topology{
    public:
        typedef typename TObject::ComplementObject ComplementObject;
        typedef typename TDigitalTopology::ReverseTopology BackgroundTopology;
        typedef ConnectedComponent<TImage, TDigitalSet, TDigitalTopology, TObject> ConnCompFG;
        typedef ConnectedComponent<TImage, TDigitalSet, BackgroundTopology, ComplementObject> ConnCompBG;

        /**
         * Constructor.
         * @param image the image, the charactestics have to be computed from
         * @param reverse reverse version of the image (see the function in evalTranspo examples)
         * @param dt digital topology to use (adjacency for foreground and background)
         */
        Topology(TImage *image, TImage *reverse, TDigitalTopology dt)
            : _image(image), _reverse(reverse), _ccFG(dt), _ccBG(dt.reverseTopology()){
        }

        virtual std::vector<int> operator()(int min, int max) = 0;

    protected:
        TImage *_image;
        TImage *_reverse;
        ConnCompFG _ccFG;
        ConnCompBG _ccBG;
        std::vector<TObject> _objFG;
        std::vector<ComplementObject> _objBG;

        /**
         * Return betti0 number
         * @param min min limit to compute betti number
         * @param max min limit to compute betti number
         * @return
         */
        int betti0(int min, int max){
            std::back_insert_iterator<std::vector<TObject>> inserter(_objFG);
            return _ccFG(*_image, inserter, min, max);
        }

        /**
         * Return betti number of the background
         * @param min min limit to compute betti number
         * @param max min limit to compute betti number
         * @return
         */
        int bettiBG(int min, int max){
            std::back_insert_iterator<std::vector<ComplementObject>> inserter(_objBG);
            return _ccBG(*_reverse, inserter, min, max);
        }
    };

    /**
    * Description of class 'Topology2D' <p>
    * \brief Aim: compute topological characteristic of a 2D image
    * @code

    typedef Topology2D<Image, KSpace, DigitalSet, DT4_8 , Object4_8> Topology2D;
    Topology2D inputTopo(&inputExtended, &inputReversed, dt4_8);

    * @endcode
    * @tparam TImage type of the image
    * @tparam TDigitalSet digital set of dimensions
    * @tparam TDigitalTopology adjacency
    * @tparam TObject type of the objects (connected components)
    */
    template<typename TImage, typename TKSpace, typename TDigitalSet, typename TDigitalTopology, typename TObject>
    class Topology2D : public Topology<TImage, TDigitalSet, TDigitalTopology, TObject>{
    public:

        /**
         * Constructor.
         * @param image the image, the charactesitics have to be computed from
         * @param reverse reverse version of the image (see the function in evalTranspo examples)
         * @param dt digital topology to use (adjacency for foreground and background)
         */
        Topology2D(TImage *image,  TImage *reverse, const TDigitalTopology &dt)
            : Topology<TImage, TDigitalSet, TDigitalTopology, TObject>(image, reverse, dt){}

        /**
         * Operator to returns topological characteristics
         * @param min
         * @param max
         * @return a vector with : {betti0, betti1}
         */
        std::vector<int> operator()(int min, int max){
            std::vector<int> bettis;
            bettis.push_back(Topology2D::betti0(min, max));
            bettis.push_back(Topology2D::bettiBG(min, max) - 1);
            return bettis;
        }

        /**
        * creates a board to display topological features
        * @param board to create in main (see evalTopoTrans for examples)
        */
        void toBoard(Board2D &board){
            typename TImage::Domain domain(Topology2D::_image->domain());
            TDigitalSet set(domain);
            GridCurve<TKSpace> c;
            unsigned int y;

            board << set << domain;

            for(unsigned int i = 0; i < Topology2D::_objFG.size(); i++){
                y = (i == 0) ? Topology2D::_objFG.size() - 1 : i - 1;
                sendToBoard(board, Topology2D::_objFG[i].pointSet(), Color::Green);
                if(Topology2D::_objFG.size() == 1){
                    c = getObjectBoundary(Topology2D::_objFG[i], Topology2D::_objFG[i].pointSet());
                }else{
                    c = getObjectBoundary(Topology2D::_objFG[i], Topology2D::_objFG[i].pointSet(), Topology2D::_objFG[y].pointSet());
                }
                board << CustomStyle(c.className(), new DGtal::CustomPenColor(Color::Red));
                board << c;
            }

            Color color[6] = {Color::Red, Color::Yellow, Color::Blue, Color::Cyan, Color::Magenta, Color::Purple};
            for(unsigned int i = 1; i < Topology2D::_objBG.size(); i++){
                y = (i == 0) ? Topology2D::_objBG.size() - 1 : i - 1;
                sendToBoard(board, Topology2D::_objBG[i].pointSet(), color[(i-1)%6]);
                if(Topology2D::_objBG.size() == 1){
                    c = getObjectBoundary(Topology2D::_objBG[i], Topology2D::_objBG[i].pointSet());
                }else{
                    c = getObjectBoundary(Topology2D::_objBG[i], Topology2D::_objBG[i].pointSet(), Topology2D::_objBG[y].pointSet());
                }
                board << CustomStyle(c.className(), new DGtal::CustomPenColor(Color::Red));
                board << c;
            }
        }

    private:
        void sendToBoard(Board2D &board, TDigitalSet &p_Object, DGtal::Color p_Color){
            board << CustomStyle(p_Object.className(), new DGtal::CustomFillColor(p_Color));
            board << p_Object;
        }

        template<typename Object>
        GridCurve<TKSpace> getObjectBoundary(Object &object, TDigitalSet set, TDigitalSet otherSet){
            typedef typename TKSpace::Point Point;
            typedef typename TKSpace::SCell SCell;
            GridCurve<TKSpace> boundaryCurve;
            SurfelAdjacency<TKSpace::dimension> sAdj(true);
            TKSpace kSpace;
            kSpace.init(set.domain().lowerBound(), set.domain().upperBound(), true);

            std::vector<SCell> boundaryCells;
            SCell aCell = Surfaces<TKSpace>::findABel(kSpace, set, *set.begin(), *otherSet.begin());
            Surfaces<TKSpace>::track2DBoundary(boundaryCells, kSpace, sAdj, set, aCell);

            for(SCell cell: boundaryCells){
                boundaryCurve.push_back(cell);
            }
            return boundaryCurve;
        }

        template<typename Object>
        GridCurve<TKSpace> getObjectBoundary(Object &object, TDigitalSet set){
            typedef typename TKSpace::Point Point;
            typedef typename TKSpace::SCell SCell;
            GridCurve<TKSpace> boundaryCurve;
            SurfelAdjacency<TKSpace::dimension> sAdj(true);
            TKSpace kSpace;
            kSpace.init(set.domain().lowerBound(), set.domain().upperBound(), true);

            std::vector<SCell> boundaryCells;
            SCell aCell = Surfaces<TKSpace>::findABel(kSpace, set, 1000);
            Surfaces<TKSpace>::track2DBoundary(boundaryCells, kSpace, sAdj, set, aCell);

            for(SCell cell: boundaryCells){
                boundaryCurve.push_back(cell);
            }
            return boundaryCurve;
        }
    };

    /**
    * Description of class 'Topology3D' <p>
    * \brief Aim: compute topological characteristic of a 3D object
    * @code

   typedef Topology3D<Image, KSpace, DigitalSet, DT6_26 , Object6_26> Topology3D;
   Topology3D inputTopo(&inputExtended, &inputReversed, dt6_26);

    * @endcode
    * @tparam TImage type of the image
    * @tparam TDigitalSet digital set of dimensions
    * @tparam TDigitalTopology adjacency
    * @tparam TObject type of the objects (connected components)
    */
    template<typename TImage, typename TKSpace, typename TDigitalSet, typename TDigitalTopology, typename TObject>
    class Topology3D : public Topology<TImage, TDigitalSet, TDigitalTopology, TObject>{
    public:

        /**
        * Constructor.
        * @param image the object, the charactesitics have to be computed from
        * @param reverse reverse version of the image (see the function in evalTranspo examples)
        * @param dt digital topology to use (adjacency for foreground and background)
        */
        Topology3D(TImage *image,  TImage *reverse, const TDigitalTopology &dt)
                : Topology<TImage, TDigitalSet, TDigitalTopology, TObject>(image, reverse, dt) {}

        /**
         * Operator to returns topological characteristics
         * @param min
         * @param max
         * @return a vector with : {betti0, betti1, betti2}
         */
        std::vector<int> operator()(int min, int max){
            std::vector<int> bettis;
            int b0 = Topology3D::betti0(min, max);
            int b2 = Topology3D::bettiBG(min, max);
            bettis.push_back(b0);
            bettis.push_back(betti1(b0, b2));
            bettis.push_back(b2);
            return bettis;
        }

    protected:
        /**
         * Compute betti1 from betti0, betti2 and euler characteristics : b1=e+b0+b2
         * @param b0 betti number 0
         * @param b2 betti number 1
         * @return betti number 2
         */
        int betti1(int b0, int b2) const{
            typedef std::map<typename TKSpace::Cell, CubicalCellData>   Map;
            typedef CubicalComplex<TKSpace, Map> CubCmplx;

            Z3i::KSpace ks;
            bool space_ok = ks.init(Topology3D::_image->domain().lowerBound(),
                                    Topology3D::_image->domain().upperBound(), true);
            if (!space_ok)
            {
                trace.error() << "Error in the Khamisky space construction."<<std::endl;
                return 2;
            }

            CubCmplx complex(ks);

            TDigitalSet set(Topology3D::_image->domain());
            complex.construct(set);
            return complex.euler() + b0 + b2;
        }
    };
}// namespace DGtal

#endif // !defined TopologyProject_h

#undef TopologyProject_RECURSES
#endif // else defined(TopologyProject_RECURSES)

