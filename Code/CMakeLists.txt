PROJECT(EvalTopoTrans)

### Required in DGtal
CMAKE_MINIMUM_REQUIRED(VERSION 3.11)
SET(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

FIND_PACKAGE(DGtal REQUIRED)
FIND_PACKAGE(CGAL REQUIRED)
FIND_PACKAGE(GUDHI 2.0.0 REQUIRED)
FIND_PACKAGE(Boost COMPONENTS iostreams system filesystem REQUIRED)

INCLUDE(${CGAL_USE_FILE})

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${GUDHI_INCLUDE_DIRS} include)
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS} ${Boost_LIBRARY_DIR})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")

ADD_EXECUTABLE(evalTopoTrans2D
        src/evalTopoTrans2D.cpp
        include/AffineTransformation.h
        include/Interpolation.h
        include/ImageTransformation.h
        include/Topology.h
        include/Homology.h
        include/utils/CLI11.hpp
        include/utils/gnuplot-iostream.h)
TARGET_LINK_LIBRARIES(evalTopoTrans2D
        ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} ${CGAL_LIBRARY})

ADD_EXECUTABLE(evalTopoTrans3D
        src/evalTopoTrans3D.cpp
        include/AffineTransformation.h
        include/Interpolation.h
        include/ImageTransformation.h
        include/Topology.h
        include/Homology.h
        include/utils/CLI11.hpp
        include/utils/gnuplot-iostream.h)
TARGET_LINK_LIBRARIES(evalTopoTrans3D
        ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} ${CGAL_LIBRARY})

ADD_EXECUTABLE(viewer3D src/viewer3D.cpp)
TARGET_LINK_LIBRARIES(viewer3D ${DGTAL_LIBRARIES})
