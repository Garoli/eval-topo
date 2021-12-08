PROJECT(EvalTopoTrans)

### Required in DGtal
CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

FIND_PACKAGE(DGtal REQUIRED)

INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS} include)
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

ADD_EXECUTABLE(evalTopoTrans2D
        src/evalTopoTrans2D.cpp
        include/AffineTransformation.h
        include/Interpolation.h
        include/ImageTransformation.h)
TARGET_LINK_LIBRARIES(evalTopoTrans2D ${DGTAL_LIBRARIES})

ADD_EXECUTABLE(evalTopoTrans3D
        src/evalTopoTrans3D.cpp
        include/AffineTransformation.h
        include/Interpolation.h
        include/ImageTransformation.h)
TARGET_LINK_LIBRARIES(evalTopoTrans3D ${DGTAL_LIBRARIES})


ADD_EXECUTABLE(viewer3D src/viewer3D.cpp)
TARGET_LINK_LIBRARIES(viewer3D ${DGTAL_LIBRARIES})