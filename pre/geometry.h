//
// Created by Reza on 6/10/22.
//

#ifndef SECTION_PROPERTIES_EXCUTABLE_GEOMETRY_H
#define SECTION_PROPERTIES_EXCUTABLE_GEOMETRY_H

#include "optional"
#include "vector"
#include "pre.h"


class geometry {
public:
    geometry(
            int geom,
            Material material,
            int tol = 12,
            std::vector<float> control_points = std::vector<float>(0)
    );
};


#endif //SECTION_PROPERTIES_EXCUTABLE_GEOMETRY_H
