#ifndef VTU_H
#define VTU_H

#include "euler_eq.h"
#include "shape_f.h"
#include "type_def.h"

void writeFRSolutionVTU(const std::string& filename,
                        const std::vector<Q9>& nodes,
                        const Mesh& mesh,
                        double gamma = GAMMA);

#endif //VTU_H
