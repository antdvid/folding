#include "drag.h"
#include <vector>

std::vector<Drag*> Drag::prototypes;

//this class need to be modified
//whenever adding a new type of drag
class DragProtoInit {
    
    PointDrag pd;
    GravityDrag gd;
    MultiplePointDrag md;
    SpinDrag sd;
    TuckDrag td;
    DragProtoInit() {
	Drag::prototypes.push_back(&pd);
	Drag::prototypes.push_back(&gd);
	Drag::prototypes.push_back(&md);
	Drag::prototypes.push_back(&sd);
	Drag::prototypes.push_back(&td);
    }
    static DragProtoInit singleton;
};

DragProtoInit DragProtoInit::singleton;
