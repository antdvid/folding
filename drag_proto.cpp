#include "drag.h"
#include <vector>

std::vector<Drag*> Drag::prototypes;

//this class need to be modified
//whenever adding a new type of drag
class DragProtoInit {
    
    PointDrag pd;
    GravityDrag gd;
    MultiplePointDrag md;
    FoldDrag sd;
    CloseUmbrellaDrag cud;
    RelaxDrag rd;
    GravityBoxDrag gbd;
    CompressDrag cpd;
    DragProtoInit() {
	Drag::prototypes.push_back(&pd);
	Drag::prototypes.push_back(&gd);
	Drag::prototypes.push_back(&md);
	Drag::prototypes.push_back(&sd);
	Drag::prototypes.push_back(&cud);
	Drag::prototypes.push_back(&rd);
	Drag::prototypes.push_back(&gbd);
	Drag::prototypes.push_back(&cpd);
    }
    static DragProtoInit singleton;
};

DragProtoInit DragProtoInit::singleton;
