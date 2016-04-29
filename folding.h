#ifndef FOLDING_H
#define FOLDING_H
#include<vector>
#include "fold_helper.h"
class Folder {
public:
    virtual void doFolding() = 0;
    virtual ~Folder(){}
protected:
    std::vector<Slice*> slices;
    std::vector<Cell*> cells;
    void computeRotateCenter();
    void inputFoldingSlices(Slice*);
    virtual void doFlatten() = 0;
    Folder(){}
};

class Folder3d:public Folder {
public:
    Folder3d(SURFACE*);
    ~Folder3d(){}
    void doFlatten();
    void doFolding();
    Folder3d(){}
};
#endif
