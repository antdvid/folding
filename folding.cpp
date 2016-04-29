#include "folding.h"

void Folder::computeRotateCenter(){

}

void Folder::inputFoldingSlices(Slice* s) {
    slices.push_back(s);
}

void Folder3d::doFlatten() {

}

void Folder3d::doFolding() {

}

Folder3d::Folder3d(SURFACE* s) {
    TRI* t;
    surf_tri_loop(s,t) {
	cells.push_back(new FoldTri(t));
    }
}
