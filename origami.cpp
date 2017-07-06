#include "origami.h"
#include <iostream>
#include <iterator>
#include <unordered_map>
#include <fstream>

static std::ostream& operator << (std::ostream& os, const std::vector<double>&);
static std::ostream& operator << (std::ostream& os, const std_matrix&);

//class OgmPoint
OgmPoint::OgmPoint(SpringVertex& sv) : SpringVertex(sv)
{
    x0.resize(3,0);
    std::copy(sv.getCoords(), sv.getCoords()+3, x0.begin());
}

//class OrigamiFold
void OrigamiFold::preprocess(std::vector<SpringVertex*>& pts)
{
     if (first)
     { 
         first = false;
	 for (size_t i = 0; i < pts.size(); ++i)
         {
              SpringVertex* sv = new OgmPoint(*pts[i]);
              delete pts[i];
              pts[i] = sv;
              pts[i]->setRegistered();
         }
	 assignFacesOgmPoints(pts); 
     }
     findNextFoldingAngle(); 
     for (int i = 0; i < creases_.size(); i++) 
	  creases_[i]->updateRotMatrix(rho_delta[i]); 
     for (auto it : faces_) 
          it->updateFoldingMatrix(); 
}

void OrigamiFold::assignFacesOgmPoints(std::vector<SpringVertex*>& pts) {
    for (size_t i = 0; i < pts.size(); i++) {
         OgmPoint* op = static_cast<OgmPoint*>(pts[i]); 
         
         for (auto it : faces_) 
              if (it->poInside(op)) {
                  op->addOgmFace(it); 
              }
    }     
}

std::vector<double> OrigamiFold::findFoldable(std::vector<double>& input, 
double& err)
{
    if (m_opt == NULL)
    {
        size_t N = rho.size();
        m_opt = new nlopt::opt(nlopt::LN_COBYLA, N);
        this->m_opt->set_stopval(1e-6);
        this->m_opt->set_ftol_abs(1e-6);
        this->m_opt->set_xtol_rel(1e-6);
        this->m_opt->set_maxeval(3000);
        nlopt::srand(1234);
        std::vector<double> lowerbounds(N, -M_PI);
        std::vector<double> upperbounds(N, M_PI);

        for (size_t i = 0; i < N; ++i)
        {
            if (rho_T[i] > 0)
                lowerbounds[i] = 0;
            if (rho_T[i] < 0)
                upperbounds[i] = 0;
        }
        m_opt->set_lower_bounds(lowerbounds);
        m_opt->set_upper_bounds(upperbounds);
        this->m_opt->set_min_objective(staticTargetFunction, this);
    }
    std::vector<double> ans = m_opt->optimize(input);
    err = targetFunction(ans);
    return ans;
}

double OrigamiFold::staticTargetFunction(
                const std::vector<double>& x,
                std::vector<double> &grad,
                void *my_func_data)
{
    OrigamiFold* obj = static_cast<OrigamiFold*>(my_func_data);
    return obj->targetFunction(x);
}

double OrigamiFold::targetFunction(
                const std::vector<double>& x)
{
    double ans = 0;
    
    for (size_t i = 0; i < creases_.size(); i++) {
         creases_[i]->updateRotMatrix(x[i]); 
    }
    for (auto it : vertices_) {
         if (!it->judgeNonVirtual()) continue; 
         Nvvertex* nv = static_cast<Nvvertex*>(it); 

         std::vector<Crease*> crsOnVertex = nv->getCreaseOnVertex();    
         std_matrix I = Math::Eye(4); 
         std_matrix M = Math::Eye(4); 
         std_matrix tempAns = Math::Mat(4,4); 

         for (auto it1 : crsOnVertex) {

              Math::MatxMat(M, it1->getRotMatrix(), tempAns);
              Math::assignMatToMat(tempAns, M);        
         } 
         Math::MatmMat(M, I, M); 
         ans += Math::Norm(M, 1); 
    }
    return ans;
}

void OrigamiFold::findNextFoldingAngle()
{
    const int maxIter = 2000;
    static int iter = 0;
    bool success = false;
    static double w = 0.8; 
    double w0, w1, w2, D;
    w1 = 0.2; w2 = 0.01; D = 0.015;
    while (!success && iter++ < maxIter)
    {
        int N = rho.size();
        std::vector<double> rho_rand(N, 0);
        std::vector<double> dir(N, 0);
        for (int i = 0; i < N; ++i)
        {
            rho_rand[i] = (rand()/(double)RAND_MAX * 2 - 1) * M_PI;
            dir[i] = (1 - w) * rho_rand[i] + w * rho_T[i];
            rho_tau[i] = rho_delta[i] + D * dir[i];
            rho_tau[i] = std::min(M_PI,
                         std::max(rho_tau[i], -M_PI));
        }

        double err = 0;
        rho = findFoldable(rho_tau, err);
        if (isValid(rho) && Math::dist3d(rho_T, rho) < Math::dist3d(rho_T, rho_delta))
        {
            rho_delta = rho;
            w = w + w1;
            success = true;
            std::cout << "rho = " << rho << std::endl;
            if (Math::dist3d(rho_delta, rho_T) < 1.0e-3) {
                iter = maxIter; 
                std::cout << "convergence reached: err = 1e-3" 
                << "folding angle = [" << rho_delta << "]" << std::endl; 
                std::cout << "Folding process is terminated" << std::endl; 
                m_t = 0; 
            }
        }
        else
        {
            w = w - w2;
            success = false;
        }
        w = std::max(std::min(w, 1.0), 0.0);
    }
    if (iter >= maxIter)
    {
        std::cout << "Maxium iteration reached\n"
        << "folding angle = [" << rho_delta << "]"
        << std::endl;
        std::cout << "err = " << targetFunction(rho_delta) << std::endl;
        std::cout << "Folding process is terminated"
        << std::endl;
        m_t = 0;
    }
}

bool OrigamiFold::isValid(const std::vector<double>& rho)
{
    for (auto a : rho)
    {
        if (fabs(a) > M_PI)
            return false;
    }
    return true;
}

void OrigamiFold::postprocess(std::vector<SpringVertex*>& pts) {}

void OrigamiFold:: setVel(SpringVertex* sv)
{
    std::vector<double> new_crds(3, 0);
    ogmComputeNewPosition(sv, new_crds);
    double* vel = sv->getVel();
    for (int i = 0; i < 3; ++i)
        vel[i] = (new_crds[i] - sv->getCoords()[i])/getTimeStepSize();

}

void OrigamiFold::setAccel(SpringVertex*){}
size_t OrigamiFold::dataSize()
{
    return totalDataSize;
}

Drag* OrigamiFold::clone(const Info & info)
{
    if (info.data().size() == 0)
    {
        std::cout << "Warning: OrigamiFold is not created, no data provided"
                  << std::endl;
        return NULL;
    }
    if (info.data().size() > 0 && info.data().size() < 2)
    {
        std::cout << "Warning: OrigamiFold is not created, "
                  << "insufficient data is given"
                  << std::endl;
    }
    std::vector<double> v = info.data();
    std::vector<std::vector<double>> points;
    std::vector<std::pair<int, int>> creases;
    std::vector<std::vector<int>> faces; 
    std::vector<std::vector<int>> vertex_crease_index; 
    std::vector<double> angles;
    size_t n_pt = (size_t)v[0];
    size_t n_crs = (size_t)v[1];
    size_t n_fs = (size_t)v[2]; 
    size_t n_fsp = (size_t)v[3]; 
    size_t n_nvv = (size_t)v[4];
    size_t n_nvc = (size_t)v[5]; 

    totalDataSize = 6 + n_pt * 3 + n_crs * 2 + n_crs * 1 + n_fs + n_fsp + 
        n_nvv + n_nvc;
    if(!validateData(info))
        return NULL;

    double *it = &(v.front()) + 6;
    for (int i = 0; i < n_pt; ++i)
    {
        points.push_back({*it, *(it+1), *(it+2)});
        it += 3;
    }
    for (int i = 0; i < n_crs; ++i)
    {
        creases.push_back(std::make_pair(int(*it), int(*(it+1))));
        it += 2;
    }
    for (int i = 0; i < n_fs; i++) {
	 std::vector<int> temp; 
         for (int j = 0; j < int(*it); j++) {
	      temp.push_back(int(*(it+j+1))); 
         }
	 it = it + int(*it)+1; 
	 faces.push_back(temp);
    }
    for (int i = 0; i < n_nvv; i++) {
         std::vector<int> temp; 
         for (int j = 0; j < int(*it); j++)
              temp.push_back(int(*(it+j+1))); 
         it = it + int(*it) + 1; 
         vertex_crease_index.push_back(temp); 
    }
    for (int i = 0; i < n_crs; ++i)
    {
        angles.push_back(*it);
        it ++;
    }
    return new OrigamiFold(points, faces, creases, vertex_crease_index, angles);
}

OrigamiFold::OrigamiFold(): m_opt(NULL) {}

OrigamiFold::OrigamiFold(const std::vector<std::vector<double>>& points,
			 const std::vector<std::vector<int>>& fs, 
                         const std::vector<std::pair<int, int>>& creases,
                         const std::vector<std::vector<int>>& 
                            vertex_crease_index, 
                         const std::vector<double>& angles)
{
    m_t = 1.4;
    m_opt = NULL;
    if (creases.size() != angles.size())
    {
        std::cerr << "Number of creases does not equal to number of angles"
                  << std::endl;
    }

    size_t N = creases.size();

    rho.resize(N, 0);
    rho_delta.resize(N, 0);
    rho_tau.resize(N, 0);
    vertices_.resize(points.size()); 
    for (int i = 0; i < points.size(); i++) 
	 vertices_[i] = NULL;  
    for (size_t i = 0; i < creases.size(); ++i)
    {
        int vindex = creases[i].first;
        int pindex = creases[i].second;

        if (vertices_[vindex] == NULL)
            vertices_[vindex] = new Nvvertex(points[vindex], vindex);
	if (vertices_[pindex] == NULL) 
	    vertices_[pindex] = new Vertex(points[pindex], pindex); 

        creases_.push_back(new Crease(vertices_[vindex], vertices_[pindex], 
            vindex, pindex, angles[i]));
    }
    for (size_t i = 0; i < points.size(); i++) {
         if (vertices_[i] == NULL) 
             vertices_[i] = new Vertex(points[i], i); 
    }

    for (auto it : creases_) 
         rho_T.push_back(it->getRotAngle()); 

    std::vector<Crease*> tempFaceCrease; 
    
    for (size_t i = 0; i < fs.size(); i++) {
	 std::vector<Vertex*> tempFaceVertex; 
         
	 tempFaceCrease.push_back(creases_[i]); 
	 for (size_t j = 0; j < fs[i].size(); j++) 
	      tempFaceVertex.push_back(vertices_[fs[i][j]]); 
	 faces_.push_back(new Face(tempFaceVertex, tempFaceCrease)); 
    }

    int count = 0; 

    for (auto it : vertices_) {
         if (!it->judgeNonVirtual()) continue; 
         Nvvertex* nv = static_cast<Nvvertex*>(it); 
         for (size_t i = 0; i < vertex_crease_index[count].size(); i++) 
              nv->insertCrease(creases_[vertex_crease_index[count][i]]); 
         count++; 
    }
}

void OrigamiFold::ogmComputeNewPosition(SpringVertex* sv, std::vector<double>& new_crds)
{
    OgmPoint* op = static_cast<OgmPoint*>(sv); 

    new_crds[0] = op->getInitialCoords()[0]; 
    new_crds[1] = op->getInitialCoords()[1];
    new_crds[2] = op->getInitialCoords()[2];
    // note: for the point on the crease, 
    // we only move it once (rotate along creases of previous face or
    // current face. The two motion is equivalent since that point is on
    // the crease line)
    Face* bface = op->getOgmFaces().back(); 

    bface->crsFoldCrds(new_crds); 
}

OrigamiFold::~OrigamiFold() {
    if (this->m_opt)
    {
        delete this->m_opt;
        this->m_opt = NULL;
    }
    for (auto v : vertices_)
         delete v; 
    for (auto c : creases_) 
         delete c; 
    for (auto f : faces_) 
         delete f; 
}

// class Crease
Crease::Crease(Vertex* v1, Vertex* v2, int idx1, int idx2, double rho) {
    v1_ = v1; 
    v2_ = v2; 
    index1 = idx1; 
    index2 = idx2; 
    rho_T = rho; 
    rot_matrix = Math::Mat(4, 4); 
    dir.resize(3); 
    Math::VecmVec(v2_->getCoords(), v1_->getCoords(), dir); 
    Math::Normalize(dir); 
}

void Crease::updateRotMatrix(double rho) {
    Math::getRotMatrix(v1_->getCoords(), getDir(), rot_matrix, rho); 
}

// class Face
Face::Face(const std::vector<Vertex*>& vv, const std::vector<Crease*>& vc) {
    vertices_ = vv; 
    crsAlongPath = vc; 
    fd_matrix = Math::Mat(4, 4);
}

bool Face::poInside(OgmPoint* op) {
    int size = vertices_.size(); 

    for (size_t i = 0; i < size-1; i++)
         if (Math::leftOn(vertices_[i]->getCoords(), 
            vertices_[i+1]->getCoords(), op->getInitialCoords())) continue; 
         else return false; 
    if (!Math::leftOn(vertices_[size-1]->getCoords(), 
            vertices_[0]->getCoords(), op->getInitialCoords()))
        return false; 
    return true; 
}

void Face::crsFoldCrds(std::vector<double>& new_crds) {
    std::vector<double> coords(new_crds); 
    std::vector<double> ans(4, 0); 

    coords.push_back(1.0);
    Math::MatxVec(fd_matrix, coords, ans);
    for (int i = 0 ;i < 3; i++) 
         new_crds[i] = ans[i]; 
}

void Face::updateFoldingMatrix() {
    std_matrix M = Math::Eye(4); 
    std_matrix tempAns = Math::Mat(4, 4); 

    for (auto it : crsAlongPath) {
         Math::MatxMat(M, it->getRotMatrix(), tempAns); 
         Math::assignMatToMat(tempAns, M); 
    }
    Math::assignMatToMat(M, fd_matrix); 
}

// class Math
void Math::MatxMat(const std::vector<std::vector<double>> & M1,
                   const std::vector<std::vector<double>> & M2,
                   std::vector<std::vector<double>> & ans)
{
    // mxn x nxp
    int m = M1.size(); 
    int n = M1[0].size(); 
    int p = M2[0].size(); 

    for (int i = 0; i < m; ++i)
    for (int j = 0; j < p; ++j)
    {
        ans[i][j] = 0.0;
        for (int k = 0; k < n; ++k)
            ans[i][j] += M1[i][k] * M2[k][j];
    }
}

void Math::MatxVec(const std::vector<std::vector<double>>& M,
                   const std::vector<double>& v,
                   std::vector<double>& ans)
{
    int m = M.size(); 
    int n = v.size(); 

    for (int i = 0; i < m; ++i)
    {
        ans[i] = 0;
        for (int j = 0; j < n; ++j)
            ans[i] += M[i][j] * v[j];
    }
}

void Math::MatpMat(const std::vector<std::vector<double>>& M1,
                 const std::vector<std::vector<double>>& M2,
                 std::vector<std::vector<double>> & ans)
{
    int m = M1.size(); 
    int n = M1[0].size(); 

    for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
    {
        ans[i][j] = M1[i][j] + M2[i][j];
    }
}

void Math::MatmMat(const std::vector<std::vector<double>>& M1,
                 const std::vector<std::vector<double>>& M2,
                 std::vector<std::vector<double>> & ans)
{
    int m = M1.size();
    int n = M1[0].size();

    for (int i = 0; i < m; ++i)
    for (int j = 0; j < n; ++j)
    {
        ans[i][j] = M1[i][j] - M2[i][j];
    }
}

void Math::VecmVec(const std::vector<double>& v1,
                 const std::vector<double>& v2,
                       std::vector<double>& ans)
{
    for (size_t i = 0; i < v1.size(); ++i)
        ans[i] = v1[i] - v2[i];
}

void Math::VecpVec(const std::vector<double>& v1,
                 const std::vector<double>& v2,
                       std::vector<double>& ans)
{
    for (size_t i = 0; i < v1.size(); ++i)
        ans[i] = v1[i] + v2[i];
}

void Math::Cross3d(const std::vector<double>& u,
                     const std::vector<double>& v,
                     std::vector<double>& ans)
{
    ans[0] = u[1]*v[2] - u[2]*v[1];
    ans[1] = u[2]*v[0] - u[0]*v[2];
    ans[2] = u[0]*v[1] - u[1]*v[0];
}

double Math::TriProd(const std::vector<double>& a,
                     const std::vector<double>& b,
                     const std::vector<double>& c)
{
    std::vector<double> crx(3, 0);
    Math::Cross3d(b,c,crx);
    return Math::Dot3d(a, crx);
}

double Math::Dot3d(const std::vector<double>& a,
                   const std::vector<double>& b)
{
    double ans = 0;
    for (size_t i = 0; i < a.size(); ++i)
    {
        ans += a[i] * b[i];
    }
    return ans;
}

double Math::Mag(const std::vector<double>& v)
{
    double ans = 0;
    for (auto vi : v)
        ans += vi*vi;
    return sqrt(ans);
}

std_matrix Math::Eye(int N)
{
    std_matrix res(N, std::vector<double>(N, 0));
    for (int i = 0; i < N; ++i)
        res[i][i] = 1.0;
    return res;
}

std_matrix Math::Mat(int m, int n)
{
   return std_matrix(m, std::vector<double>(n, 0));
}


double Math::Norm(const std_matrix& m, int p)
{
   if (p == 0)
        return 1.0;
   double res = 0;
   for (auto v : m)
   for (auto a : v)
   {
        res += std::pow(fabs(a), p);
   }
   return std::pow(res, 1/p);
}

double Math::Norm(const std::vector<double>& v, int p)
{
    if (p == 0)
        return 1.0;
    double res = 0;
    for (auto a : v)
        res += std::pow(fabs(a), p);
    return std::pow(res, 1/p);
}


double Math::dist3d(const std::vector<double>& p1,
                    const std::vector<double>& p2)
{
    double ans = 0;
    for (size_t i = 0; i < p1.size(); ++i)
    {
        ans += (p1[i]-p2[i])*(p1[i]-p2[i]);
    }
    return sqrt(ans);
}

void Math::Normalize(std::vector<double>& v)
{
   double vm = Mag(v);
   if (vm == 0)
        return;
   for (size_t i = 0; i < v.size(); ++i)
        v[i] /= vm;
}

bool Math::leftOn(const std::vector<double>& p1, const std::vector<double>& p2,
     const std::vector<double>& q) {
    std::vector<double> a(3, 0);
 
    Math::VecmVec(p2, p1,a);

    std::vector<double> b(3, 0);
    
    Math::VecmVec(q, p1, b);
 
    std::vector<double> ans(3, 0); 

    Math::Cross3d(a, b, ans); 
    if (fabs(ans[2]) < 1.0e-12 || ans[2] > 0) 
        return true; 
    else return false; 
}

void Math::getRotMatrix(const std::vector<double>& point, 
    const std::vector<double>& dir, std_matrix& mat, double rho) {
    double u = dir[0], v = dir[1], w = dir[2]; 
    double a = point[0], b = point[1], c = point[2]; 

    mat[0][0] = u*u+(1-u*u)*cos(rho); 
    mat[0][1] = u*v*(1-cos(rho))-w*sin(rho); 
    mat[0][2] = u*w*(1-cos(rho))+v*sin(rho); 
    mat[0][3] = (a*(v*v+w*w)-u*(b*v+c*w))*(1-cos(rho))+(b*w-c*v)*sin(rho); 

    mat[1][0] = u*v*(1-cos(rho))+w*sin(rho); 
    mat[1][1] = v*v+(1-v*v)*cos(rho); 
    mat[1][2] = v*w*(1-cos(rho))-u*sin(rho); 
    mat[1][3] = (b*(u*u+w*w)-v*(a*u+c*w))*(1-cos(rho))+(c*u-a*w)*sin(rho); 

    mat[2][0] = u*w*(1-cos(rho)) - v*sin(rho);
    mat[2][1] = v*w*(1-cos(rho)) + u*sin(rho);
    mat[2][2] = w*w+(1-w*w)*cos(rho);
    mat[2][3] = (c*(u*u+v*v)-w*(a*u+b*v))*(1-cos(rho))+(a*v-b*u)*sin(rho); 
    
    mat[3][0] = mat[3][1] = mat[3][2] = 0.0; 
    mat[3][3] = 1.0; 
}

void Math::assignMatToMat(const std_matrix& M1, std_matrix& M2) {
    int m = M1.size(); 
    int n = M1[0].size(); 

    for (int i = 0; i < m; i++)
    for (int j = 0; j < n; j++) 
         M2[i][j] = M1[i][j]; 
} 

std::ostream& operator << (std::ostream& os, const std::vector<double>& v)
{
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(os, " "));
    return os;
}

std::ostream& operator << (std::ostream& os, const std_matrix& m)
{
    os << "[" << std::endl;
    for (auto v : m)
    {
        std::copy(v.begin(), v.end(), std::ostream_iterator<double>(os, " "));
        os << std::endl;
    }
    os << "]" << std::endl;
    return os;
}
