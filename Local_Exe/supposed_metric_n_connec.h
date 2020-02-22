#include <string>


class supposed_metric_n_connec {
    public:
    // from GUI input
    std::string Registre_file;  //file for communication register
    double convergeCrit;        // convergence criterion  
    double cfl,mach,AoA;
    int iterMax,Order,RK_step;                // Max iteration number, Scheme order and Runge-Kutta stages
    // from reader
    int ndime;
    int nelem;                  // real cells
    int npoints;
    int nhalo;
    int ncell;  //nelem+nhalo
    double** coord;
    int**   elem2node;
    int*    elem2vtk;
    int**   BoundIndex;
    // from connectivity
    int**   elem2elem;
    int**   node2elem;
    int**   face2elem;
    int**   elem2face;
    int*    face2facelocal;
    int**   node2node;
    int     nface;
    int**   face2point;
    int**   local2ZoneLocal;    //local element number TO other zone's local index
    // from Metric
    double**     face2norm;
    double***    cent2face;
    double*     face2area;
    double*     elem2vol;
    double**    elem2spectral;
    //double*     face2norm;
};