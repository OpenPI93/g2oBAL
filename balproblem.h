#ifndef BALPROBLEM_H
#define BALPROBLEM_H

#include "bal_type.h"
#include <vector>
#include <string>
#include <pangolin/pangolin.h>

using std::vector;
using std::string;

class balCamera
{
public:
    balCamera(double* data){
        Eigen::Vector3d so3;
        so3 << data[0], data[1], data[2];
        Eigen::Matrix3d R = Sophus::SO3d::exp(so3).matrix();
        Sophus::SE3d T(R, Eigen::Vector3d(data[3], data[4], data[5]));
        Eigen::Matrix<double, 6, 1> se3 = T.log();
        cam << se3(0, 0), se3(1, 0), se3(2, 0), se3(3, 0), se3(4, 0), se3(5, 0), data[6], data[7], data[8];
    }
    Eigen::Matrix<double, 9, 1> cam;
};

class balPoint
{
public:
    balPoint(double* data) : pt(data){}
    Eigen::Vector3d pt;
};

class balEdge
{
public:
    balEdge(int _cam, int _pt, double* data) : cam(_cam), pt(_pt), uv(data){}
    int cam;
    int pt;
    Eigen::Vector2d uv;
};

class balProblem
{
public:
    balProblem(const string& filename);
    void buildProblem();
    void getInfo();
    void solveProblem(int iter = 100);
    void writeToPLY(const string& filename);
    void showByPangolin();
protected:
    vector<balCamera> cameras;
    vector<balPoint> points;
    vector<balEdge> edges;
    g2o::SparseOptimizer optimizer;
};

#endif // BALPROBLEM_H
