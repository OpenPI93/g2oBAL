#ifndef BAL_TYPE_H
#define BAL_TYPE_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sophus/se3.hpp>
#include <Eigen/Core>

#include "g2o/stuff/sampler.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/robust_kernel_impl.h"
#include "g2o/core/batch_stats.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/optimization_algorithm_dogleg.h"
#include <g2o/solvers/csparse/linear_solver_csparse.h>

#include "g2o/solvers/cholmod/linear_solver_cholmod.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/solvers/eigen/linear_solver_eigen.h"
#include "g2o/solvers/pcg/linear_solver_pcg.h"
#include "g2o/types/sba/types_six_dof_expmap.h"

#include "g2o/solvers/structure_only/structure_only_solver.h"


/////////////////////////////// Vertex ///////////////////////////////////////////
/**
 * @brief The bal3DVertex class
 * base 3D point
*/
class bal3DVertex : public g2o::BaseVertex<3, Eigen::Vector3d >
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    bal3DVertex(){}
    bool read(std::istream &is){}
    bool write(std::ostream &os) const{}
    virtual void setToOriginImpl()
    {
        _estimate.fill(0);
    }
    virtual void oplusImpl(const double* _update){
        Eigen::Map< const Eigen::Vector3d > update(_update);
        _estimate += update;
    }
};

/**
 * @brief The balCameraVertex class
 * Eigen::Matrix<double, 9, 1>(R, t, f, k1 ,k2)
 */
class balCameraVertex : public g2o::BaseVertex<9, Eigen::Matrix<double, 9, 1> >
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    balCameraVertex(){}
    bool read(std::istream &is){}
    bool write(std::ostream &os) const{}
    virtual void setToOriginImpl()
    {
        _estimate.fill(0);
    }
    virtual void oplusImpl(const double* _update)
    {
        Eigen::Map<const Eigen::Matrix<double, 9, 1> > update(_update);
        Eigen::Matrix<double, 6, 1> se3;
        se3 << update(0, 0), update(1, 0), update(2, 0), update(3, 0), update(4, 0), update(5, 0);
        Sophus::SE3d T = Sophus::SE3d::exp(se3);
        double f = update(6, 0) + _estimate(6, 0), k1 = update(7, 0) + _estimate(7, 0), k2 = update(8, 0) + _estimate(8, 0);

        Eigen::Matrix<double, 6, 1> pose;
        pose << _estimate(0, 0), _estimate(1, 0), _estimate(2, 0), _estimate(3, 0), _estimate(4, 0), _estimate(5, 0);

        Sophus::SE3d SE3_pose = Sophus::SE3d::exp(pose);
        SE3_pose = T * SE3_pose;

        pose = SE3_pose.log();

        _estimate << pose(0, 0), pose(1, 0), pose(2, 0), pose(3, 0), pose(4, 0), pose(5, 0), f, k1, k2;
    }

};

/**
 * @brief The baseCameraVertex class
 * Eigen::Matrix<double, 6, 1>(R, t)
 */
class baseCameraVertex : public g2o::BaseVertex<6, Eigen::Matrix<double, 6, 1> >
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    baseCameraVertex(){}
    bool read(std::istream &is){}
    bool write(std::ostream &os) const{}
    virtual void setToOriginImpl(){
        _estimate.fill(0);
    }
    virtual void oplusImpl(const double* _update)
    {
        Eigen::Map<const Eigen::Matrix<double, 6, 1> > update(_update);
        Sophus::SE3d T = Sophus::SE3d::exp(update);
        Sophus::SE3d pose = Sophus::SE3d::exp(_estimate);
        pose = T * pose;
        _estimate = pose.log();

    }

};

////////////////////////////////////end of vertix///////////////////////////////////////////////




///////////////////////////////////////  Edge   ///////////////////////////////////////////////
class balCamPtEdge : public g2o::BaseBinaryEdge<2, Eigen::Vector2d, balCameraVertex, bal3DVertex>
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    balCamPtEdge(){}
    bool read(std::istream &is){}
    bool write(std::ostream &os) const{}
    virtual void computeError()
    {
        const balCameraVertex* v0 = static_cast<const balCameraVertex*>(_vertices[0]);
        const bal3DVertex* v1 = static_cast<const bal3DVertex*>(_vertices[1]);

        Eigen::Matrix<double, 9, 1>  cam = v0->estimate();
        Eigen::Matrix<double, 6, 1> se3;
        se3 << cam(0, 0), cam(1, 0), cam(2, 0), cam(3, 0), cam(4, 0), cam(5, 0);
        double f = cam(6, 0);
        double k1 = cam(7, 0);
        double k2 = cam(8, 0);

        Sophus::SE3d T = Sophus::SE3d::exp(se3);
        Eigen::Vector3d pt(v1->estimate());

        Eigen::Vector3d cur_pt = T * pt;
        cur_pt = -cur_pt / cur_pt(2, 0);

        double r_2 = cur_pt(0, 0) * cur_pt(0, 0) + cur_pt(1, 0) * cur_pt(1, 0);
        double rp = 1.0 + k1 * r_2 + k2 * r_2 * r_2;

        Eigen::Vector2d cur_uv;
        cur_uv << cur_pt(0, 0) * f * rp, cur_pt(1, 0) * f * rp;

        Eigen::Vector2d obs(_measurement);

        _error = obs - cur_uv;
    }

    virtual void linearizeOplus()override
    {
        const balCameraVertex* v0 = static_cast<const balCameraVertex*>(_vertices[0]);
        const bal3DVertex* v1 = static_cast<const bal3DVertex*>(_vertices[1]);

        Eigen::Matrix<double, 9, 1> cam(v0->estimate());
        Eigen::Matrix<double, 6, 1> se3;
        se3 << cam(0, 0), cam(1, 0), cam(2, 0), cam(3, 0), cam(4, 0), cam(5, 0);
        Sophus::SE3d T = Sophus::SE3d::exp(se3);

        Eigen::Vector3d pt = v1->estimate();
        Eigen::Vector3d cur_pt = T * pt;
        double x = cur_pt(0, 0);
        double y = cur_pt(1, 0);
        double z = cur_pt(2, 0);
        double x_2 = x * x;
        double x_4 = x_2 * x_2;
        double y_2 = y * y;
        double y_4 = y_2 * y_2;
        double invz = 1.0 / z;
        double invz_2 = invz * invz;
        double invz_4 = invz_2 * invz_2;
        double f = cam(6, 0);
        double k1 = cam(7, 0);
        double k2 = cam(8, 0);

        cur_pt = -cur_pt / cur_pt(2, 0);

        double r_2 = (x * x + y * y) * invz_2;
        double rp = 1 + k1 * r_2 + k2 * r_2 * r_2;

        Eigen::Matrix<double, 2, 3, Eigen::ColMajor> tmp; //this is dh/dp which h means observation equation

        tmp(0, 0) = -f * invz * (1 + k1 * (3 * x_2 + y_2) * invz_2 + k2 * (5 * x_4 + 6 * x_2 * y_2 + y_4) * invz_4);
        tmp(0, 1) = -f * 2 * x * y * invz * invz_2 * (k1 + 2 * k2 * r_2);
        tmp(0, 2) = f * x * invz_2 * (1 + 3 * k1 * r_2 + 5 * k2 * r_2 * r_2);

        tmp(1, 0) = -f * 2 * x * y * invz * invz_2 * (k1 + 2 * k2 * r_2);
        tmp(1, 1) = -f * invz * (1 + k1 * (x_2 + 3 * y_2) * invz_2 + k2 * (5 * y_4 + 6 * x_2 * y_2 + x_4) * invz_4);
        tmp(1, 2) = f * y * invz_2 * (1 + 3 * k1 * r_2 + 5 * k2 * r_2 * r_2);

        _jacobianOplusXj = -tmp * T.rotationMatrix();

        _jacobianOplusXi(0, 0) = tmp(0, 0);
        _jacobianOplusXi(0, 1) = tmp(0, 1);
        _jacobianOplusXi(0, 2) = tmp(0, 2);
        _jacobianOplusXi(0, 3) = tmp(0, 2) * y - tmp(0, 1) * z;
        _jacobianOplusXi(0, 4) = tmp(0, 0) * z - tmp(0, 2) * x;
        _jacobianOplusXi(0, 5) = tmp(0, 1) * x - tmp(0, 0) * y;
        _jacobianOplusXi(0, 6) = -x * invz * rp;
        _jacobianOplusXi(0, 7) = -f * invz * x * r_2;
        _jacobianOplusXi(0, 8) = -f * invz * x * r_2 * r_2;

        _jacobianOplusXi(1, 0) = tmp(1, 0);
        _jacobianOplusXi(1, 1) = tmp(1, 1);
        _jacobianOplusXi(1, 2) = tmp(1, 2);
        _jacobianOplusXi(1, 3) = tmp(1, 2) * y - tmp(1, 1) * z;
        _jacobianOplusXi(1, 4) = tmp(1, 0) * z - tmp(1, 2) * x;
        _jacobianOplusXi(1, 5) = tmp(1, 1) * x - tmp(1, 0) * y;
        _jacobianOplusXi(1, 6) = -y * invz * rp;
        _jacobianOplusXi(1, 7) = -f * invz * y * r_2;
        _jacobianOplusXi(1, 8) = -f * invz * y * r_2 * r_2;

        _jacobianOplusXi = -_jacobianOplusXi;
    }
};

#endif // BAL_TYPE_H
