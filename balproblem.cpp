#include "balproblem.h"

#include <iostream>
#include <fstream>
#include <unistd.h>
using namespace std;

balProblem::balProblem(const string& filename)//get the datas
{
    ifstream ifs;
    ifs.open(filename, fstream::in);
    if(!ifs.is_open())return;

    int camera_num, point_num, edge_num;
    ifs >> camera_num >> point_num >> edge_num;

    for(int i = 0 ; i < edge_num ; ++i)
    {
        int c, p;
        double uv[2];
        ifs >> c >> p >> uv[0] >> uv[1];
        edges.push_back(balEdge(c, p, uv));
    }

    for(int i = 0 ; i < camera_num ; ++i)
    {
        double data[9];
        for(int j = 0 ; j < 9 ; ++j)
            ifs >> data[j];
        cameras.push_back(balCamera(data));
    }

    for(int i = 0 ; i < point_num ; ++i)
    {
        double data[3];
        ifs >> data[0] >> data[1] >> data[2];
        points.push_back((balPoint(data)));
    }

    cout << "data is ready\nthere are " << cameras.size() << " cameras, " << points.size() << " points and " << edges.size() << " edges\n";
}

void balProblem::buildProblem()
{
    typedef g2o::BlockSolver< g2o::BlockSolverTraits<9,3> > Block;  // pose 维度为 9, landmark 维度为 3
    Block::LinearSolverType* linearSolver = new g2o::LinearSolverCholmod<Block::PoseMatrixType>(); // 线性方程求解器, 使用稀疏计算方法
    dynamic_cast<g2o::LinearSolverCholmod<Block::PoseMatrixType>*>(linearSolver)->setBlockOrdering(true);//排序以保持稀疏性


    Block* solver_ptr = new Block ( linearSolver );     // 矩阵块求解器
    g2o::OptimizationAlgorithmWithHessian* solver = new g2o::OptimizationAlgorithmLevenberg ( solver_ptr );

    optimizer.setAlgorithm ( solver );

    int camera_num = cameras.size();
    int point_num = points.size();
    int edge_num = edges.size();

    int index = 0;
    for(auto data:cameras)
    {
        balCameraVertex* camera = new balCameraVertex();
        camera->setEstimate(data.cam);
        camera->setId(index);
        optimizer.addVertex(camera);
        ++index;
    }

    for(auto data:points)
    {
        bal3DVertex* point = new bal3DVertex();
        point->setEstimate(data.pt);
        point->setId(index);
        point->setMarginalized(true);
        optimizer.addVertex(point);
        ++index;
    }

    index = 0;
    for(auto data:edges)
    {
        balCamPtEdge* edge = new balCamPtEdge();//0: camera; 1:point
        edge->setVertex(0, dynamic_cast< balCameraVertex* >( optimizer.vertex(data.cam) ));
        edge->setVertex(1, dynamic_cast< bal3DVertex* >(optimizer.vertex(data.pt + camera_num)));
        edge->setMeasurement(data.uv);
        edge->setId(index);
        edge->setInformation(Eigen::Matrix2d::Identity());
        optimizer.addEdge(edge);
        ++index;
    }

    cout << "the edge and vertex is ready, let's solve the problem\n";
}

void balProblem::solveProblem(int iter)
{

    optimizer.setVerbose(true);//verbose information during optimization
    optimizer.initializeOptimization();

    writeToPLY("./start.ply");
    optimizer.optimize(iter);

    cout << "mission completed\nnow you can choose the way to show the point cloud\n" ;

}

void balProblem::writeToPLY(const string &filename)
{
    std::ofstream of(filename.c_str());

      int camera_num = cameras.size();
      int point_num = points.size();

      of<< "ply"
        << '\n' << "format ascii 1.0"
        << '\n' << "element vertex " << camera_num + point_num
        << '\n' << "property float x"
        << '\n' << "property float y"
        << '\n' << "property float z"
        << '\n' << "property uchar red"
        << '\n' << "property uchar green"
        << '\n' << "property uchar blue"
        << '\n' << "end_header" << std::endl;

        // Export extrinsic data (i.e. camera centers) as green points.

        for(int i = 0; i < camera_num; ++i){
            const balCameraVertex* mycamera = dynamic_cast< balCameraVertex* >(optimizer.vertex(i));
            Eigen::Matrix<double, 9, 1> cam = mycamera->estimate();
            Eigen::Matrix<double, 6, 1> se3;
            se3 << cam(0, 0), cam(1, 0), cam(2, 0), cam(3, 0), cam(4, 0), cam(5, 0);
            Sophus::SE3d pose = Sophus::SE3d::exp(se3);
            Eigen::Matrix4d mat_T = pose.matrix();

            of << mat_T(0, 3) << ' ' << mat_T(1, 3) << ' ' << mat_T(2, 3)
             << " 0 255 0" << '\n';
        }

        // Export the structure (i.e. 3D Points) as white points.

        for(int i = 0; i < point_num; ++i){
          const bal3DVertex* mypoint = dynamic_cast<const bal3DVertex*>(optimizer.vertex(i + camera_num));
          Eigen::Vector3d point_data = mypoint->estimate();
          for(int j = 0; j < 3; ++j){
            of << point_data(j, 0) << ' ';
          }
          of << "255 255 255\n";
        }
        of.close();
        cout << "now you can open the .ply file\n";
}

void balProblem::showByPangolin()
{
    pangolin::CreateWindowAndBind("Trajectory Viewer", 1024, 768);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    pangolin::OpenGlRenderState s_cam(
             pangolin::ProjectionMatrix(1024, 768, 500, 500, 512, 389, 0.1, 1000),
             pangolin::ModelViewLookAt(0, -0.1, -1.8, 0, 0, 0, 0.0, -1.0, 0.0)
        );

    pangolin::View &d_cam = pangolin::CreateDisplay()
            .SetBounds(0.0, 1.0, pangolin::Attach::Pix(175), 1.0, -1024.0f / 768.0f)
            .SetHandler(new pangolin::Handler3D(s_cam));


    while (pangolin::ShouldQuit() == false) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        d_cam.Activate(s_cam);
        glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

            // points
        glPointSize(1);
        glBegin(GL_POINTS);
        int camera_num = cameras.size();
        for (size_t i = 0; i < points.size(); i++) {
            const bal3DVertex* mypoint = dynamic_cast<const bal3DVertex*>(optimizer.vertex(i + camera_num));
            Eigen::Vector3d point_data = mypoint->estimate();
            glColor3f(1.0, 1.0, 1.0);
            glVertex3d(point_data(0, 0), -point_data(1, 0), -point_data(2, 0) );
        }
        glEnd();

        pangolin::FinishFrame();
        usleep(5000);   // sleep 5 ms
    }

}


