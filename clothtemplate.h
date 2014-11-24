#ifndef CLOTH_H
#define CLOTH_H

#include <Eigen/Core>
#include <vector>

class Mesh;

class ClothTemplate
{
public:
    ClothTemplate(std::string &meshFilename);
    ~ClothTemplate();

    void computeMassInv();
    void computeConstants();
    const Mesh &getMesh() const { return *m_; }

    Eigen::MatrixXd mass;
    Eigen::MatrixXd massinv;
    Eigen::Matrix<double,9,6> matC;
    std::vector<Eigen::Matrix<double,4,9> > matAs;
    std::vector<Eigen::Matrix2d> matGs;
    std::vector<Eigen::Matrix<double,2,3> > matEs;

private:
    Mesh *m_;
};

#endif // CLOTH_H
