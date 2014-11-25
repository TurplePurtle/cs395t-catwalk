#ifndef CLOTH_H
#define CLOTH_H

#include <Eigen/Core>
#include <vector>

class Mesh;

struct ClothHinge
{
    double coeff;
    int i;
    int j;
    int k;
    int l;
    int F0;
    int F1;

    ClothHinge(int i, int j, int k, int l, int F0, int F1, double coeff)
        : i(i), j(j), k(k), l(l), F0(F0), F1(F1), coeff(coeff)
    {}
};

class ClothTemplate
{
public:
    ClothTemplate(std::string &meshFilename);
    ~ClothTemplate();

    void computeMassInv();
    void computeConstants();
    void buildHingeList();
    const Mesh &getMesh() const { return *m_; }

    Eigen::MatrixXd mass;
    Eigen::MatrixXd massinv;
    Eigen::Matrix<double,9,6> matC;
    std::vector<Eigen::Matrix<double,4,9> > matAs;
    std::vector<Eigen::Matrix2d> matGs;
    std::vector<Eigen::Matrix<double,2,3> > matEs;

    std::vector<ClothHinge> hinges;

private:
    Mesh *m_;
};

#endif // CLOTH_H
