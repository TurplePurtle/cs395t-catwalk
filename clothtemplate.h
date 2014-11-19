#ifndef CLOTH_H
#define CLOTH_H

#include <Eigen/Core>

class Mesh;

class ClothTemplate
{
public:
    ClothTemplate(std::string &meshFilename);
    ~ClothTemplate();

    void computeMassInv();
    const Mesh &getMesh() const { return *m_; }

    Eigen::MatrixXd massinv;

private:
    Mesh *m_;
};

#endif // CLOTH_H
