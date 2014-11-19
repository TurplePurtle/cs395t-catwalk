#ifndef CLOTHINSTANCE_H
#define CLOTHINSTANCE_H

#include <Eigen/Core>
#include <vector>

class ClothTemplate;

class ClothInstance
{
public:
    ClothInstance(ClothTemplate &temp, const Eigen::Vector3d &trans);

    void computeNormals();
    void render();

    const ClothTemplate &getTemplate() { return temp_; }

private:
    const ClothTemplate &temp_;
    Eigen::VectorXd q_;
    Eigen::VectorXd v_;
    Eigen::VectorXd vertNormals_;
    std::vector<Eigen::Vector3d> faceNormals_;
    std::vector<double> faceAreas_;
};

#endif // CLOTHINSTANCE_H
