#ifndef CLOTH_H
#define CLOTH_H

#include <Eigen/Core>

class Mesh;

class ClothTemplate
{
public:
    ClothTemplate(std::string &meshFilename);
    ~ClothTemplate();

    const Mesh &getMesh() const { return *m_; }

private:
    Mesh *m_;
};

#endif // CLOTH_H
