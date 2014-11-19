#ifndef CLOTHINSTANCE_H
#define CLOTHINSTANCE_H

#include <Eigen/Core>

class ClothTemplate;

class ClothInstance
{
public:
    ClothInstance();

    void render();

    const ClothTemplate &getTemplate() { return &template_; }

private:
    const ClothTemplate &template_;
};

#endif // CLOTHINSTANCE_H
