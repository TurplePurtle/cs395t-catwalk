#include "clothtemplate.h"
#include "mesh.h"

using namespace std;
using namespace Eigen;

ClothTemplate::ClothTemplate(std::string &meshFilename)
{
    m_ = new Mesh(meshFilename);
}

ClothTemplate::~ClothTemplate()
{
    delete m_;
}
