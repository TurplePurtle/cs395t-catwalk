#include "clothtemplate.h"
#include "mesh.h"

using namespace std;
using namespace Eigen;

ClothTemplate::ClothTemplate(std::string &meshFilename)
{
    m_ = new Mesh(meshFilename);
    computeMassInv();
}

ClothTemplate::~ClothTemplate()
{
    delete m_;
}

void ClothTemplate::computeMassInv()
{
    int nv = m_->getNumVerts();
    int nf = m_->getNumFaces();

    mass.resize(3*nv,3*nv);
    massinv.resize(3*nv,3*nv);
    mass.setZero();
    massinv.setZero();

    for (int i=0; i<nf; i++)
    {
        const Vector3i verts = m_->getFace(i);
        double area = m_->getFaceArea(i);

        for (int j=0; j<3; j++)
        {
            mass(3*verts[j]  , 3*verts[j]  ) = area;
            mass(3*verts[j]+1, 3*verts[j]+1) = area;
            mass(3*verts[j]+2, 3*verts[j]+2) = area;
        }
    }

    mass /= 3.0;
    massinv.diagonal() = 1.0 / mass.diagonal().array();
}
