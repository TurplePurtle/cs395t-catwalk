#include "clothtemplate.h"
#include "mesh.h"
#include <Eigen/LU>

using namespace std;
using namespace Eigen;

template<typename Derived>
inline bool is_finite(const Eigen::MatrixBase<Derived>& x)
{
    return ( (x - x).array() == (x - x).array()).all();
}

ClothTemplate::ClothTemplate(std::string &meshFilename)
{
    m_ = new Mesh(meshFilename);
    computeMassInv();
    computeConstants();
    buildHingeList();
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

void ClothTemplate::computeConstants()
{
    matC.block<3,3>(0,0) = -Matrix3d::Identity();
    matC.block<3,3>(0,3) = -Matrix3d::Identity();
    matC.block<6,6>(3,0) = MatrixXd::Identity(6,6);

    int n = m_->getNumFaces();

    matAs.clear();
    matGs.clear();
    matEs.clear();
    matAs.reserve(n);
    matGs.reserve(n);
    matEs.reserve(n);

    for (int i=0; i<n; i++)
    {
        Vector3i verts = m_->getFace(i);
        Vector3d pi = m_->getVert(verts[0]);
        Vector3d e1 = m_->getVert(verts[1]) - pi;
        Vector3d e2 = m_->getVert(verts[2]) - pi;

        Matrix<double,2,3> e;
        e << e1.transpose(), e2.transpose();
        assert(is_finite(e));
        matEs.push_back(e);

        Matrix2d g;
        g << e1.dot(e1), e1.dot(e2), e1.dot(e2), e2.dot(e2);
        assert(is_finite(g));
        matGs.push_back(g);

        Matrix<double,2,3> B;
        B = g.inverse() * e;

        Matrix<double,9,4> A;
        A <<    B(0,0)*B(0,0), B(0,0)*B(1,0), B(0,0)*B(1,0), B(1,0)*B(1,0),
                B(0,0)*B(0,1), B(0,0)*B(1,1), B(1,0)*B(0,1), B(1,0)*B(1,1),
                B(0,0)*B(0,2), B(0,0)*B(1,2), B(1,0)*B(0,2), B(1,0)*B(1,2),
                B(0,0)*B(0,1), B(1,0)*B(0,1), B(0,0)*B(1,1), B(1,0)*B(1,1),
                B(0,1)*B(0,1), B(0,1)*B(1,1), B(0,1)*B(1,1), B(1,1)*B(1,1),
                B(0,1)*B(0,2), B(0,1)*B(1,2), B(1,1)*B(0,2), B(1,1)*B(1,2),
                B(0,0)*B(0,2), B(1,0)*B(0,2), B(0,0)*B(1,2), B(1,0)*B(1,2),
                B(0,1)*B(0,2), B(1,1)*B(0,2), B(0,1)*B(1,2), B(1,1)*B(1,2),
                B(0,2)*B(0,2), B(0,2)*B(1,2), B(0,2)*B(1,2), B(1,2)*B(1,2);

        assert(is_finite(A));
        matAs.push_back(A.transpose());
    }
}

void ClothTemplate::buildHingeList()
{
    int n = m_->getNumFaces();

    hinges.clear();
    hinges.reserve(n);

    for (int i=0; i<n; i++)
    {
        for (int j=i+1; j<n; j++)
        {
            Vector3i vs0 = m_->getFace(i);
            Vector3i vs1 = m_->getFace(j);

            int cvs[3];
            int common = 0;

            for (int i1=0; i1<3; i1++)
                for (int i2=0; i2<3; i2++)
                    if (vs0[i1] == vs1[i2] && vs0[i1] >= 0)
                    {
                        cvs[common++] = vs0[i1];
                        vs0[i1] = -1;
                        vs1[i2] = -1;
                    }

            assert(common < 3);

            if (common == 2)
            {
                int vk = -1;
                for (int i1=0; i1<3; i1++)
                    if (vs0[i1] != -1)
                    {
                        vk = vs0[i1];
                        break;
                    }
                int vl = -1;
                for (int i1=0; i1<3; i1++)
                    if (vs1[i1] != -1)
                    {
                        vl = vs1[i1];
                        break;
                    }
                double coeff = (m_->getVert(cvs[0]) - m_->getVert(cvs[1])).squaredNorm() / (m_->getFaceArea(i) + m_->getFaceArea(j));
                hinges.push_back(ClothHinge(cvs[0], cvs[1], vk, vl, i, j, coeff));
            }
        }
    }
}
