#include "simulation.h"
#include <QGLWidget>
#include "simparameters.h"
#include <iostream>
#include <Eigen/Geometry>
#include <QDebug>
#include "SOIL.h"
#include "rigidbodytemplate.h"
#include "rigidbodyinstance.h"
#include "vectormath.h"
#include <Eigen/Dense>
#include "mesh.h"
#include "signeddistancefield.h"
#include "clothtemplate.h"
#include "clothinstance.h"
#include <Eigen/LU>

const double PI = 3.1415926535898;

using namespace Eigen;
using namespace std;

Simulation::Simulation(const SimParameters &params) : params_(params), time_(0), floorTex_(0)
{
    loadRigidBodies();
    bodyInstance_ = NULL;
    string clothname("resources/square.obj");
    clothTemplate_ = new ClothTemplate(clothname);
    cloth_ = NULL;
    clearScene();
}

Simulation::~Simulation()
{
    delete bodyInstance_;
    delete bodyTemplate_;
    delete cloth_;
    delete clothTemplate_;
}

void Simulation::initializeGL()
{
    loadFloorTexture();
}

void Simulation::loadRigidBodies()
{
    string objname("resources/2by4.obj");
    bodyTemplate_ = new RigidBodyTemplate(objname);
    string sdfname("resources/2by4.sdf");
    bodyTemplate_->computeSDF(sdfname.c_str());
}

void Simulation::loadFloorTexture()
{
    floorTex_ = SOIL_load_OGL_texture("resources/grid.jpg", SOIL_LOAD_AUTO, SOIL_CREATE_NEW_ID, SOIL_FLAG_INVERT_Y |  SOIL_FLAG_NTSC_SAFE_RGB | SOIL_FLAG_COMPRESS_TO_DXT | SOIL_FLAG_MIPMAPS);
    if(floorTex_ != 0)
    {
        glBindTexture(GL_TEXTURE_2D, floorTex_);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
    }
}


void Simulation::renderFloor()
{
    renderLock_.lock();

    glColor4f(1.0, 1.0, 1.0, 1.0);

    if(floorTex_)
    {
        glBindTexture(GL_TEXTURE_2D, floorTex_);
        glEnable(GL_TEXTURE_2D);
    }
    else
    {
        glColor3f(0.5, 0.5, 0.5);
        glDisable(GL_TEXTURE_2D);
    }

    double texsize = 5.0;
    double gridsize = 1000.0;

    double texmax = gridsize/texsize;

    Vector3d tangent1(1.0,0,0);
    Vector3d tangent2(0,-1.0,0);

    Vector3d corner;

    glBegin(GL_QUADS);
    {
        glTexCoord2f(texmax, texmax);
        glNormal3f(0, 0, 1.0);
        corner = -gridsize*tangent1 + gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(texmax, -texmax);
        glNormal3f(0, 0, 1.0);
        corner = -gridsize*tangent1 - gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(-texmax, -texmax);
        glNormal3f(0, 0, 1.0);
        corner = gridsize*tangent1 - gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);

        glTexCoord2f(-texmax, texmax);
        glNormal3f(0, 0, 1.0);
        corner = gridsize*tangent1 + gridsize*tangent2;
        glVertex3f(corner[0], corner[1], corner[2]);
    }
    glDisable(GL_TEXTURE_2D);
    glEnd();
    renderLock_.unlock();
}

void Simulation::renderObjects()
{
    renderLock_.lock();
    {
        bodyInstance_->render();
        cloth_->computeNormals();
        cloth_->render();
    }
    renderLock_.unlock();
}

void Simulation::takeSimulationStep()
{
    const double dt = params_.timeStep;
    time_ += dt;

    VectorXd Fcloth;
    Fcloth.resize(cloth_->q.size());
    Fcloth.setZero();

    renderLock_.lock();
    {
        bodyInstance_->c += dt*bodyInstance_->cvel;
        cloth_->q += cloth_->v * dt;
    }
    renderLock_.unlock();

    bodyInstance_->theta = VectorMath::axisAngle(VectorMath::rotationMatrix(dt*bodyInstance_->w)*VectorMath::rotationMatrix(bodyInstance_->theta));

    computeClothForces(Fcloth);
    cloth_->v += cloth_->getTemplate().massinv * Fcloth * dt;

    if (params_.pinCorner)
        cloth_->v.segment<3>(0).setZero();
}

void Simulation::computeClothForces(VectorXd &F)
{
    if (params_.activeForces & SimParameters::F_GRAVITY)
        for (int i=0; i<cloth_->q.size()/3; i++)
            F(3*i+2) += cloth_->getTemplate().mass(3*i+2,3*i+2) * params_.gravityG;

    if (params_.activeForces & SimParameters::F_DAMPING)
        F -= params_.dampingCoeff * cloth_->getTemplate().mass * cloth_->v;

    if (params_.activeForces & SimParameters::F_STRETCHING)
    {
        int n = clothTemplate_->getMesh().getNumFaces();

        for (int i=0; i<n; i++)
        {
            Vector3i verts = cloth_->getTemplate().getMesh().getFace(i);
            Vector3d pi = cloth_->q.segment<3>(3*verts[0]);
            Vector3d e1 = cloth_->q.segment<3>(3*verts[1]) - pi;
            Vector3d e2 = cloth_->q.segment<3>(3*verts[2]) - pi;

            Matrix<double,2,3> e;
            e << e1.transpose(), e2.transpose();

            Matrix2d g = e * e.transpose();

            Matrix<double,6,4> DifG;
            DifG.block<3,1>(0,0) = 2*e1;
            DifG.block<3,1>(3,0).setZero();
            DifG.block<3,1>(0,1) = e2;
            DifG.block<3,1>(3,1) = e1;
            DifG.block<3,1>(0,2) = e2;
            DifG.block<3,1>(3,2) = e1;
            DifG.block<3,1>(0,3).setZero();
            DifG.block<3,1>(3,3) = 2*e2;

            Matrix2d dg = g - clothTemplate_->matGs[i];
            Matrix<double,2,3> B = clothTemplate_->matGs[i].inverse() * clothTemplate_->matEs[i];
            Matrix3d strainTensor = 0.5 * B.transpose()*dg*B;
            Matrix<double,9,1> strainTensorv;
            strainTensorv.segment<3>(0) = strainTensor.block<1,3>(0,0);
            strainTensorv.segment<3>(3) = strainTensor.block<1,3>(1,0);
            strainTensorv.segment<3>(6) = strainTensor.block<1,3>(2,0);

            Matrix<double,9,1> Ftri = -2*params_.stretchingK*clothTemplate_->getMesh().getFaceArea(i)*clothTemplate_->matC*DifG*clothTemplate_->matAs[i]*strainTensorv;

            F.segment<3>(3*verts[0]) += Ftri.segment<3>(0);
            F.segment<3>(3*verts[1]) += Ftri.segment<3>(3);
            F.segment<3>(3*verts[2]) += Ftri.segment<3>(6);
        }
    }

    if (params_.activeForces & SimParameters::F_BENDING)
    {
        int n = cloth_->getTemplate().hinges.size();

        for (int i=0; i<n; i++)
        {
            const ClothHinge &h = cloth_->getTemplate().hinges[i];

            const Vector3d pi = cloth_->q.segment<3>(3*h.i);
            const Vector3d pj = cloth_->q.segment<3>(3*h.j);
            const Vector3d pk = cloth_->q.segment<3>(3*h.k);
            const Vector3d pl = cloth_->q.segment<3>(3*h.l);
            const Vector3d n0 = (pj - pi).cross(pk - pi);
            const Vector3d n1 = (pl - pi).cross(pj - pi);
            const Vector3d n0xn1 = n0.cross(n1);

            if (n0xn1.squaredNorm() != 0.0)
            {
                double n0xn1_norm = n0xn1.norm();
                double theta = 2*atan2(n0xn1_norm, n0.norm()*n1.norm() + n0.dot(n1));
                double k = -2*params_.bendingK*h.coeff*theta;
                Matrix<double,1,3> C1 = (n0xn1 / n0xn1_norm).cross(n1 / n1.squaredNorm()).transpose();
                Matrix<double,1,3> C0 = (n0xn1 / n0xn1_norm).cross(n0 / n0.squaredNorm()).transpose();

                Matrix3d Dn0 = VectorMath::crossProductMatrix(pk - pj);
                Matrix3d Dn1 = VectorMath::crossProductMatrix(pj - pl);
                F.segment<3>(3*h.i) += k*(C1*Dn1 - C0*Dn0);

                Dn0 = VectorMath::crossProductMatrix(pi - pk);
                Dn1 = VectorMath::crossProductMatrix(pl - pi);
                F.segment<3>(3*h.j) += k*(C1*Dn1 - C0*Dn0);

                Dn0 = VectorMath::crossProductMatrix(pj - pi);
                F.segment<3>(3*h.k) += k*(-C0*Dn0);

                Dn1 = VectorMath::crossProductMatrix(pi - pj);
                F.segment<3>(3*h.l) += k*(C1*Dn1);
            }
        }
    }
}

void Simulation::clearScene()
{
    renderLock_.lock();
    {
        delete bodyInstance_;
        Vector3d pos(5, 0, 3);
        Vector3d zero(0,0,0);
        bodyInstance_ = new RigidBodyInstance(*bodyTemplate_, pos, zero, 1.0);

        delete cloth_;
        Vector3d trans(5, 0, 4);
        cloth_ = new ClothInstance(*clothTemplate_, trans);
    }
    renderLock_.unlock();
}

void Simulation::accelerateBody(double vx, double vy, double vz, double wx, double wy, double wz)
{
    bodyInstance_->cvel += Vector3d(vx,vy,vz);
    bodyInstance_->w += Vector3d(wx,wy,wz);
}
