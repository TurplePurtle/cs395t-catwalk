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
