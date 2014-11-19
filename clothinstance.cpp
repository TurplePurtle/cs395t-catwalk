#include "clothinstance.h"
#include "clothtemplate.h"
#include "mesh.h"
#include <QGLWidget>
#include <vector>
#include <Eigen/Geometry>

using namespace Eigen;
using namespace std;

ClothInstance::ClothInstance(ClothTemplate &temp, const Vector3d &trans) : temp_(temp)
{
    const Mesh &mesh = temp_.getMesh();
    const int n = mesh.getNumVerts();

    q_.resize(3*n);

    for (int i=0; i<n; i++)
    {
        q_.segment<3>(3*i) = mesh.getVert(i) + trans;
    }
}

void ClothInstance::computeNormals()
{
    faceNormals_.clear();
    faceAreas_.clear();

    const Mesh &mesh = temp_.getMesh();
    vector<double> verttotalarea;

    vertNormals_.resize(q_.size());
    vertNormals_.setZero();

    for(int i=0; i < (int)q_.size()/3; i++)
    {
        verttotalarea.push_back(0);
    }

    for(int i=0; i < (int)mesh.getNumFaces(); i++)
    {
        Vector3i fverts = mesh.getFace(i);
        Vector3d pts[3];
        for(int j=0; j<3; j++)
            pts[j] = q_.segment<3>(3*fverts[j]);

        Vector3d normal = (pts[1]-pts[0]).cross(pts[2]-pts[0]);
        double norm = normal.norm();
        faceAreas_.push_back(norm/2.0);
        faceNormals_.push_back(normal/norm);
        for(int j=0; j<3; j++)
        {
            vertNormals_.segment<3>(3*fverts[j]) += normal;
            verttotalarea[fverts[j]] += norm;
        }
    }

    for(int i=0; i<(int)q_.size()/3; i++)
        vertNormals_.segment<3>(3*i) /= verttotalarea[i];
}

void ClothInstance::render()
{
    computeNormals();

    glShadeModel(GL_SMOOTH);
    glDisable(GL_TEXTURE_2D);
    glEnable(GL_LIGHTING);
    glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
    glEnable(GL_COLOR_MATERIAL);
    glColor4d(0.6, 0.9, 0.3, 1.0);

    glPushMatrix();
    {
//        GLdouble xform[16];
//        for(int i=0; i<3; i++)
//        {
//            for(int j=0; j<3; j++)
//                xform[4*j+i] = 0;
//            xform[4*i+3] = 0;
//            xform[12+i] = 0;
//        }
//        xform[15] = 1.0;
//        glMultMatrixd(xform);
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_NORMAL_ARRAY);

        glVertexPointer(3, GL_DOUBLE, 0, &q_[0]);
        glNormalPointer(GL_DOUBLE, 0, &vertNormals_[0]);

        glDrawElements(GL_TRIANGLES, q_.size(), GL_UNSIGNED_INT, temp_.getMesh().getFacePointer());

        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    glPopMatrix();
}
