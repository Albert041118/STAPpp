/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! Four-node isoparametric quadrilateral element for 2D elasticity
class CQ4 : public CElement
{
public:

//! Constructor
    CQ4();

//! Destructor
    ~CQ4();

//! Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//! Write element data to stream
    virtual void Write(COutputter& output);

//! Generate location matrix using ux and uy of each node only
    virtual void GenerateLocationMatrix();

//! Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

//! Calculate averaged element stress: sigma_x, sigma_y, tau_xy
    virtual void ElementStress(double* stress, double* Displacement);

private:
//! Build the elasticity matrix for plane stress or plane strain
    void ElasticityMatrix(double D[3][3]) const;

//! Calculate shape function derivatives with respect to xi and eta
    void ShapeFunctionDerivatives(double xi, double eta, double dNdxi[4], double dNdeta[4]) const;

//! Calculate the strain-displacement matrix
    bool StrainDisplacementMatrix(double xi, double eta, double B[3][8], double& detJ) const;
};
