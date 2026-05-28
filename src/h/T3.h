/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*****************************************************************************/

#pragma once

#include "Element.h"

using namespace std;

//! Three-node constant strain triangular element for 2D elasticity
class CT3 : public CElement
{
public:

//! Constructor
    CT3();

//! Destructor
    ~CT3();

//! Read element data from stream Input
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);

//! Write element data to stream
    virtual void Write(COutputter& output);

//! Generate location matrix using ux and uy of each node only
    virtual void GenerateLocationMatrix();

//! Calculate element stiffness matrix
    virtual void ElementStiffness(double* Matrix);

//! Calculate constant element stress: sigma_x, sigma_y, tau_xy
    virtual void ElementStress(double* stress, double* Displacement);

private:
//! Build the elasticity matrix for plane stress or plane strain
    void ElasticityMatrix(double D[3][3]) const;

//! Calculate the constant strain-displacement matrix and element area
    bool StrainDisplacementMatrix(double B[3][6], double& area) const;
};
