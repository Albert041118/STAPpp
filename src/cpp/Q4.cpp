/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*****************************************************************************/

#include "Q4.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace std;

CQ4::CQ4()
{
    NEN_ = 4;
    nodes_ = new CNode*[NEN_];

    ND_ = 8;    // ux, uy for 4 nodes
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

CQ4::~CQ4()
{
}

bool CQ4::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;
    unsigned int N[4];

    Input >> N[0] >> N[1] >> N[2] >> N[3] >> MSet;

    ElementMaterial_ = dynamic_cast<CQ4Material*>(MaterialSets) + MSet - 1;

    for (unsigned int i = 0; i < 4; ++i)
        nodes_[i] = &NodeList[N[i] - 1];

    return true;
}

void CQ4::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber
           << setw(9) << nodes_[2]->NodeNumber
           << setw(9) << nodes_[3]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

void CQ4::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
    {
        LocationMatrix_[i++] = nodes_[N]->bcode[0];
        LocationMatrix_[i++] = nodes_[N]->bcode[1];
    }
}

void CQ4::ElasticityMatrix(double D[3][3]) const
{
    CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);

    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            D[i][j] = 0.0;

    const double E = material->E;
    const double nu = material->nu;

    if (material->mode == 1) // plane stress
    {
        const double c = E / (1.0 - nu * nu);
        D[0][0] = c;
        D[0][1] = c * nu;
        D[1][0] = c * nu;
        D[1][1] = c;
        D[2][2] = c * (1.0 - nu) / 2.0;
    }
    else if (material->mode == 2) // plane strain
    {
        const double c = E / ((1.0 + nu) * (1.0 - 2.0 * nu));
        D[0][0] = c * (1.0 - nu);
        D[0][1] = c * nu;
        D[1][0] = c * nu;
        D[1][1] = c * (1.0 - nu);
        D[2][2] = c * (1.0 - 2.0 * nu) / 2.0;
    }
    else
    {
        cerr << "*** Error *** Invalid Q4 material mode " << material->mode
             << ". Use 1 for plane stress or 2 for plane strain." << endl;
        exit(5);
    }
}

void CQ4::ShapeFunctionDerivatives(double xi, double eta, double dNdxi[4], double dNdeta[4]) const
{
    dNdxi[0] = -0.25 * (1.0 - eta);
    dNdxi[1] =  0.25 * (1.0 - eta);
    dNdxi[2] =  0.25 * (1.0 + eta);
    dNdxi[3] = -0.25 * (1.0 + eta);

    dNdeta[0] = -0.25 * (1.0 - xi);
    dNdeta[1] = -0.25 * (1.0 + xi);
    dNdeta[2] =  0.25 * (1.0 + xi);
    dNdeta[3] =  0.25 * (1.0 - xi);
}

bool CQ4::StrainDisplacementMatrix(double xi, double eta, double B[3][8], double& detJ) const
{
    double dNdxi[4], dNdeta[4];
    ShapeFunctionDerivatives(xi, eta, dNdxi, dNdeta);

    double dx_dxi = 0.0, dy_dxi = 0.0, dx_deta = 0.0, dy_deta = 0.0;
    for (unsigned int i = 0; i < 4; ++i)
    {
        dx_dxi  += dNdxi[i]  * nodes_[i]->XYZ[0];
        dy_dxi  += dNdxi[i]  * nodes_[i]->XYZ[1];
        dx_deta += dNdeta[i] * nodes_[i]->XYZ[0];
        dy_deta += dNdeta[i] * nodes_[i]->XYZ[1];
    }

    detJ = dx_dxi * dy_deta - dx_deta * dy_dxi;

    if (detJ <= 0.0)
        return false;

    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 8; ++j)
            B[i][j] = 0.0;

    for (unsigned int i = 0; i < 4; ++i)
    {
        const double dNdx = ( dy_deta * dNdxi[i] - dy_dxi * dNdeta[i]) / detJ;
        const double dNdy = (-dx_deta * dNdxi[i] + dx_dxi * dNdeta[i]) / detJ;

        const unsigned int ux = 2 * i;
        const unsigned int uy = 2 * i + 1;

        B[0][ux] = dNdx;
        B[1][uy] = dNdy;
        B[2][ux] = dNdy;
        B[2][uy] = dNdx;
    }

    return true;
}

void CQ4::ElementStiffness(double* Matrix)
{
    const unsigned int size = SizeOfStiffnessMatrix();
    for (unsigned int i = 0; i < size; ++i)
        Matrix[i] = 0.0;

    double K[8][8];
    for (unsigned int i = 0; i < 8; ++i)
        for (unsigned int j = 0; j < 8; ++j)
            K[i][j] = 0.0;

    double D[3][3];
    ElasticityMatrix(D);

    CQ4Material* material = dynamic_cast<CQ4Material*>(ElementMaterial_);
    const double thickness = material->thickness;
    const double a = 1.0 / sqrt(3.0);
    const double gauss[2] = {-a, a};

    for (unsigned int ig = 0; ig < 2; ++ig)
    {
        for (unsigned int jg = 0; jg < 2; ++jg)
        {
            double B[3][8], detJ;
            if (!StrainDisplacementMatrix(gauss[ig], gauss[jg], B, detJ))
            {
                cerr << "*** Error *** Non-positive Jacobian determinant in Q4 element. "
                     << "Check that the element nodes are ordered counter-clockwise." << endl;
                exit(5);
            }

            double DB[3][8];
            for (unsigned int i = 0; i < 3; ++i)
                for (unsigned int j = 0; j < 8; ++j)
                {
                    DB[i][j] = 0.0;
                    for (unsigned int k = 0; k < 3; ++k)
                        DB[i][j] += D[i][k] * B[k][j];
                }

            for (unsigned int i = 0; i < 8; ++i)
                for (unsigned int j = 0; j < 8; ++j)
                    for (unsigned int k = 0; k < 3; ++k)
                        K[i][j] += B[k][i] * DB[k][j] * detJ * thickness;
        }
    }

    for (unsigned int j = 0; j < ND_; ++j)
    {
        const unsigned int offset = (j + 1) * j / 2;
        for (unsigned int i = 0; i <= j; ++i)
            Matrix[offset + j - i] = K[i][j];
    }
}

void CQ4::ElementStress(double* stress, double* Displacement)
{
    for (unsigned int i = 0; i < 3; ++i)
        stress[i] = 0.0;

    double u[8];
    for (unsigned int i = 0; i < 8; ++i)
        u[i] = LocationMatrix_[i] ? Displacement[LocationMatrix_[i] - 1] : 0.0;

    double D[3][3];
    ElasticityMatrix(D);

    double B[3][8], detJ;
    if (!StrainDisplacementMatrix(0.0, 0.0, B, detJ))
    {
        cerr << "*** Error *** Non-positive Jacobian determinant in Q4 element. "
             << "Check that the element nodes are ordered counter-clockwise." << endl;
        exit(5);
    }

    double strain[3] = {0.0, 0.0, 0.0};
    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 8; ++j)
            strain[i] += B[i][j] * u[j];

    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            stress[i] += D[i][j] * strain[j];
}
