/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*****************************************************************************/

#include "T3.h"

#include <cstdlib>
#include <iomanip>
#include <iostream>

using namespace std;

CT3::CT3()
{
    NEN_ = 3;
    nodes_ = new CNode*[NEN_];

    ND_ = 6;    // ux, uy for 3 nodes
    LocationMatrix_ = new unsigned int[ND_];

    ElementMaterial_ = nullptr;
}

CT3::~CT3()
{
}

bool CT3::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int MSet;
    unsigned int N[3];

    Input >> N[0] >> N[1] >> N[2] >> MSet;

    ElementMaterial_ = dynamic_cast<CT3Material*>(MaterialSets) + MSet - 1;

    for (unsigned int i = 0; i < 3; ++i)
        nodes_[i] = &NodeList[N[i] - 1];

    return true;
}

void CT3::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber
           << setw(9) << nodes_[2]->NodeNumber
           << setw(12) << ElementMaterial_->nset << endl;
}

void CT3::GenerateLocationMatrix()
{
    unsigned int i = 0;
    for (unsigned int N = 0; N < NEN_; N++)
    {
        LocationMatrix_[i++] = nodes_[N]->bcode[0];
        LocationMatrix_[i++] = nodes_[N]->bcode[1];
    }
}

void CT3::ElasticityMatrix(double D[3][3]) const
{
    CT3Material* material = dynamic_cast<CT3Material*>(ElementMaterial_);

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
        cerr << "*** Error *** Invalid T3 material mode " << material->mode
             << ". Use 1 for plane stress or 2 for plane strain." << endl;
        exit(5);
    }
}

bool CT3::StrainDisplacementMatrix(double B[3][6], double& area) const
{
    const double x1 = nodes_[0]->XYZ[0];
    const double y1 = nodes_[0]->XYZ[1];
    const double x2 = nodes_[1]->XYZ[0];
    const double y2 = nodes_[1]->XYZ[1];
    const double x3 = nodes_[2]->XYZ[0];
    const double y3 = nodes_[2]->XYZ[1];

    const double twoArea = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1);
    area = 0.5 * twoArea;

    if (area <= 0.0)
        return false;

    const double b[3] = {y2 - y3, y3 - y1, y1 - y2};
    const double c[3] = {x3 - x2, x1 - x3, x2 - x1};

    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 6; ++j)
            B[i][j] = 0.0;

    for (unsigned int i = 0; i < 3; ++i)
    {
        const unsigned int ux = 2 * i;
        const unsigned int uy = 2 * i + 1;

        B[0][ux] = b[i] / twoArea;
        B[1][uy] = c[i] / twoArea;
        B[2][ux] = c[i] / twoArea;
        B[2][uy] = b[i] / twoArea;
    }

    return true;
}

void CT3::ElementStiffness(double* Matrix)
{
    const unsigned int size = SizeOfStiffnessMatrix();
    for (unsigned int i = 0; i < size; ++i)
        Matrix[i] = 0.0;

    double B[3][6], area;
    if (!StrainDisplacementMatrix(B, area))
    {
        cerr << "*** Error *** Non-positive area in T3 element. "
             << "Check that the element nodes are ordered counter-clockwise." << endl;
        exit(5);
    }

    double D[3][3];
    ElasticityMatrix(D);

    CT3Material* material = dynamic_cast<CT3Material*>(ElementMaterial_);
    const double thickness = material->thickness;

    double DB[3][6];
    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 6; ++j)
        {
            DB[i][j] = 0.0;
            for (unsigned int k = 0; k < 3; ++k)
                DB[i][j] += D[i][k] * B[k][j];
        }

    double K[6][6];
    for (unsigned int i = 0; i < 6; ++i)
        for (unsigned int j = 0; j < 6; ++j)
        {
            K[i][j] = 0.0;
            for (unsigned int k = 0; k < 3; ++k)
                K[i][j] += B[k][i] * DB[k][j] * area * thickness;
        }

    for (unsigned int j = 0; j < ND_; ++j)
    {
        const unsigned int offset = (j + 1) * j / 2;
        for (unsigned int i = 0; i <= j; ++i)
            Matrix[offset + j - i] = K[i][j];
    }
}

void CT3::ElementStress(double* stress, double* Displacement)
{
    for (unsigned int i = 0; i < 3; ++i)
        stress[i] = 0.0;

    double u[6];
    for (unsigned int i = 0; i < 6; ++i)
        u[i] = LocationMatrix_[i] ? Displacement[LocationMatrix_[i] - 1] : 0.0;

    double B[3][6], area;
    if (!StrainDisplacementMatrix(B, area))
    {
        cerr << "*** Error *** Non-positive area in T3 element. "
             << "Check that the element nodes are ordered counter-clockwise." << endl;
        exit(5);
    }

    double D[3][3];
    ElasticityMatrix(D);

    double strain[3] = {0.0, 0.0, 0.0};
    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 6; ++j)
            strain[i] += B[i][j] * u[j];

    for (unsigned int i = 0; i < 3; ++i)
        for (unsigned int j = 0; j < 3; ++j)
            stress[i] += D[i][j] * strain[j];
}
