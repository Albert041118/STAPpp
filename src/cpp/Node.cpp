/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "Node.h"

CNode::CNode(double X, double Y, double Z)
{
    XYZ[0] = X;		// Coordinates of the node
    XYZ[1] = Y;
    XYZ[2] = Z;
    
    bcode[0] = 0;	// Boundary codes
    bcode[1] = 0;
    bcode[2] = 0;
    bcode[3] = 1;
    bcode[4] = 1;
    bcode[5] = 1;
};

//	Read element data from stream Input
bool CNode::Read(ifstream& Input)
{
    string line;
    do {
        if (!std::getline(Input, line))
            return false;
    } while (line.find_first_not_of(" \t\r\n") == string::npos);

    std::stringstream ss(line);
    std::vector<double> v;
    double x;
    while (ss >> x)
        v.push_back(x);

    if (v.size() >= 10)
    {
        NodeNumber = static_cast<unsigned int>(v[0]);
        for (unsigned int i = 0; i < NDF; ++i)
            bcode[i] = static_cast<unsigned int>(v[1+i]);
        XYZ[0] = v[7]; XYZ[1] = v[8]; XYZ[2] = v[9];
    }
    else if (v.size() >= 7)
    {
        // Legacy STAP++ .dat: node b1 b2 b3 x y z.
        // Rotational DOFs are fixed by default so old Bar/Q4/T3 cases do not gain
        // zero-stiffness equations after NDF is extended to 6.
        NodeNumber = static_cast<unsigned int>(v[0]);
        bcode[0] = static_cast<unsigned int>(v[1]);
        bcode[1] = static_cast<unsigned int>(v[2]);
        bcode[2] = static_cast<unsigned int>(v[3]);
        bcode[3] = bcode[4] = bcode[5] = 1;
        XYZ[0] = v[4]; XYZ[1] = v[5]; XYZ[2] = v[6];
    }
    else
    {
        cerr << "*** Error *** Invalid node input line: " << line << endl;
        return false;
    }

	return true;
}

//	Output nodal point data to stream
void CNode::Write(COutputter& output)
{
	output << setw(9) << NodeNumber << setw(5) << bcode[0] << setw(5) << bcode[1] << setw(5) << bcode[2]
		   << setw(18) << XYZ[0] << setw(15) << XYZ[1] << setw(15) << XYZ[2] << endl;
}

//	Output equation numbers of nodal point to stream
void CNode::WriteEquationNo(COutputter& output)
{
	output << setw(9) << NodeNumber << "       ";

	for (unsigned int dof = 0; dof < CNode::NDF; dof++)	// Loop over for DOFs of node np
	{
		output << setw(5) << bcode[dof];
	}

	output << endl;
}

//	Write nodal displacement
void CNode::WriteNodalDisplacement(COutputter& output, double* Displacement)
{
	output << setw(5) << NodeNumber << "        ";

	for (unsigned int j = 0; j < NDF; j++)
	{
		if (bcode[j] == 0)
		{
			output << setw(18) << 0.0;
		}
		else
		{
			output << setw(18) << Displacement[bcode[j] - 1];
		}
	}

	output << endl;
}
