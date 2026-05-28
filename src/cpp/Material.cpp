/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#include "Material.h"

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

//	Read material data from stream Input
bool CBarMaterial::Read(ifstream& Input)
{
	Input >> nset;	// Number of property set

	Input >> E >> Area;	// Young's modulus and section area

	return true;
}

//	Write material data to Stream
void CBarMaterial::Write(COutputter& output)
{
	output << setw(16) << E << setw(16) << Area << endl;
}

//  Read Q4 material data from stream Input
bool CQ4Material::Read(ifstream& Input)
{
    Input >> nset;  // Number of property set
    Input >> E >> nu >> thickness >> mode;

    if (mode != 1 && mode != 2)
    {
        cerr << "*** Error *** Invalid Q4 material mode " << mode
             << ". Use 1 for plane stress or 2 for plane strain." << endl;
        return false;
    }

    if (thickness <= 0.0)
    {
        cerr << "*** Error *** Q4 material thickness must be positive." << endl;
        return false;
    }

    return true;
}

//  Write Q4 material data to Stream
void CQ4Material::Write(COutputter& output)
{
    output << setw(16) << E
           << setw(16) << nu
           << setw(16) << thickness
           << setw(8) << mode << endl;
}

//  Read T3 material data from stream Input
bool CT3Material::Read(ifstream& Input)
{
    Input >> nset;  // Number of property set
    Input >> E >> nu >> thickness >> mode;

    if (mode != 1 && mode != 2)
    {
        cerr << "*** Error *** Invalid T3 material mode " << mode
             << ". Use 1 for plane stress or 2 for plane strain." << endl;
        return false;
    }

    if (thickness <= 0.0)
    {
        cerr << "*** Error *** T3 material thickness must be positive." << endl;
        return false;
    }

    return true;
}

//  Write T3 material data to Stream
void CT3Material::Write(COutputter& output)
{
    output << setw(16) << E
           << setw(16) << nu
           << setw(16) << thickness
           << setw(8) << mode << endl;
}
