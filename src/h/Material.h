/*****************************************************************************/
/*  STAP++ : A C++ FEM code sharing the same input data file with STAP90     */
/*     Computational Dynamics Laboratory                                     */
/*     School of Aerospace Engineering, Tsinghua University                  */
/*                                                                           */
/*     Release 1.11, November 22, 2017                                       */
/*                                                                           */
/*     http://www.comdyn.cn/                                                 */
/*****************************************************************************/

#pragma once

#include "Outputter.h"

using namespace std;

//!	Material base class which only define one data member
/*!	All type of material classes should be derived from this base class */
class CMaterial
{
public:

	unsigned int nset;	//!< Number of set
	
	double E;  //!< Young's modulus

public:

//! Virtual deconstructor
    virtual ~CMaterial() {};

//!	Read material data from stream Input
	virtual bool Read(ifstream& Input) = 0;

//!	Write material data to Stream
    virtual void Write(COutputter& output) = 0;

};

//!	Material class for bar element
class CBarMaterial : public CMaterial
{
public:

	double Area;	//!< Sectional area of a bar element

public:
	
//!	Read material data from stream Input
	virtual bool Read(ifstream& Input);

//!	Write material data to Stream
	virtual void Write(COutputter& output);
};

//! Material class for Q4 plane elasticity element
class CQ4Material : public CMaterial
{
public:

    double nu;          //!< Poisson's ratio
    double thickness;   //!< Element thickness
    unsigned int mode;  //!< 1: plane stress, 2: plane strain

public:

//! Read material data from stream Input
    virtual bool Read(ifstream& Input);

//! Write material data to Stream
    virtual void Write(COutputter& output);
};

//! Material class for T3 plane elasticity element
class CT3Material : public CMaterial
{
public:

    double nu;          //!< Poisson's ratio
    double thickness;   //!< Element thickness
    unsigned int mode;  //!< 1: plane stress, 2: plane strain

public:

//! Read material data from stream Input
    virtual bool Read(ifstream& Input);

//! Write material data to Stream
    virtual void Write(COutputter& output);
};

//! Generic material/section data for Bridge H8, Beam and Shell elements
class CBridgeMaterial : public CMaterial
{
public:
    double nu;
    double p[8];   //!< section parameters: shell t; beam A/Iy/Iz/J/orient; H8 unused

public:
    virtual bool Read(ifstream& Input);
    virtual void Write(COutputter& output);
};
