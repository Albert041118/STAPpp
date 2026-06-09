#pragma once

#include "Element.h"

class CH8 : public CElement
{
public:
    CH8();
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
    virtual void Write(COutputter& output);
    virtual void GenerateLocationMatrix();
    virtual void ElementStiffness(double* Matrix);
    virtual void ElementStress(double* stress, double* Displacement);
};

class CBeam3D : public CElement
{
public:
    CBeam3D();
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
    virtual void Write(COutputter& output);
    virtual void GenerateLocationMatrix();
    virtual void ElementStiffness(double* Matrix);
    virtual void ElementStress(double* stress, double* Displacement);
};

class CShell4R : public CElement
{
public:
    CShell4R();
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
    virtual void Write(COutputter& output);
    virtual void GenerateLocationMatrix();
    virtual void ElementStiffness(double* Matrix);
    virtual void ElementStress(double* stress, double* Displacement);
};

class CTieSpring : public CElement
{
private:
    unsigned int DofCount_;
public:
    CTieSpring();
    virtual bool Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList);
    virtual void Write(COutputter& output);
    virtual void GenerateLocationMatrix();
    virtual void ElementStiffness(double* Matrix);
    virtual void ElementStress(double* stress, double* Displacement);
};
