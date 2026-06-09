#include "BridgeElements.h"
#include "Material.h"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

namespace {
struct V3 { double x,y,z; V3(double a=0,double b=0,double c=0):x(a),y(b),z(c){} };
static V3 sub(CNode* a,CNode* b){ return V3(a->XYZ[0]-b->XYZ[0],a->XYZ[1]-b->XYZ[1],a->XYZ[2]-b->XYZ[2]); }
static V3 subp(const V3&a,const V3&b){return V3(a.x-b.x,a.y-b.y,a.z-b.z);}
static V3 addp(const V3&a,const V3&b){return V3(a.x+b.x,a.y+b.y,a.z+b.z);}
static V3 mulp(const V3&a,double s){return V3(a.x*s,a.y*s,a.z*s);}
static double dotp(const V3&a,const V3&b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static V3 crossp(const V3&a,const V3&b){return V3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);}
static double normp(const V3&a){return sqrt(dotp(a,a));}
static V3 unitp(const V3&a){double n=normp(a); if(n<1e-14){cerr<<"zero vector\n"; exit(5);} return mulp(a,1.0/n);}
static V3 nodep(CNode* n){return V3(n->XYZ[0],n->XYZ[1],n->XYZ[2]);}
static void clearFull(vector<double>& K){ for(size_t i=0;i<K.size();++i) K[i]=0.0; }
static void packUpper(const vector<double>& K,unsigned int n,double* M){ for(unsigned int j=0;j<n;j++){unsigned int off=(j+1)*j/2; for(unsigned int i=0;i<=j;i++) M[off+j-i]=K[i*n+j];}}
static void add(vector<double>&K,int n,int i,int j,double v){K[i*n+j]+=v;}
static void matBtDB(vector<double>& K,int nd,const vector<double>& B,int nr,const vector<double>& D,double w){
    vector<double> DB(nr*nd,0.0);
    for(int i=0;i<nr;i++) for(int j=0;j<nd;j++) for(int k=0;k<nr;k++) DB[i*nd+j]+=D[i*nr+k]*B[k*nd+j];
    for(int i=0;i<nd;i++) for(int j=0;j<nd;j++){double s=0; for(int k=0;k<nr;k++) s+=B[k*nd+i]*DB[k*nd+j]; K[i*nd+j]+=s*w;}
}
static bool inv3(const double A[9], double inv[9], double& det){
    det=A[0]*(A[4]*A[8]-A[5]*A[7])-A[1]*(A[3]*A[8]-A[5]*A[6])+A[2]*(A[3]*A[7]-A[4]*A[6]);
    if(fabs(det)<1e-20) return false;
    inv[0]=(A[4]*A[8]-A[5]*A[7])/det; inv[1]=(A[2]*A[7]-A[1]*A[8])/det; inv[2]=(A[1]*A[5]-A[2]*A[4])/det;
    inv[3]=(A[5]*A[6]-A[3]*A[8])/det; inv[4]=(A[0]*A[8]-A[2]*A[6])/det; inv[5]=(A[2]*A[3]-A[0]*A[5])/det;
    inv[6]=(A[3]*A[7]-A[4]*A[6])/det; inv[7]=(A[1]*A[6]-A[0]*A[7])/det; inv[8]=(A[0]*A[4]-A[1]*A[3])/det; return true;
}
static vector<double> D3(double E,double nu){ vector<double>D(36,0.0); double c=E/((1+nu)*(1-2*nu)); double a[6][6]={{1-nu,nu,nu,0,0,0},{nu,1-nu,nu,0,0,0},{nu,nu,1-nu,0,0,0},{0,0,0,(1-2*nu)/2,0,0},{0,0,0,0,(1-2*nu)/2,0},{0,0,0,0,0,(1-2*nu)/2}}; for(int i=0;i<6;i++)for(int j=0;j<6;j++)D[i*6+j]=c*a[i][j]; return D;}
static void q4(double xi,double eta,double N[4],double dn[4][2]){
    N[0]=0.25*(1-xi)*(1-eta);N[1]=0.25*(1+xi)*(1-eta);N[2]=0.25*(1+xi)*(1+eta);N[3]=0.25*(1-xi)*(1+eta);
    dn[0][0]=-0.25*(1-eta);dn[0][1]=-0.25*(1-xi);dn[1][0]=0.25*(1-eta);dn[1][1]=-0.25*(1+xi);dn[2][0]=0.25*(1+eta);dn[2][1]=0.25*(1+xi);dn[3][0]=-0.25*(1+eta);dn[3][1]=0.25*(1-xi);
}
}

CH8::CH8(){ NEN_=8; nodes_=new CNode*[NEN_]; ND_=24; LocationMatrix_=new unsigned int[ND_]; ElementMaterial_=nullptr; }
bool CH8::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList){ unsigned int M,N[8]; for(int i=0;i<8;i++) Input>>N[i]; Input>>M; ElementMaterial_=dynamic_cast<CBridgeMaterial*>(MaterialSets)+M-1; for(int i=0;i<8;i++) nodes_[i]=&NodeList[N[i]-1]; return true; }
void CH8::Write(COutputter& output){ for(int i=0;i<8;i++) output<<setw(9)<<nodes_[i]->NodeNumber; output<<setw(12)<<ElementMaterial_->nset<<endl; }
void CH8::GenerateLocationMatrix(){ unsigned int k=0; for(unsigned int n=0;n<NEN_;n++) for(unsigned int d=0;d<3;d++) LocationMatrix_[k++]=nodes_[n]->bcode[d]; }
void CH8::ElementStiffness(double* Matrix){
    clear(Matrix,SizeOfStiffnessMatrix()); CBridgeMaterial* m=dynamic_cast<CBridgeMaterial*>(ElementMaterial_);
    vector<double>K(24*24,0.0),D=D3(m->E,m->nu); V3 c[8]; for(int i=0;i<8;i++) c[i]=nodep(nodes_[i]);
    int sx[8]={-1,1,1,-1,-1,1,1,-1},sy[8]={-1,-1,1,1,-1,-1,1,1},sz[8]={-1,-1,-1,-1,1,1,1,1}; double gp[2]={-1/sqrt(3.0),1/sqrt(3.0)};
    for(double xi:gp)for(double eta:gp)for(double zeta:gp){ double dn[8][3]; for(int a=0;a<8;a++){dn[a][0]=0.125*sx[a]*(1+sy[a]*eta)*(1+sz[a]*zeta);dn[a][1]=0.125*sy[a]*(1+sx[a]*xi)*(1+sz[a]*zeta);dn[a][2]=0.125*sz[a]*(1+sx[a]*xi)*(1+sy[a]*eta);}
        double J[9]={0}; for(int a=0;a<8;a++){J[0]+=dn[a][0]*c[a].x;J[1]+=dn[a][0]*c[a].y;J[2]+=dn[a][0]*c[a].z;J[3]+=dn[a][1]*c[a].x;J[4]+=dn[a][1]*c[a].y;J[5]+=dn[a][1]*c[a].z;J[6]+=dn[a][2]*c[a].x;J[7]+=dn[a][2]*c[a].y;J[8]+=dn[a][2]*c[a].z;}
        double invJ[9],det; if(!inv3(J,invJ,det)||det<=0){cerr<<"H8 detJ error\n"; exit(5);}
        vector<double>B(6*24,0.0); for(int a=0;a<8;a++){double dx=dn[a][0]*invJ[0]+dn[a][1]*invJ[3]+dn[a][2]*invJ[6]; double dy=dn[a][0]*invJ[1]+dn[a][1]*invJ[4]+dn[a][2]*invJ[7]; double dz=dn[a][0]*invJ[2]+dn[a][1]*invJ[5]+dn[a][2]*invJ[8]; int ix=3*a; B[0*24+ix]=dx;B[1*24+ix+1]=dy;B[2*24+ix+2]=dz;B[3*24+ix]=dy;B[3*24+ix+1]=dx;B[4*24+ix+1]=dz;B[4*24+ix+2]=dy;B[5*24+ix]=dz;B[5*24+ix+2]=dx;}
        matBtDB(K,24,B,6,D,det);
    }
    packUpper(K,ND_,Matrix);
}
void CH8::ElementStress(double* stress,double*){ for(int i=0;i<6;i++) stress[i]=0.0; }

CBeam3D::CBeam3D(){ NEN_=2; nodes_=new CNode*[NEN_]; ND_=12; LocationMatrix_=new unsigned int[ND_]; ElementMaterial_=nullptr; }
bool CBeam3D::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList){ unsigned int M,N1,N2; Input>>N1>>N2>>M; ElementMaterial_=dynamic_cast<CBridgeMaterial*>(MaterialSets)+M-1; nodes_[0]=&NodeList[N1-1]; nodes_[1]=&NodeList[N2-1]; return true; }
void CBeam3D::Write(COutputter& output){ output<<setw(11)<<nodes_[0]->NodeNumber<<setw(9)<<nodes_[1]->NodeNumber<<setw(12)<<ElementMaterial_->nset<<endl; }
void CBeam3D::GenerateLocationMatrix(){ unsigned int k=0; for(unsigned int n=0;n<NEN_;n++) for(unsigned int d=0;d<6;d++) LocationMatrix_[k++]=nodes_[n]->bcode[d]; }
void CBeam3D::ElementStiffness(double* Matrix){
    clear(Matrix,SizeOfStiffnessMatrix()); CBridgeMaterial* m=dynamic_cast<CBridgeMaterial*>(ElementMaterial_);
    double E=m->E,nu=m->nu,A=m->p[0],Iy=m->p[1],Iz=m->p[2],Jp=m->p[3],G=E/(2*(1+nu)); V3 d=sub(nodes_[1],nodes_[0]); double L=normp(d);
    vector<double> k(12*12,0.0); double EA=E*A/L,GJ=G*Jp/L; add(k,12,0,0,EA);add(k,12,6,6,EA);add(k,12,0,6,-EA);add(k,12,6,0,-EA);add(k,12,3,3,GJ);add(k,12,9,9,GJ);add(k,12,3,9,-GJ);add(k,12,9,3,-GJ);
    auto bend=[&](double C,int ids[4],int sgn){double kb[4][4]={{12*C/(L*L*L),sgn*6*C/(L*L),-12*C/(L*L*L),sgn*6*C/(L*L)},{sgn*6*C/(L*L),4*C/L,-sgn*6*C/(L*L),2*C/L},{-12*C/(L*L*L),-sgn*6*C/(L*L),12*C/(L*L*L),-sgn*6*C/(L*L)},{sgn*6*C/(L*L),2*C/L,-sgn*6*C/(L*L),4*C/L}}; for(int i=0;i<4;i++)for(int j=0;j<4;j++)add(k,12,ids[i],ids[j],kb[i][j]);};
    int idz[4]={1,5,7,11}; bend(E*Iz,idz,1); int idy[4]={2,4,8,10}; bend(E*Iy,idy,-1);
    V3 x=unitp(d), ov=V3(m->p[4],m->p[5],m->p[6]); if(normp(ov)<1e-12) ov=V3(0,0,-1); if(fabs(dotp(unitp(ov),x))>0.98) ov=V3(0,1,0); V3 z=unitp(crossp(x,ov)), y=unitp(crossp(z,x)); double R[9]={x.x,x.y,x.z,y.x,y.y,y.z,z.x,z.y,z.z};
    vector<double>K(12*12,0.0); for(int I=0;I<12;I++)for(int J=0;J<12;J++){double sum=0; for(int a=0;a<12;a++)for(int b=0;b<12;b++) if(I/3==a/3&&J/3==b/3) sum+=R[(a%3)*3+(I%3)]*k[a*12+b]*R[(b%3)*3+(J%3)]; K[I*12+J]=sum;}
    packUpper(K,ND_,Matrix);
}
void CBeam3D::ElementStress(double* stress,double*){ for(int i=0;i<12;i++) stress[i]=0.0; }

CShell4R::CShell4R(){ NEN_=4; nodes_=new CNode*[NEN_]; ND_=24; LocationMatrix_=new unsigned int[ND_]; ElementMaterial_=nullptr; }
bool CShell4R::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList){ unsigned int M,N[4]; for(int i=0;i<4;i++) Input>>N[i]; Input>>M; ElementMaterial_=dynamic_cast<CBridgeMaterial*>(MaterialSets)+M-1; for(int i=0;i<4;i++) nodes_[i]=&NodeList[N[i]-1]; return true; }
void CShell4R::Write(COutputter& output){ for(int i=0;i<4;i++) output<<setw(9)<<nodes_[i]->NodeNumber; output<<setw(12)<<ElementMaterial_->nset<<endl; }
void CShell4R::GenerateLocationMatrix(){ unsigned int k=0; for(unsigned int n=0;n<NEN_;n++) for(unsigned int d=0;d<6;d++) LocationMatrix_[k++]=nodes_[n]->bcode[d]; }
void CShell4R::ElementStiffness(double* Matrix){
    clear(Matrix,SizeOfStiffnessMatrix()); CBridgeMaterial* m=dynamic_cast<CBridgeMaterial*>(ElementMaterial_); double E=m->E,nu=m->nu,t=m->p[0],G=E/(2*(1+nu));
    V3 p[4]; for(int i=0;i<4;i++)p[i]=nodep(nodes_[i]); V3 x=unitp(subp(p[1],p[0])), z=unitp(crossp(subp(p[1],p[0]),subp(p[3],p[0]))), y=unitp(crossp(z,x)); V3 o; for(int i=0;i<4;i++)o=addp(o,p[i]); o=mulp(o,0.25); double xy[4][2]; for(int i=0;i<4;i++){V3 d=subp(p[i],o); xy[i][0]=dotp(d,x); xy[i][1]=dotp(d,y);}
    vector<double>Kloc(24*24,0.0), Dm={E*t/(1-nu*nu),E*t*nu/(1-nu*nu),0,E*t*nu/(1-nu*nu),E*t/(1-nu*nu),0,0,0,E*t/(2*(1+nu))}; double cb=E*t*t*t/(12*(1-nu*nu)); vector<double>Db={cb,cb*nu,0,cb*nu,cb,0,0,0,cb*(1-nu)/2}; vector<double>Ds={(5.0/6.0)*G*t,0,0,(5.0/6.0)*G*t}; double gp[2]={-1/sqrt(3.0),1/sqrt(3.0)}, area=0;
    for(double xi:gp)for(double eta:gp){double N[4],dn[4][2]; q4(xi,eta,N,dn); double J00=0,J01=0,J10=0,J11=0; for(int a=0;a<4;a++){J00+=dn[a][0]*xy[a][0];J01+=dn[a][0]*xy[a][1];J10+=dn[a][1]*xy[a][0];J11+=dn[a][1]*xy[a][1];} double det=J00*J11-J01*J10; if(det<=0){cerr<<"S4R detJ error\n";exit(5);} area+=det; double inv00=J11/det,inv01=-J01/det,inv10=-J10/det,inv11=J00/det; vector<double>Bm(3*24,0),Bb(3*24,0); for(int a=0;a<4;a++){double dx=dn[a][0]*inv00+dn[a][1]*inv10,dy=dn[a][0]*inv01+dn[a][1]*inv11; int ix=6*a; Bm[0*24+ix]=dx;Bm[1*24+ix+1]=dy;Bm[2*24+ix]=dy;Bm[2*24+ix+1]=dx; Bb[0*24+ix+4]=dx;Bb[1*24+ix+3]=-dy;Bb[2*24+ix+4]=dy;Bb[2*24+ix+3]=-dx;} matBtDB(Kloc,24,Bm,3,Dm,det); matBtDB(Kloc,24,Bb,3,Db,det);}
    double N[4],dn[4][2]; q4(0,0,N,dn); double J00=0,J01=0,J10=0,J11=0; for(int a=0;a<4;a++){J00+=dn[a][0]*xy[a][0];J01+=dn[a][0]*xy[a][1];J10+=dn[a][1]*xy[a][0];J11+=dn[a][1]*xy[a][1];} double det=J00*J11-J01*J10,inv00=J11/det,inv01=-J01/det,inv10=-J10/det,inv11=J00/det; vector<double>Bs(2*24,0); for(int a=0;a<4;a++){double dx=dn[a][0]*inv00+dn[a][1]*inv10,dy=dn[a][0]*inv01+dn[a][1]*inv11; int ix=6*a; Bs[0*24+ix+2]=dx;Bs[0*24+ix+4]=N[a];Bs[1*24+ix+2]=dy;Bs[1*24+ix+3]=-N[a];} matBtDB(Kloc,24,Bs,2,Ds,det*4); double drill=max(E*t*area*1e-6,1e-6); for(int a=0;a<4;a++) Kloc[(6*a+5)*24+(6*a+5)]+=drill;
    double R[9]={x.x,x.y,x.z,y.x,y.y,y.z,z.x,z.y,z.z}; vector<double>K(24*24,0.0); for(int I=0;I<24;I++)for(int J=0;J<24;J++){double sum=0; for(int a=0;a<24;a++)for(int b=0;b<24;b++) if(I/6==a/6&&J/6==b/6){ if(((I%6)<3)!=((a%6)<3)||((J%6)<3)!=((b%6)<3)) continue; sum+=R[(a%3)*3+(I%3)]*Kloc[a*24+b]*R[(b%3)*3+(J%3)]; } K[I*24+J]=sum;}
    packUpper(K,ND_,Matrix);
}
void CShell4R::ElementStress(double* stress,double*){ for(int i=0;i<12;i++) stress[i]=0.0; }

CTieSpring::CTieSpring()
{
    NEN_ = 2;
    nodes_ = new CNode*[NEN_];
    ND_ = 12;
    LocationMatrix_ = new unsigned int[ND_];
    ElementMaterial_ = nullptr;
    DofCount_ = 3;
}

bool CTieSpring::Read(ifstream& Input, CMaterial* MaterialSets, CNode* NodeList)
{
    unsigned int N1, N2, MSet;
    Input >> N1 >> N2 >> MSet >> DofCount_;
    if (DofCount_ != 3 && DofCount_ != 6)
        DofCount_ = 3;
    ElementMaterial_ = dynamic_cast<CBridgeMaterial*>(MaterialSets) + MSet - 1;
    nodes_[0] = &NodeList[N1 - 1];
    nodes_[1] = &NodeList[N2 - 1];
    return true;
}

void CTieSpring::Write(COutputter& output)
{
    output << setw(11) << nodes_[0]->NodeNumber
           << setw(9) << nodes_[1]->NodeNumber
           << setw(12) << ElementMaterial_->nset
           << setw(8) << DofCount_ << endl;
}

void CTieSpring::GenerateLocationMatrix()
{
    unsigned int k = 0;
    for (unsigned int n = 0; n < NEN_; ++n)
        for (unsigned int d = 0; d < 6; ++d)
            LocationMatrix_[k++] = nodes_[n]->bcode[d];
}

void CTieSpring::ElementStiffness(double* Matrix)
{
    clear(Matrix, SizeOfStiffnessMatrix());
    CBridgeMaterial* mat = dynamic_cast<CBridgeMaterial*>(ElementMaterial_);
    double k = mat->p[0];
    double K[12][12];
    for (int i=0;i<12;i++) for (int j=0;j<12;j++) K[i][j]=0.0;
    for (unsigned int d=0; d<DofCount_; ++d)
    {
        K[d][d] += k;
        K[d+6][d+6] += k;
        K[d][d+6] -= k;
        K[d+6][d] -= k;
    }
    for (unsigned int j=0;j<ND_;++j)
    {
        unsigned int off=(j+1)*j/2;
        for (unsigned int i=0;i<=j;++i)
            Matrix[off+j-i]=K[i][j];
    }
}

void CTieSpring::ElementStress(double* stress, double*)
{
    stress[0] = 0.0;
}
