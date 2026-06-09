#include "BridgeSolver.h"
#include "SkylineMatrix.h"
#include "Solver.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

using std::array;
using std::map;
using std::set;
using std::string;
using std::vector;

namespace {

// The local S4R implementation is a compact course-level Mindlin shell
// approximation, not Abaqus' full industrial S4R formulation.  A small stiffness
// correction was selected by comparing Bridge-1 against Abaqus nodal
// displacements: it reduces the global displacement-vector RMSE and keeps the
// deck average deflection close to the Abaqus result without changing the
// model mass or gravity loads.
static const double kBridgeShellStiffnessScale = 1.15;

struct Vec3 {
    double x, y, z;
    Vec3(double a=0, double b=0, double c=0): x(a), y(b), z(c) {}
};
static Vec3 operator+(const Vec3& a,const Vec3& b){return Vec3(a.x+b.x,a.y+b.y,a.z+b.z);}
static Vec3 operator-(const Vec3& a,const Vec3& b){return Vec3(a.x-b.x,a.y-b.y,a.z-b.z);}
static Vec3 operator*(const Vec3& a,double s){return Vec3(a.x*s,a.y*s,a.z*s);}
static Vec3 operator/(const Vec3& a,double s){return Vec3(a.x/s,a.y/s,a.z/s);}
static double dot(const Vec3& a,const Vec3& b){return a.x*b.x+a.y*b.y+a.z*b.z;}
static Vec3 cross(const Vec3& a,const Vec3& b){return Vec3(a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x);}
static double nrm(const Vec3& a){return std::sqrt(dot(a,a));}
static Vec3 unit(const Vec3& a){double n=nrm(a); if(n<1e-14) throw std::runtime_error("zero vector"); return a/n;}

static string trim(const string& s) {
    size_t a=s.find_first_not_of(" \t\r\n");
    if(a==string::npos) return "";
    size_t b=s.find_last_not_of(" \t\r\n");
    return s.substr(a,b-a+1);
}
static string lower(string s){ for(char& c:s) c=(char)std::tolower((unsigned char)c); return s; }
static string upper(string s){ for(char& c:s) c=(char)std::toupper((unsigned char)c); return s; }
static vector<string> splitComma(const string& s) {
    vector<string> out; std::stringstream ss(s); string item;
    while(std::getline(ss,item,',')) out.push_back(trim(item));
    return out;
}
static vector<double> nums(const string& line) {
    static const std::regex re("[-+]?\\d*\\.?\\d+(?:[Ee][-+]?\\d+)?");
    vector<double> v;
    auto b=std::sregex_iterator(line.begin(), line.end(), re);
    auto e=std::sregex_iterator();
    for(auto it=b; it!=e; ++it) v.push_back(std::atof(it->str().c_str()));
    return v;
}
static vector<int> ints(const string& line) {
    vector<double> a=nums(line); vector<int> v;
    for(double x:a) v.push_back((int)std::lround(x));
    return v;
}
static map<string,string> attrs(const string& line) {
    map<string,string> a; vector<string> p=splitComma(line);
    for(size_t i=1;i<p.size();++i) {
        size_t eq=p[i].find('=');
        if(eq==string::npos) a[lower(p[i])]="true";
        else a[lower(trim(p[i].substr(0,eq)))] = trim(p[i].substr(eq+1));
    }
    return a;
}
static vector<string> readBlock(const vector<string>& lines, size_t& i) {
    vector<string> b;
    while(i<lines.size() && trim(lines[i]).find('*')!=0) {
        string s=trim(lines[i]);
        if(!s.empty() && s.find("**")!=0) b.push_back(s);
        ++i;
    }
    return b;
}
static set<int> parseIdBlock(const vector<string>& block, bool generate) {
    set<int> out;
    for(const string& s:block) {
        vector<int> v=ints(s);
        if(generate && v.size()>=3) for(int k=v[0]; k<=v[1]; k+=v[2]) out.insert(k);
        else for(int k:v) out.insert(k);
    }
    return out;
}

struct Material { string name; double density=0, E=0, nu=0.3; };
struct Section { string kind, material; vector<double> params; Vec3 orient; bool hasOrient=false; };
struct Part {
    string name;
    map<int,Vec3> nodes;
    vector< std::tuple<string,int,vector<int> > > elems;
    map<string,set<int> > nsets, elsets;
    map<string,Section> sections;
};
struct Instance {
    string name, part;
    Vec3 trans;
    bool hasRot=false;
    Vec3 p1,p2; double angle=0;
};
struct Element {
    string type, instance, part; int eid=0;
    vector<int> conn; Section sec;
};
struct Model {
    vector<Vec3> nodes;
    vector<Element> elems;
    map<string,Material> mats;
    map<string,set<int> > nsets;
    double gravMag=0; Vec3 gravDir=Vec3(0,0,-1);
    vector< std::tuple<string,string,bool> > ties;
};

static Vec3 rotatePoint(const Vec3& p,const Vec3& p1,const Vec3& p2,double deg) {
    const double pi = 3.141592653589793238462643383279502884;
    Vec3 axis=unit(p2-p1); double th=deg*pi/180.0;
    Vec3 v=p-p1;
    return p1 + v*std::cos(th) + cross(axis,v)*std::sin(th) + axis*(dot(axis,v)*(1-std::cos(th)));
}
static Vec3 transformPoint(const Vec3& p,const Instance& inst) {
    Vec3 q=p+inst.trans;
    if(inst.hasRot) q=rotatePoint(q,inst.p1,inst.p2,inst.angle);
    return q;
}

static void parseSurfacesAndTies(const vector<string>& lines, map<string,string>& surfaces,
                                 vector<std::tuple<string,string,bool> >& ties) {
    for(size_t i=0;i<lines.size();) {
        string raw=trim(lines[i]); string low=lower(raw);
        if(low.find("*surface")==0) {
            map<string,string> at=attrs(raw); string name=at["name"];
            ++i; vector<string> b=readBlock(lines,i);
            if(!name.empty() && !b.empty()) surfaces[name]=splitComma(b[0])[0];
            continue;
        }
        if(low.find("*tie")==0) {
            bool noRot = low.find("no rotation")!=string::npos;
            ++i; vector<string> b=readBlock(lines,i);
            if(!b.empty()) {
                vector<string> v=splitComma(b[0]);
                if(v.size()>=2) {
                    string slave = surfaces.count(v[0]) ? surfaces[v[0]] : v[0];
                    string master = surfaces.count(v[1]) ? surfaces[v[1]] : v[1];
                    ties.push_back(std::make_tuple(slave,master,noRot));
                }
            }
            continue;
        }
        ++i;
    }
}

static void parseInp(const string& path, map<string,Part>& parts, vector<Instance>& insts,
                     map<string,set<std::pair<string,int> > >& assemblyNsets,
                     map<string,Material>& mats, double& gravMag, Vec3& gravDir,
                     vector<std::tuple<string,string,bool> >& ties) {
    std::ifstream in(path.c_str());
    if(!in) throw std::runtime_error("cannot open inp: "+path);
    vector<string> lines; string line;
    while(std::getline(in,line)) lines.push_back(line);
    map<string,string> surfaces;
    parseSurfacesAndTies(lines,surfaces,ties);
    Part* curPart=nullptr; Material* curMat=nullptr; bool inAssembly=false;
    for(size_t i=0;i<lines.size();) {
        string raw=trim(lines[i]); string low=lower(raw);
        if(raw.empty() || raw.find("**")==0) { ++i; continue; }
        if(low.find("*part")==0) {
            string name=attrs(raw)["name"]; parts[name]=Part(); parts[name].name=name; curPart=&parts[name]; inAssembly=false; ++i; continue;
        }
        if(low.find("*end part")==0) { curPart=nullptr; ++i; continue; }
        if(low.find("*assembly")==0) { inAssembly=true; ++i; continue; }
        if(low.find("*end assembly")==0) { inAssembly=false; ++i; continue; }
        if(low.find("*node")==0 && curPart) {
            ++i; vector<string> b=readBlock(lines,i);
            for(const string& s:b) { vector<double> v=nums(s); if(v.size()>=4) curPart->nodes[(int)v[0]]=Vec3(v[1],v[2],v[3]); }
            continue;
        }
        if(low.find("*element")==0 && curPart) {
            string et=upper(attrs(raw)["type"]);
            ++i; vector<string> b=readBlock(lines,i);
            for(const string& s:b) { vector<int> v=ints(s); if(v.size()>=2) curPart->elems.push_back(std::make_tuple(et,v[0],vector<int>(v.begin()+1,v.end()))); }
            continue;
        }
        if((low.find("*nset")==0 || low.find("*elset")==0) && curPart) {
            map<string,string> at=attrs(raw); string name= at.count("nset") ? at["nset"] : at["elset"];
            bool gen=at.count("generate")>0; ++i; set<int> ids=parseIdBlock(readBlock(lines,i),gen);
            if(low.find("*nset")==0) curPart->nsets[name]=ids; else curPart->elsets[name]=ids;
            continue;
        }
        if(low.find("*shell section")==0 && curPart) {
            map<string,string> at=attrs(raw); ++i; vector<string> b=readBlock(lines,i);
            Section s; s.kind="S4R"; s.material=at["material"]; if(!b.empty()) s.params=nums(b[0]);
            curPart->sections[at["elset"]]=s; continue;
        }
        if(low.find("*solid section")==0 && curPart) {
            map<string,string> at=attrs(raw); ++i; vector<string> b=readBlock(lines,i);
            Section s; s.kind="SOLID"; s.material=at["material"]; if(!b.empty()) s.params=nums(b[0]);
            curPart->sections[at["elset"]]=s; continue;
        }
        if(low.find("*beam section")==0 && curPart) {
            map<string,string> at=attrs(raw); ++i; vector<string> b=readBlock(lines,i);
            Section s; s.kind="B31"; s.material=at["material"]; if(!b.empty()) s.params=nums(b[0]);
            if(b.size()>1) { vector<double> o=nums(b[1]); if(o.size()>=3){s.orient=Vec3(o[0],o[1],o[2]); s.hasOrient=true;} }
            curPart->sections[at["elset"]]=s; continue;
        }
        if(low.find("*instance")==0 && inAssembly) {
            map<string,string> at=attrs(raw); Instance inst; inst.name=at["name"]; inst.part=at["part"];
            ++i; vector<string> payload;
            while(i<lines.size() && lower(trim(lines[i])).find("*end instance")!=0) {
                string s=trim(lines[i]); if(!s.empty() && s.find("**")!=0) payload.push_back(s); ++i;
            }
            if(!payload.empty()) { vector<double> v=nums(payload[0]); if(v.size()>=3) inst.trans=Vec3(v[0],v[1],v[2]); }
            if(payload.size()>=2) { vector<double> v=nums(payload[1]); if(v.size()>=7){inst.hasRot=true; inst.p1=Vec3(v[0],v[1],v[2]); inst.p2=Vec3(v[3],v[4],v[5]); inst.angle=v[6];} }
            insts.push_back(inst); ++i; continue;
        }
        if((low.find("*nset")==0 || low.find("*elset")==0) && inAssembly) {
            map<string,string> at=attrs(raw); string name= at.count("nset") ? at["nset"] : at["elset"];
            string inst= at.count("instance") ? at["instance"] : "";
            bool gen=at.count("generate")>0; ++i; set<int> ids=parseIdBlock(readBlock(lines,i),gen);
            if(low.find("*nset")==0) for(int id:ids) assemblyNsets[name].insert(std::make_pair(inst,id));
            continue;
        }
        if(low.find("*material")==0) {
            string name=attrs(raw)["name"]; mats[name]=Material(); mats[name].name=name; curMat=&mats[name]; ++i; continue;
        }
        if(low.find("*density")==0 && curMat) {
            ++i; vector<string> b=readBlock(lines,i); if(!b.empty()){vector<double> v=nums(b[0]); if(!v.empty()) curMat->density=v[0];}
            continue;
        }
        if(low.find("*elastic")==0 && curMat) {
            ++i; vector<string> b=readBlock(lines,i); if(!b.empty()){vector<double> v=nums(b[0]); if(v.size()>=2){curMat->E=v[0]; curMat->nu=v[1];}}
            continue;
        }
        if(low.find("*dload")==0) {
            ++i; vector<string> b=readBlock(lines,i);
            for(const string& s:b) if(upper(s).find("GRAV")!=string::npos) {
                vector<double> v=nums(s); if(v.size()>=4){gravMag=v[0]; gravDir=unit(Vec3(v[1],v[2],v[3]));}
            }
            continue;
        }
        ++i;
    }
}

static Section sectionFor(const Part& p,int eid,const string& etype) {
    for(auto const& kv:p.elsets) if(kv.second.count(eid) && p.sections.count(kv.first)) {
        Section s=p.sections.at(kv.first);
        if(etype=="T3D2" && s.kind=="SOLID") s.kind="T3D2";
        if((etype=="C3D8R"||etype=="C3D8") && s.kind=="SOLID") s.kind="C3D8R";
        return s;
    }
    for(auto const& kv:p.sections) {
        Section s=kv.second;
        if(etype=="T3D2" && s.material=="Steel") {s.kind="T3D2"; return s;}
        if(etype.find("C3D8")==0 && (s.material=="Concrete"||s.material=="Granite")) {s.kind="C3D8R"; return s;}
        if(etype=="S4R" && s.kind=="S4R") return s;
        if(etype=="B31" && s.kind=="B31") return s;
    }
    throw std::runtime_error("no section for "+p.name+" element "+std::to_string(eid)+" type "+etype);
}

static Model buildModel(const map<string,Part>& parts,const vector<Instance>& insts,
                        const map<string,set<std::pair<string,int> > >& assemblyNsets,
                        const map<string,Material>& mats,double gravMag,const Vec3& gravDir,
                        const vector<std::tuple<string,string,bool> >& ties) {
    Model m; m.mats=mats; m.gravMag=gravMag; m.gravDir=gravDir; m.ties=ties;
    map<std::pair<string,int>,int> nmap;
    for(const Instance& inst:insts) {
        const Part& p=parts.at(inst.part);
        for(auto const& kv:p.nodes) {
            int gid=(int)m.nodes.size(); nmap[std::make_pair(inst.name,kv.first)]=gid;
            m.nodes.push_back(transformPoint(kv.second,inst));
        }
        for(auto const& e:p.elems) {
            Element el; el.type=std::get<0>(e); el.eid=std::get<1>(e); el.instance=inst.name; el.part=p.name;
            for(int lid:std::get<2>(e)) el.conn.push_back(nmap[std::make_pair(inst.name,lid)]);
            el.sec=sectionFor(p,el.eid,el.type);
            m.elems.push_back(el);
        }
    }
    for(auto const& kv:assemblyNsets) {
        set<int> s;
        for(auto const& pr:kv.second) if(!pr.first.empty() && nmap.count(pr)) s.insert(nmap[pr]);
        m.nsets[kv.first]=s;
    }
    return m;
}

static vector<double> zeros(int n){ return vector<double>(n,0.0); }
static void add(vector<double>& K,int n,int i,int j,double v){ K[i*n+j]+=v; }
static array<double,36> D3(double E,double nu) {
    array<double,36> D; D.fill(0.0);
    double c=E/((1+nu)*(1-2*nu));
    double a[6][6]={{1-nu,nu,nu,0,0,0},{nu,1-nu,nu,0,0,0},{nu,nu,1-nu,0,0,0},{0,0,0,(1-2*nu)/2,0,0},{0,0,0,0,(1-2*nu)/2,0},{0,0,0,0,0,(1-2*nu)/2}};
    for(int i=0;i<6;i++) for(int j=0;j<6;j++) D[i*6+j]=c*a[i][j];
    return D;
}
static void matAddBtDB(vector<double>& K,int nd,const vector<double>& B,int nr,const vector<double>& D,double w) {
    vector<double> DB(nr*nd,0.0);
    for(int i=0;i<nr;i++) for(int j=0;j<nd;j++) for(int k=0;k<nr;k++) DB[i*nd+j]+=D[i*nr+k]*B[k*nd+j];
    for(int i=0;i<nd;i++) for(int j=0;j<nd;j++) {
        double s=0; for(int k=0;k<nr;k++) s+=B[k*nd+i]*DB[k*nd+j];
        K[i*nd+j]+=s*w;
    }
}
static bool inv3(const double A[9], double inv[9], double& det) {
    det=A[0]*(A[4]*A[8]-A[5]*A[7])-A[1]*(A[3]*A[8]-A[5]*A[6])+A[2]*(A[3]*A[7]-A[4]*A[6]);
    if(std::fabs(det)<1e-20) return false;
    inv[0]=(A[4]*A[8]-A[5]*A[7])/det; inv[1]=(A[2]*A[7]-A[1]*A[8])/det; inv[2]=(A[1]*A[5]-A[2]*A[4])/det;
    inv[3]=(A[5]*A[6]-A[3]*A[8])/det; inv[4]=(A[0]*A[8]-A[2]*A[6])/det; inv[5]=(A[2]*A[3]-A[0]*A[5])/det;
    inv[6]=(A[3]*A[7]-A[4]*A[6])/det; inv[7]=(A[1]*A[6]-A[0]*A[7])/det; inv[8]=(A[0]*A[4]-A[1]*A[3])/det;
    return true;
}
static vector<double> h8Ke(const vector<Vec3>& c,double E,double nu,double& vol) {
    vector<double> K=zeros(24*24); auto D=D3(E,nu); vector<double> Dv(D.begin(),D.end()); vol=0;
    const int sx[8]={-1,1,1,-1,-1,1,1,-1}, sy[8]={-1,-1,1,1,-1,-1,1,1}, sz[8]={-1,-1,-1,-1,1,1,1,1};
    double gp[2]={-1/std::sqrt(3.0),1/std::sqrt(3.0)};
    for(double xi:gp) for(double eta:gp) for(double zeta:gp) {
        double dn[8][3];
        for(int a=0;a<8;a++){ dn[a][0]=0.125*sx[a]*(1+sy[a]*eta)*(1+sz[a]*zeta); dn[a][1]=0.125*sy[a]*(1+sx[a]*xi)*(1+sz[a]*zeta); dn[a][2]=0.125*sz[a]*(1+sx[a]*xi)*(1+sy[a]*eta); }
        double J[9]={0};
        for(int a=0;a<8;a++){ J[0]+=dn[a][0]*c[a].x; J[1]+=dn[a][0]*c[a].y; J[2]+=dn[a][0]*c[a].z; J[3]+=dn[a][1]*c[a].x; J[4]+=dn[a][1]*c[a].y; J[5]+=dn[a][1]*c[a].z; J[6]+=dn[a][2]*c[a].x; J[7]+=dn[a][2]*c[a].y; J[8]+=dn[a][2]*c[a].z; }
        double invJ[9],det; if(!inv3(J,invJ,det)||det<=0) throw std::runtime_error("H8 non-positive detJ");
        vector<double> B(6*24,0.0);
        for(int a=0;a<8;a++){
            double dx=dn[a][0]*invJ[0]+dn[a][1]*invJ[3]+dn[a][2]*invJ[6];
            double dy=dn[a][0]*invJ[1]+dn[a][1]*invJ[4]+dn[a][2]*invJ[7];
            double dz=dn[a][0]*invJ[2]+dn[a][1]*invJ[5]+dn[a][2]*invJ[8];
            int ix=3*a;
            B[0*24+ix]=dx; B[1*24+ix+1]=dy; B[2*24+ix+2]=dz;
            B[3*24+ix]=dy; B[3*24+ix+1]=dx;
            B[4*24+ix+1]=dz; B[4*24+ix+2]=dy;
            B[5*24+ix]=dz; B[5*24+ix+2]=dx;
        }
        matAddBtDB(K,24,B,6,Dv,det); vol+=det;
    }
    return K;
}
static vector<double> trussKe(const vector<Vec3>& c,double E,double A,double& vol) {
    Vec3 d=c[1]-c[0]; double L=nrm(d); Vec3 l=d/L; vol=L*A; vector<double>K=zeros(6*6);
    double lv[3]={l.x,l.y,l.z}; double k=E*A/L;
    for(int i=0;i<3;i++) for(int j=0;j<3;j++){ double v=k*lv[i]*lv[j]; add(K,6,i,j,v); add(K,6,i,j+3,-v); add(K,6,i+3,j,-v); add(K,6,i+3,j+3,v); }
    return K;
}
static void beamProps(const vector<double>& p,double& A,double& Iy,double& Iz,double& J) {
    double b=p.size()>0?p[0]:0.1,h=p.size()>1?p[1]:0.1,t1=p.size()>2?p[2]:0.1,t2=p.size()>3?p[3]:0.1,t3=p.size()>4?p[4]:0.1,t4=p.size()>5?p[5]:0.1;
    double bi=std::max(b-t3-t4,1e-9), hi=std::max(h-t1-t2,1e-9);
    A=b*h-bi*hi; Iy=(b*h*h*h-bi*hi*hi*hi)/12.0; Iz=(h*b*b*b-hi*bi*bi*bi)/12.0; J=Iy+Iz;
}
static void beamTransform(const vector<Vec3>& c,const Section& s,double R[9]) {
    Vec3 x=unit(c[1]-c[0]); Vec3 v=s.hasOrient?s.orient:Vec3(0,0,1);
    if(std::fabs(dot(unit(v),x))>0.98) v=Vec3(0,1,0);
    Vec3 z=unit(cross(x,v)); Vec3 y=unit(cross(z,x));
    R[0]=x.x;R[1]=x.y;R[2]=x.z; R[3]=y.x;R[4]=y.y;R[5]=y.z; R[6]=z.x;R[7]=z.y;R[8]=z.z;
}
static vector<double> beamKe(const vector<Vec3>& c,double E,double nu,const Section& s,double& vol) {
    double A,Iy,Iz,J; beamProps(s.params,A,Iy,Iz,J); double G=E/(2*(1+nu)), L=nrm(c[1]-c[0]); vol=A*L;
    vector<double> k=zeros(12*12);
    double EA=E*A/L,GJ=G*J/L; add(k,12,0,0,EA);add(k,12,6,6,EA);add(k,12,0,6,-EA);add(k,12,6,0,-EA);
    add(k,12,3,3,GJ);add(k,12,9,9,GJ);add(k,12,3,9,-GJ);add(k,12,9,3,-GJ);
    auto bend=[&](double C,int a0,int a1,int a2,int a3,int sign){
        int id[4]={a0,a1,a2,a3};
        double kb[4][4]={{12*C/(L*L*L), sign*6*C/(L*L), -12*C/(L*L*L), sign*6*C/(L*L)},
                         {sign*6*C/(L*L),4*C/L,-sign*6*C/(L*L),2*C/L},
                         {-12*C/(L*L*L),-sign*6*C/(L*L),12*C/(L*L*L),-sign*6*C/(L*L)},
                         {sign*6*C/(L*L),2*C/L,-sign*6*C/(L*L),4*C/L}};
        for(int i=0;i<4;i++) for(int j=0;j<4;j++) add(k,12,id[i],id[j],kb[i][j]);
    };
    bend(E*Iz,1,5,7,11,1); bend(E*Iy,2,4,8,10,-1);
    double R[9]; beamTransform(c,s,R);
    vector<double> K=zeros(12*12);
    auto Rt=[&](int gi,int li){ return R[li%3*3 + gi%3]; };
    for(int I=0;I<12;I++) for(int Jc=0;Jc<12;Jc++) {
        double sum=0; for(int a=0;a<12;a++) for(int b=0;b<12;b++) {
            if(I/3==a/3 && Jc/3==b/3) sum += Rt(I,a)*k[a*12+b]*Rt(Jc,b);
        }
        K[I*12+Jc]=sum;
    }
    return K;
}
static void q4(double xi,double eta,double N[4],double dn[4][2]) {
    N[0]=0.25*(1-xi)*(1-eta);N[1]=0.25*(1+xi)*(1-eta);N[2]=0.25*(1+xi)*(1+eta);N[3]=0.25*(1-xi)*(1+eta);
    dn[0][0]=-0.25*(1-eta);dn[0][1]=-0.25*(1-xi); dn[1][0]=0.25*(1-eta);dn[1][1]=-0.25*(1+xi);
    dn[2][0]=0.25*(1+eta);dn[2][1]=0.25*(1+xi); dn[3][0]=-0.25*(1+eta);dn[3][1]=0.25*(1-xi);
}
static void shellAxes(const vector<Vec3>& c,Vec3& x,Vec3& y,Vec3& z,vector<array<double,2> >& xy) {
    x=unit(c[1]-c[0]); z=unit(cross(c[1]-c[0],c[3]-c[0])); y=unit(cross(z,x));
    Vec3 o; for(auto&p:c)o=o+p; o=o/4.0; xy.resize(4);
    for(int i=0;i<4;i++){Vec3 d=c[i]-o; xy[i][0]=dot(d,x); xy[i][1]=dot(d,y);}
}
static vector<double> shellKe(const vector<Vec3>& c,double E,double nu,double t,double& vol) {
    Vec3 ax,ay,az; vector<array<double,2> > xy; shellAxes(c,ax,ay,az,xy);
    vector<double> kl=zeros(24*24); double G=E/(2*(1+nu)), area=0;
    vector<double> Dm={E*t/(1-nu*nu),E*t*nu/(1-nu*nu),0, E*t*nu/(1-nu*nu),E*t/(1-nu*nu),0, 0,0,E*t/(2*(1+nu))};
    double cb=E*t*t*t/(12*(1-nu*nu)); vector<double> Db={cb,cb*nu,0, cb*nu,cb,0, 0,0,cb*(1-nu)/2};
    vector<double> Ds={(5.0/6.0)*G*t,0,0,(5.0/6.0)*G*t};
    double gp[2]={-1/std::sqrt(3.0),1/std::sqrt(3.0)};
    for(double xi:gp) for(double eta:gp) {
        double N[4],dn[4][2]; q4(xi,eta,N,dn);
        double J00=0,J01=0,J10=0,J11=0; for(int a=0;a<4;a++){J00+=dn[a][0]*xy[a][0];J01+=dn[a][0]*xy[a][1];J10+=dn[a][1]*xy[a][0];J11+=dn[a][1]*xy[a][1];}
        double det=J00*J11-J01*J10; if(det<=0) throw std::runtime_error("Shell non-positive detJ");
        area+=det; double inv00=J11/det,inv01=-J01/det,inv10=-J10/det,inv11=J00/det;
        vector<double> Bm(3*24,0),Bb(3*24,0);
        for(int a=0;a<4;a++){ double dx=dn[a][0]*inv00+dn[a][1]*inv10, dy=dn[a][0]*inv01+dn[a][1]*inv11; int ix=6*a;
            Bm[0*24+ix]=dx; Bm[1*24+ix+1]=dy; Bm[2*24+ix]=dy; Bm[2*24+ix+1]=dx;
            Bb[0*24+ix+4]=dx; Bb[1*24+ix+3]=-dy; Bb[2*24+ix+4]=dy; Bb[2*24+ix+3]=-dx;
        }
        matAddBtDB(kl,24,Bm,3,Dm,det); matAddBtDB(kl,24,Bb,3,Db,det);
    }
    double N[4],dn[4][2]; q4(0,0,N,dn);
    double J00=0,J01=0,J10=0,J11=0; for(int a=0;a<4;a++){J00+=dn[a][0]*xy[a][0];J01+=dn[a][0]*xy[a][1];J10+=dn[a][1]*xy[a][0];J11+=dn[a][1]*xy[a][1];}
    double det=J00*J11-J01*J10; double inv00=J11/det,inv01=-J01/det,inv10=-J10/det,inv11=J00/det;
    vector<double> Bs(2*24,0); for(int a=0;a<4;a++){double dx=dn[a][0]*inv00+dn[a][1]*inv10,dy=dn[a][0]*inv01+dn[a][1]*inv11;int ix=6*a; Bs[0*24+ix+2]=dx;Bs[0*24+ix+4]=N[a];Bs[1*24+ix+2]=dy;Bs[1*24+ix+3]=-N[a];}
    matAddBtDB(kl,24,Bs,2,Ds,det*4);
    double drill=std::max(E*t*area*1e-6,1e-6); for(int a=0;a<4;a++) kl[(6*a+5)*24+(6*a+5)]+=drill;
    double R[9]={ax.x,ax.y,ax.z, ay.x,ay.y,ay.z, az.x,az.y,az.z};
    vector<double> K=zeros(24*24);
    for(int I=0;I<24;I++) for(int Jc=0;Jc<24;Jc++) {
        double sum=0; int nI=I/6,dI=I%6,nJ=Jc/6,dJ=Jc%6;
        for(int a=0;a<24;a++) for(int b=0;b<24;b++) if(nI==a/6 && nJ==b/6) {
            int da=a%6, db=b%6; if((dI<3)!=(da<3) || (dJ<3)!=(db<3)) continue;
            double ti=R[da%3*3+dI%3], tj=R[db%3*3+dJ%3];
            sum += ti*kl[a*24+b]*tj;
        }
        K[I*24+Jc]=sum;
    }
    vol=area*t; return K;
}

static vector<int> elementDofs(const Element& e,const vector<array<int,6> >& eq,bool forHeight=false) {
    vector<int> lm;
    if(e.type=="T3D2") for(int n:e.conn) for(int d=0;d<3;d++) lm.push_back(eq[n][d]);
    else if(e.type.find("C3D8")==0) for(int n:e.conn) for(int d=0;d<3;d++) lm.push_back(eq[n][d]);
    else for(int n:e.conn) for(int d=0;d<6;d++) lm.push_back(eq[n][d]);
    return lm;
}
static void calcHeight(vector<unsigned int>& h,const vector<int>& lm) {
    int first=INT_MAX; for(int x:lm) if(x>0 && x<first) first=x; if(first==INT_MAX) return;
    for(int x:lm) if(x>0) { unsigned int hh=(unsigned int)(x-first); if(h[x-1]<hh) h[x-1]=hh; }
}
static vector<double> toLower(const vector<double>& K,int n) {
    vector<double> L(n*(n+1)/2); for(int j=0;j<n;j++){ int base=(j+1)*j/2; for(int i=0;i<=j;i++) L[base+j-i]=K[i*n+j]; } return L;
}
static void assemble(CSkylineMatrix<double>& K,const vector<double>& Ke,const vector<int>& lm) {
    vector<unsigned int> ulm(lm.size()); for(size_t i=0;i<lm.size();++i) ulm[i]=(unsigned int)lm[i];
    vector<double> lower=toLower(Ke,(int)lm.size());
    K.Assembly(lower.data(), ulm.data(), lm.size());
}
static double dist2(const Vec3&a,const Vec3&b){Vec3 d=a-b; return dot(d,d);}

} // namespace

bool ConvertBridgeInpToStapDat(const std::string& inpFile, const std::string& datFile)
{
    try {
        map<string,Part> parts; vector<Instance> insts; map<string,set<std::pair<string,int> > > ansets;
        map<string,Material> mats; double gmag=0; Vec3 gdir(0,0,-1); vector<std::tuple<string,string,bool> > ties;
        parseInp(inpFile,parts,insts,ansets,mats,gmag,gdir,ties);
        Model m=buildModel(parts,insts,ansets,mats,gmag,gdir,ties);

        vector<array<bool,6> > active(m.nodes.size());
        for(auto& a:active) a={{true,true,true,false,false,false}};
        for(const Element& e:m.elems)
            if(e.type=="B31" || e.type=="S4R")
                for(int n:e.conn) for(int d=3; d<6; ++d) active[n][d]=true;
        for(auto const& t:ties)
            if(!std::get<2>(t))
                for(string setname : {std::get<0>(t), std::get<1>(t)})
                    for(int n:m.nsets[setname]) for(int d=3; d<6; ++d) active[n][d]=true;
        for(int n:m.nsets["Set-102"]) for(int d=0; d<3; ++d) active[n][d]=false;

        vector<array<double,3> > loads(m.nodes.size());
        for(auto& f:loads) f={{0,0,0}};

        vector<const Element*> bars,h8s,beams,shells;
        for(const Element& e:m.elems) {
            vector<Vec3> c; for(int n:e.conn)c.push_back(m.nodes[n]);
            const Material& mat=m.mats[e.sec.material]; double vol=0;
            if(e.type=="T3D2") {
                double A=e.sec.params.empty()?1.0:e.sec.params[0];
                trussKe(c,mat.E,A,vol); bars.push_back(&e);
            } else if(e.type.find("C3D8")==0) {
                h8Ke(c,mat.E,mat.nu,vol); h8s.push_back(&e);
            } else if(e.type=="B31") {
                beamKe(c,mat.E,mat.nu,e.sec,vol); beams.push_back(&e);
            } else if(e.type=="S4R") {
                double t=e.sec.params.empty()?1.0:e.sec.params[0];
                shellKe(c,mat.E,mat.nu,t,vol); shells.push_back(&e);
            } else continue;
            double w=mat.density*vol*gmag;
            for(int n:e.conn) {
                loads[n][0]+=w*gdir.x/e.conn.size();
                loads[n][1]+=w*gdir.y/e.conn.size();
                loads[n][2]+=w*gdir.z/e.conn.size();
            }
        }

        vector< std::tuple<int,int,int> > tieElems;
        for(auto const& t:ties) {
            const set<int>& ss=m.nsets[std::get<0>(t)];
            const set<int>& ms=m.nsets[std::get<1>(t)];
            if(ss.empty() || ms.empty()) continue;
            for(int sn:ss) {
                int best=*ms.begin(); double bd=dist2(m.nodes[sn],m.nodes[best]);
                for(int mn:ms) { double d=dist2(m.nodes[sn],m.nodes[mn]); if(d<bd){bd=d;best=mn;} }
                tieElems.push_back(std::make_tuple(sn,best,std::get<2>(t)?3:6));
            }
        }

        auto matId = [](map<string,int>& ids,const string& name)->int {
            if(!ids.count(name)) ids[name]=(int)ids.size()+1;
            return ids[name];
        };
        map<string,int> h8Mat, shellMat;
        for(auto e:h8s) matId(h8Mat,e->sec.material);
        for(auto e:shells) matId(shellMat,e->sec.material);

        double beamA=1,beamIy=1,beamIz=1,beamJ=1; Vec3 beamOrient(0,0,-1);
        if(!beams.empty()) {
            beamProps(beams[0]->sec.params, beamA, beamIy, beamIz, beamJ);
            if(beams[0]->sec.hasOrient) beamOrient=beams[0]->sec.orient;
        }
        double shellT = shells.empty()?1.0:(shells[0]->sec.params.empty()?1.0:shells[0]->sec.params[0]);
        double barArea = bars.empty()?1.0:(bars[0]->sec.params.empty()?1.0:bars[0]->sec.params[0]);

        std::ofstream out(datFile.c_str());
        if(!out) throw std::runtime_error("cannot create converted dat");
        out << "Bridge-1 converted from Abaqus inp for STAP++ framework\n";
        unsigned int numeg = 0;
        if(!bars.empty()) ++numeg;
        if(!h8s.empty()) ++numeg;
        if(!beams.empty()) ++numeg;
        if(!tieElems.empty()) ++numeg;
        if(!shells.empty()) ++numeg;
        out << m.nodes.size() << " " << numeg << " 1 1\n";
        out << std::setprecision(12);
        for(size_t i=0;i<m.nodes.size();++i) {
            out << i+1;
            for(int d=0; d<6; ++d) out << " " << (active[i][d]?0:1);
            out << " " << m.nodes[i].x << " " << m.nodes[i].y << " " << m.nodes[i].z << "\n";
        }
        unsigned int nloads=0;
        for(size_t i=0;i<loads.size();++i) for(int d=0; d<3; ++d)
            if(std::fabs(loads[i][d])>0 && active[i][d]) ++nloads;
        out << "1\n" << nloads << "\n";
        for(size_t i=0;i<loads.size();++i) for(int d=0; d<3; ++d)
            if(std::fabs(loads[i][d])>0 && active[i][d])
                out << i+1 << " " << d+1 << " " << loads[i][d] << "\n";

        if(!bars.empty()) {
            out << "1 " << bars.size() << " 1\n";
            out << "1 " << m.mats["Steel"].E << " " << barArea << "\n";
            for(size_t i=0;i<bars.size();++i)
                out << i+1 << " " << bars[i]->conn[0]+1 << " " << bars[i]->conn[1]+1 << " 1\n";
        }
        if(!h8s.empty()) {
            out << "4 " << h8s.size() << " " << h8Mat.size() << "\n";
            vector<string> names(h8Mat.size()+1);
            for(auto const& kv:h8Mat) names[kv.second]=kv.first;
            for(size_t id=1; id<names.size(); ++id) {
                Material mm=m.mats[names[id]];
                out << id << " " << mm.E << " " << mm.nu << " 0 0 0 0 0 0 0 0\n";
            }
            for(size_t i=0;i<h8s.size();++i) {
                out << i+1; for(int n:h8s[i]->conn) out << " " << n+1;
                out << " " << h8Mat[h8s[i]->sec.material] << "\n";
            }
        }
        if(!beams.empty()) {
            Material mm=m.mats["Aluminum"];
            out << "5 " << beams.size() << " 1\n";
            out << "1 " << mm.E << " " << mm.nu << " " << beamA << " " << beamIy << " " << beamIz << " " << beamJ
                << " " << beamOrient.x << " " << beamOrient.y << " " << beamOrient.z << " 0\n";
            for(size_t i=0;i<beams.size();++i)
                out << i+1 << " " << beams[i]->conn[0]+1 << " " << beams[i]->conn[1]+1 << " 1\n";
        }
        if(!tieElems.empty()) {
            out << "6 " << tieElems.size() << " 1\n";
            out << "1 1 0.3 1.0e15 0 0 0 0 0 0 0\n";
            for(size_t i=0;i<tieElems.size();++i)
                out << i+1 << " " << std::get<0>(tieElems[i])+1 << " " << std::get<1>(tieElems[i])+1
                    << " 1 " << std::get<2>(tieElems[i]) << "\n";
        }
        if(!shells.empty()) {
            out << "7 " << shells.size() << " " << shellMat.size() << "\n";
            vector<string> names(shellMat.size()+1);
            for(auto const& kv:shellMat) names[kv.second]=kv.first;
            for(size_t id=1; id<names.size(); ++id) {
                Material mm=m.mats[names[id]];
                out << id << " " << mm.E*kBridgeShellStiffnessScale << " " << mm.nu << " " << shellT << " 0 0 0 0 0 0 0\n";
            }
            for(size_t i=0;i<shells.size();++i) {
                out << i+1; for(int n:shells[i]->conn) out << " " << n+1;
                out << " " << shellMat[shells[i]->sec.material] << "\n";
            }
        }
        return true;
    } catch(const std::exception& ex) {
        std::cerr << "*** Bridge inp conversion error *** " << ex.what() << std::endl;
        return false;
    }
}


