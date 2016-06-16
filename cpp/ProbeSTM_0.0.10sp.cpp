
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Vec3.cpp"

// ================= CONSTANTS

const double Go = 39.47841760435743; //(2*pi)^2
const double eV = 0.036749034;
const double aB = 1.889725989;
const double rev_dpi = 0.1591549431; // 1/(2 Pi)
const double four_pi = 12.566371; // pi*4
const double N_p = 1.7320508; // Sqrt(3)

// ================= GLOBAL VARIABLES


static double decay = 1;
static double Norm = 1;
Vec3d * dr; //=new Vec3d[NoAt];
double* radial; //= new double[NoAt];
double* rev_rr; //= new double[NoAt];

// ================= INNER FUNCTIONS

 inline double sqr(double x){
 return x*x;
 }

 inline double Lor(double V, double eig, double eta){
 double f = rev_dpi*eta/((V-eig)*(V-eig)+(0.25*eta*eta));
 //printf("lor: %f\n", f);
 return f;
 }

// ================= orbital conductances (for s tip only now)
                          
 inline double ssp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0];					// s  orb. of sample
 f += coe[1]*dR.y*rev_rr*N_p;		// py orb. of sample
 f += coe[2]*dR.z*rev_rr*N_p;		// pz orb. of sample
 f += coe[3]*dR.x*rev_rr*N_p;		// px orb. of sample
 	//printf("%f %f %f %f\n",  coe[0], coe[1]*-1.0*dR.y*rev_rr, coe[2]*dR.z*rev_rr, coe[3]*dR.x*rev_rr);
 return f;
 }

//=================== single orb. conductances ============================
// ******************           s-tip          ****************************
 // ================== sqrt(G) - single mol. orb. (sp) vs. s-tip over atoms
 double Gsatomsp(int NoAt, Vec3d * dr, double* rev_rr, double* radial,  double* coes){
 double f = 0.0;
  for(int iat=0; iat<NoAt; iat++){
	f += radial[iat]*ssp( coes+(4*iat), rev_rr[iat], dr[iat] );
	}
 return f;
 }

// ==================== single point calculations =====================================

 // ================== single point dI/dV calculation s - sp
 double dIdVssp_vec( const Vec3d& r, int NoAt, int NoOrb, double V, double eta,double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    //dri = dri * aB;
	dr[iat]= dri;
	double rri = dri.norm();
	radial[iat] = exp(-(rri*decay));  
	rev_rr[iat] = 1/rri;
	}
 for (int i=0; i<NoOrb; i++){
	f += Lor(V,eig[i],eta)*sqr( Gsatomsp(NoAt,dr,rev_rr,radial, coesin+(i*NoAt*4) ) );
	}
 f *= Go*Norm;
 return f;
 }

// =====================================================
// ==========   Export these functions ( to Python )
// ========================================================

extern "C"{

 // ================== procedure dIdV -s vs. sp sample for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_dIdVssp( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* cur){
 	decay = sqrt( abs(WF+WF)*eV);
	Norm = four_pi*decay;
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * Ratin = (Vec3d *) Rat_;
	dr = new Vec3d[NoAt];
	radial = new double[NoAt];
	rev_rr = new double[NoAt];
 //printf("inside a function sp\n");
	for (int s=0; s<Npoints; s++){
 	 //printf("%d\n", s);
	 cur[s]=dIdVssp_vec( R[s], NoAt, NoOrb, V, eta, eig, Ratin, coesin);
	}
	delete dr; delete radial; delete rev_rr;
 }

}

//******************************************************************************************************************
//******************* !!!!!! PX PART !!!!!!!!! *********************************************************************
//******************************************************************************************************************

// ================= orbital conductances (for px tip only now)

 inline double pxsp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.x*decay;													// s  orb. of sample
 f += coe[1]*N_p*dR.x*dR.y*( decay*rev_rr + sqr(rev_rr) );						// py orb. of sample
 f += coe[2]*N_p*dR.x*dR.z*( decay*rev_rr + sqr(rev_rr) );						// pz orb. of sample
 f += coe[3]*N_p*( -1 + decay*rev_rr*sqr(dR.x) + sqr(rev_rr)*sqr(dR.x) );		// px orb. of sample
 	//printf("%f %f %f %f\n",  coe[0], coe[1]*-1.0*dR.y*rev_rr, coe[2]*dR.z*rev_rr, coe[3]*dR.x*rev_rr);
 return f;
 }

//=================== single orb. conductances ============================
// ******************           px-tip          ****************************
 // ================== sqrt(G) - single mol. orb. (sp) vs. s-tip over atoms
 double Gpxatomsp(int NoAt, Vec3d * dr, double* rev_rr, double* radial,  double* coes){
 //printf("inside Gs \n");
 double f = 0.0;
  for(int iat=0; iat<NoAt; iat++){
	f += radial[iat]*rev_rr[iat]*pxsp( coes+(4*iat), rev_rr[iat], dr[iat] );
 	//printf("%f\n", f);
	}
 //printf("%f\n", f);
 return f;
 }

// ==================== single point calculations =====================================

 // ================== single point dI/dV calculation s - sp
 double dIdVpxsp_vec( const Vec3d& r, int NoAt, int NoOrb, double V, double eta,double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    //dri = dri *aB;
	dr[iat]= dri;
	double rri = dri.norm();
	radial[iat] = exp(-(rri*decay));  
	rev_rr[iat] = 1/rri;
	}
 for (int i=0; i<NoOrb; i++){
 	//printf("inside a loop \n");	
	f += Lor(V,eig[i],eta)*sqr( Gpxatomsp(NoAt,dr,rev_rr,radial, coesin+(i*NoAt*4) ) );
	}
 f *= Go*Norm;
 return f;
 }

// =====================================================
// ==========   Export these functions ( to Python )
// ========================================================

extern "C"{

 // ================== procedure dIdV -s vs. sp sample for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_dIdVpxsp( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* cur){
 	decay = sqrt( abs(WF+WF)*eV);
	Norm = four_pi*decay;
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * Ratin = (Vec3d *) Rat_;
	dr = new Vec3d[NoAt];
	radial = new double[NoAt];
	rev_rr = new double[NoAt];
 //printf("inside a function sp\n");
	for (int s=0; s<Npoints; s++){
 	 //printf("%d\n", s);
	 cur[s]=sqr(N_p/decay)*dIdVpxsp_vec( R[s], NoAt, NoOrb, V, eta, eig, Ratin, coesin);
	}
	delete dr; delete radial; delete rev_rr;
 }

}

//******************************************************************************************************************
//******************* !!!!!! PY PART !!!!!!!!! *********************************************************************
//******************************************************************************************************************

// ================= orbital conductances (for py tip only now)

 inline double pysp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.y*decay;													// s  orb. of sample
 f += coe[1]*N_p*( -1 + decay*rev_rr*sqr(dR.y) + sqr(rev_rr)*sqr(dR.y) );		// py orb. of sample
 f += coe[2]*N_p*dR.y*dR.z*( decay*rev_rr + sqr(rev_rr) );						// pz orb. of sample
 f += coe[3]*N_p*dR.y*dR.x*( decay*rev_rr + sqr(rev_rr) );						// px orb. of sample
 	//printf("%f %f %f %f\n",  coe[0], coe[1]*-1.0*dR.y*rev_rr, coe[2]*dR.z*rev_rr, coe[3]*dR.x*rev_rr);
 return f;
 }

//=================== single orb. conductances ============================
// ******************           py-tip          ****************************
 // ================== sqrt(G) - single mol. orb. (sp) vs. s-tip over atoms
 double Gpyatomsp(int NoAt, Vec3d * dr, double* rev_rr, double* radial,  double* coes){
 //printf("inside Gs \n");
 double f = 0.0;
  for(int iat=0; iat<NoAt; iat++){
	f += radial[iat]*rev_rr[iat]*pysp( coes+(4*iat), rev_rr[iat], dr[iat] );
 	//printf("%f\n", f);
	}
 //printf("%f\n", f);
 return f;
 }

// ==================== single point calculations =====================================

 // ================== single point dI/dV calculation s - sp
 double dIdVpysp_vec( const Vec3d& r, int NoAt, int NoOrb, double V, double eta,double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    //dri = dri *aB;
	dr[iat]= dri;
	double rri = dri.norm();
	radial[iat] = exp(-(rri*decay));  
	rev_rr[iat] = 1/rri;
	}
 for (int i=0; i<NoOrb; i++){
 	//printf("inside a loop \n");	
	f += Lor(V,eig[i],eta)*sqr( Gpyatomsp(NoAt,dr,rev_rr,radial, coesin+(i*NoAt*4) ) );
	}
 f *= Go*Norm;
 return f;
 }

// =====================================================
// ==========   Export these functions ( to Python )
// ========================================================

extern "C"{

 // ================== procedure dIdV -s vs. sp sample for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_dIdVpysp( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* cur){
 	decay = sqrt( abs(WF+WF)*eV);
	Norm = four_pi*decay;
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * Ratin = (Vec3d *) Rat_;
	dr = new Vec3d[NoAt];
	radial = new double[NoAt];
	rev_rr = new double[NoAt];
 //printf("inside a function sp\n");
	for (int s=0; s<Npoints; s++){
 	 //printf("%d\n", s);
	 cur[s]=sqr(N_p/decay)*dIdVpysp_vec( R[s], NoAt, NoOrb, V, eta, eig, Ratin, coesin);
	}
	delete dr; delete radial; delete rev_rr;
 }

}

//******************************************************************************************************************
//******************* !!!!!! PX & PY fast PART !!!!!!!!! *********************************************************************
//******************************************************************************************************************

// ==================== single point calculations =====================================

 // ================== single point dI/dV calculation pxy - sp
 double dIdVpxysp_vec( const Vec3d& r, int NoAt, int NoOrb, double V, double eta,double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    // dri = dri *aB;
	dr[iat]= dri;
	double rri = dri.norm();
	radial[iat] = exp(-(rri*decay));  
	rev_rr[iat] = 1/rri;
	}
 for (int i=0; i<NoOrb; i++){
 	//printf("inside a loop \n");	
	f += Lor(V,eig[i],eta)*(sqr( Gpyatomsp(NoAt,dr,rev_rr,radial, coesin+(i*NoAt*4) ) ) + sqr( Gpxatomsp(NoAt,dr,rev_rr,radial, coesin+(i*NoAt*4) ) ));
	}
 f *= Go*Norm;
 return f;
 }

// =====================================================
// ==========   Export these functions ( to Python )
// ========================================================

extern "C"{

 // ================== procedure dIdV -pxy vs. sp sample for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_dIdVpxysp( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* cur){
 	decay = sqrt( abs(WF+WF)*eV);
	Norm = four_pi*decay;
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * Ratin = (Vec3d *) Rat_;
	dr = new Vec3d[NoAt];
	radial = new double[NoAt];
	rev_rr = new double[NoAt];
 //printf("inside a function sp\n");
	for (int s=0; s<Npoints; s++){
 	 //printf("%d\n", s);
	 cur[s]=sqr(N_p/decay)*dIdVpxysp_vec( R[s], NoAt, NoOrb, V, eta, eig, Ratin, coesin);
	}
	delete dr; delete radial; delete rev_rr;
 }

}

//******************************************************************************************************************
//******************* !!!!!! PZ PART !!!!!!!!! *********************************************************************
//******************************************************************************************************************

// ================= orbital conductances (for py tip only now)

 inline double pzsp(double* coe, double rev_rr, const Vec3d& dR){
 //printf("inside ssp \n");
 double f = coe[0]*dR.z*decay;													// s  orb. of sample
 f += coe[1]*N_p*dR.z*dR.y*( decay*rev_rr + sqr(rev_rr) );						// py orb. of sample
 f += coe[2]*N_p*( -1 + decay*rev_rr*sqr(dR.z) + sqr(rev_rr)*sqr(dR.z) );		// pz orb. of sample
 f += coe[3]*N_p*dR.z*dR.x*( decay*rev_rr + sqr(rev_rr) );						// px orb. of sample
 	//printf("%f %f %f %f\n",  coe[0], coe[1]*-1.0*dR.y*rev_rr, coe[2]*dR.z*rev_rr, coe[3]*dR.x*rev_rr);
 return f;
 }

//=================== single orb. conductances ============================
// ******************           py-tip          ****************************
 // ================== sqrt(G) - single mol. orb. (sp) vs. s-tip over atoms
 double Gpzatomsp(int NoAt, Vec3d * dr, double* rev_rr, double* radial,  double* coes){
 //printf("inside Gs \n");
 double f = 0.0;
  for(int iat=0; iat<NoAt; iat++){
	f += radial[iat]*rev_rr[iat]*pzsp( coes+(4*iat), rev_rr[iat], dr[iat] );
 	//printf("%f\n", f);
	}
 //printf("%f\n", f);
 return f;
 }

// ==================== single point calculations =====================================

 // ================== single point dI/dV calculation s - sp
 double dIdVpzsp_vec( const Vec3d& r, int NoAt, int NoOrb, double V, double eta,double* eig, Vec3d * Ratin, double* coesin){
 //printf("inside a function \n");
 double f = 0.0;
 for(int iat=0; iat<NoAt; iat++){
	//printf("inside first iat pre calc \n");
	Vec3d dri;
	dri.set_sub( r, Ratin[iat] );
	dri.mul(aB);
    //dri = dri *aB;
	dr[iat]= dri;
	double rri = dri.norm();
	radial[iat] = exp(-(rri*decay));  
	rev_rr[iat] = 1/rri;
	}
 for (int i=0; i<NoOrb; i++){
 	//printf("inside a loop \n");	
	f += Lor(V,eig[i],eta)*sqr( Gpzatomsp(NoAt,dr,rev_rr,radial, coesin+(i*NoAt*4) ) );
	}
 f *= Go*Norm;
 //printf("%f\n", decay);
 return f;
 }

// =====================================================
// ==========   Export these functions ( to Python )
// ========================================================

extern "C"{

 // ================== procedure dIdV -s vs. sp sample for a stack of data (basicly 3D of coordinate = 4D) 
 void proc_dIdVpzsp( int NoAt, int NoOrb, int Npoints, double V, double WF, double eta, double* eig, double* R_, double* Rat_, double* coesin, double* cur){
 	decay = sqrt( abs(WF+WF)*eV);
	//printf("%f\n", decay);
	Norm = four_pi*decay;
	Vec3d * R = ( Vec3d * ) R_;
	Vec3d * Ratin = (Vec3d *) Rat_;
	dr = new Vec3d[NoAt];
	radial = new double[NoAt];
	rev_rr = new double[NoAt];
 //printf("inside a function sp\n");
	for (int s=0; s<Npoints; s++){
 	 //printf("%d\n", s);
	 cur[s]=sqr(N_p/decay)*dIdVpzsp_vec( R[s], NoAt, NoOrb, V, eta, eig, Ratin, coesin);
	}
	delete dr; delete radial; delete rev_rr;
 }

}
