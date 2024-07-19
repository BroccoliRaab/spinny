#ifndef PGA3D_H
#define PGA3D_H

// 3D Projective Geometric Algebra
// Written by a generator written by enki.
// Modified for C by Broccoliraab
#include <stdio.h>
#include <math.h>

/*
#define PI 3.14159265358979323846

static const char* basis[] = { "1","e0","e1","e2","e3","e01","e02","e03","e12","e31","e23","e021","e013","e032","e123","e0123" }

class PGA3D {
  public:
    PGA3D ()  { std::fill( mvec, mvec + sizeof( mvec )/4, 0.0f ); }
    PGA3D (float f, int idx=0) { std::fill( mvec, mvec + sizeof( mvec )/4, 0.0f ); mvec[idx] = f; }
    float& operator [] (size_t idx) { return mvec[idx]; }
    const float& operator [] (size_t idx) const { return mvec[idx]; }
    PGA3D log () { int n=0; for (int i=0,j=0;i<16;i++) if (mvec[i]!=0.0f) { n++; printf("%s%0.7g%s",(j>0)?" + ":"",mvec[i],(i==0)?"":basis[i]); j++; }if (n==0) printf("0");  printf("\n"); return *this; }
    PGA3D Conjugate(); 
    PGA3D Involute();
    float norm();
    float inorm();
    PGA3D normalized();
  private:  
    float mvec[16];
}
*/
typedef struct multivector_s multivector_t;
struct multivector_s
{
    float mvec[16];
};
//***********************
// PGA3D.Reverse : res = ~a
// Reverse the order of the basis blades.
//***********************
static inline multivector_t PGA3D_reverse (const multivector_t mv) {
  multivector_t r;
  const float * a = &mv.mvec[0];
  float * res = &r.mvec[0];
  res[0]=a[0];
  res[1]=a[1];
  res[2]=a[2];
  res[3]=a[3];
  res[4]=a[4];
  res[5]=-a[5];
  res[6]=-a[6];
  res[7]=-a[7];
  res[8]=-a[8];
  res[9]=-a[9];
  res[10]=-a[10];
  res[11]=-a[11];
  res[12]=-a[12];
  res[13]=-a[13];
  res[14]=-a[14];
  res[15]=a[15];
  return r;
}

//***********************
// PGA3D.Dual : res = !a
// Poincare duality operator.
//***********************
static inline multivector_t PGA3D_dual (const multivector_t mv) {
  multivector_t r;
  const float * a = &mv.mvec[0];
  float * res = &r.mvec[0];
  res[0]=a[15];
  res[1]=a[14];
  res[2]=a[13];
  res[3]=a[12];
  res[4]=a[11];
  res[5]=a[10];
  res[6]=a[9];
  res[7]=a[8];
  res[8]=a[7];
  res[9]=a[6];
  res[10]=a[5];
  res[11]=a[4];
  res[12]=a[3];
  res[13]=a[2];
  res[14]=a[1];
  res[15]=a[0];
  return r;
}

//***********************
// PGA3D.Conjugate : res = a.Conjugate()
// Clifford Conjugation
//***********************
static inline multivector_t PGA3D_conjugate (const multivector_t mv) {
  multivector_t r;
  const float * a = &mv.mvec[0];
  float * res = &r.mvec[0];
  res[0]=mv.mvec[0];
  res[1]=-mv.mvec[1];
  res[2]=-mv.mvec[2];
  res[3]=-mv.mvec[3];
  res[4]=-mv.mvec[4];
  res[5]=-mv.mvec[5];
  res[6]=-mv.mvec[6];
  res[7]=-mv.mvec[7];
  res[8]=-mv.mvec[8];
  res[9]=-mv.mvec[9];
  res[10]=-mv.mvec[10];
  res[11]=mv.mvec[11];
  res[12]=mv.mvec[12];
  res[13]=mv.mvec[13];
  res[14]=mv.mvec[14];
  res[15]=mv.mvec[15];
  return r;
}

//***********************
// PGA3D.Involute : res = a.Involute()
// Main involution
//***********************
static inline multivector_t PGA3D_involute (const multivector_t mv) {
  multivector_t r;
  const float * a = &mv.mvec[0];
  float * res = &r.mvec[0];
  res[0]=mv.mvec[0];
  res[1]=-mv.mvec[1];
  res[2]=-mv.mvec[2];
  res[3]=-mv.mvec[3];
  res[4]=-mv.mvec[4];
  res[5]=mv.mvec[5];
  res[6]=mv.mvec[6];
  res[7]=mv.mvec[7];
  res[8]=mv.mvec[8];
  res[9]=mv.mvec[9];
  res[10]=mv.mvec[10];
  res[11]=-mv.mvec[11];
  res[12]=-mv.mvec[12];
  res[13]=-mv.mvec[13];
  res[14]=-mv.mvec[14];
  res[15]=mv.mvec[15];
  return r;
}

//***********************
// PGA3D.Mul : res = a * b 
// The geometric product.
//***********************
static inline multivector_t PGA3D_mul (const multivector_t amv,const multivector_t bmv) {
  multivector_t r;
  const float * a = &amv.mvec[0];
  const float * b = &bmv.mvec[0];
  float * res = &r.mvec[0];
  res[0]=b[0]*a[0]+b[2]*a[2]+b[3]*a[3]+b[4]*a[4]-b[8]*a[8]-b[9]*a[9]-b[10]*a[10]-b[14]*a[14];
  res[1]=b[1]*a[0]+b[0]*a[1]-b[5]*a[2]-b[6]*a[3]-b[7]*a[4]+b[2]*a[5]+b[3]*a[6]+b[4]*a[7]+b[11]*a[8]+b[12]*a[9]+b[13]*a[10]+b[8]*a[11]+b[9]*a[12]+b[10]*a[13]+b[15]*a[14]-b[14]*a[15];
  res[2]=b[2]*a[0]+b[0]*a[2]-b[8]*a[3]+b[9]*a[4]+b[3]*a[8]-b[4]*a[9]-b[14]*a[10]-b[10]*a[14];
  res[3]=b[3]*a[0]+b[8]*a[2]+b[0]*a[3]-b[10]*a[4]-b[2]*a[8]-b[14]*a[9]+b[4]*a[10]-b[9]*a[14];
  res[4]=b[4]*a[0]-b[9]*a[2]+b[10]*a[3]+b[0]*a[4]-b[14]*a[8]+b[2]*a[9]-b[3]*a[10]-b[8]*a[14];
  res[5]=b[5]*a[0]+b[2]*a[1]-b[1]*a[2]-b[11]*a[3]+b[12]*a[4]+b[0]*a[5]-b[8]*a[6]+b[9]*a[7]+b[6]*a[8]-b[7]*a[9]-b[15]*a[10]-b[3]*a[11]+b[4]*a[12]+b[14]*a[13]-b[13]*a[14]-b[10]*a[15];
  res[6]=b[6]*a[0]+b[3]*a[1]+b[11]*a[2]-b[1]*a[3]-b[13]*a[4]+b[8]*a[5]+b[0]*a[6]-b[10]*a[7]-b[5]*a[8]-b[15]*a[9]+b[7]*a[10]+b[2]*a[11]+b[14]*a[12]-b[4]*a[13]-b[12]*a[14]-b[9]*a[15];
  res[7]=b[7]*a[0]+b[4]*a[1]-b[12]*a[2]+b[13]*a[3]-b[1]*a[4]-b[9]*a[5]+b[10]*a[6]+b[0]*a[7]-b[15]*a[8]+b[5]*a[9]-b[6]*a[10]+b[14]*a[11]-b[2]*a[12]+b[3]*a[13]-b[11]*a[14]-b[8]*a[15];
  res[8]=b[8]*a[0]+b[3]*a[2]-b[2]*a[3]+b[14]*a[4]+b[0]*a[8]+b[10]*a[9]-b[9]*a[10]+b[4]*a[14];
  res[9]=b[9]*a[0]-b[4]*a[2]+b[14]*a[3]+b[2]*a[4]-b[10]*a[8]+b[0]*a[9]+b[8]*a[10]+b[3]*a[14];
  res[10]=b[10]*a[0]+b[14]*a[2]+b[4]*a[3]-b[3]*a[4]+b[9]*a[8]-b[8]*a[9]+b[0]*a[10]+b[2]*a[14];
  res[11]=b[11]*a[0]-b[8]*a[1]+b[6]*a[2]-b[5]*a[3]+b[15]*a[4]-b[3]*a[5]+b[2]*a[6]-b[14]*a[7]-b[1]*a[8]+b[13]*a[9]-b[12]*a[10]+b[0]*a[11]+b[10]*a[12]-b[9]*a[13]+b[7]*a[14]-b[4]*a[15];
  res[12]=b[12]*a[0]-b[9]*a[1]-b[7]*a[2]+b[15]*a[3]+b[5]*a[4]+b[4]*a[5]-b[14]*a[6]-b[2]*a[7]-b[13]*a[8]-b[1]*a[9]+b[11]*a[10]-b[10]*a[11]+b[0]*a[12]+b[8]*a[13]+b[6]*a[14]-b[3]*a[15];
  res[13]=b[13]*a[0]-b[10]*a[1]+b[15]*a[2]+b[7]*a[3]-b[6]*a[4]-b[14]*a[5]-b[4]*a[6]+b[3]*a[7]+b[12]*a[8]-b[11]*a[9]-b[1]*a[10]+b[9]*a[11]-b[8]*a[12]+b[0]*a[13]+b[5]*a[14]-b[2]*a[15];
  res[14]=b[14]*a[0]+b[10]*a[2]+b[9]*a[3]+b[8]*a[4]+b[4]*a[8]+b[3]*a[9]+b[2]*a[10]+b[0]*a[14];
  res[15]=b[15]*a[0]+b[14]*a[1]+b[13]*a[2]+b[12]*a[3]+b[11]*a[4]+b[10]*a[5]+b[9]*a[6]+b[8]*a[7]+b[7]*a[8]+b[6]*a[9]+b[5]*a[10]-b[4]*a[11]-b[3]*a[12]-b[2]*a[13]-b[1]*a[14]+b[0]*a[15];
  return r;
}

//***********************
// PGA3D.Wedge : res = a ^ b 
// The outer product. (MEET)
//***********************
static inline multivector_t PGA3D_wedge (const multivector_t amv,const multivector_t bmv) {
  multivector_t r;
  const float * a = &amv.mvec[0];
  const float * b = &bmv.mvec[0];
  float * res = &r.mvec[0];
  res[0]=b[0]*a[0];
  res[1]=b[1]*a[0]+b[0]*a[1];
  res[2]=b[2]*a[0]+b[0]*a[2];
  res[3]=b[3]*a[0]+b[0]*a[3];
  res[4]=b[4]*a[0]+b[0]*a[4];
  res[5]=b[5]*a[0]+b[2]*a[1]-b[1]*a[2]+b[0]*a[5];
  res[6]=b[6]*a[0]+b[3]*a[1]-b[1]*a[3]+b[0]*a[6];
  res[7]=b[7]*a[0]+b[4]*a[1]-b[1]*a[4]+b[0]*a[7];
  res[8]=b[8]*a[0]+b[3]*a[2]-b[2]*a[3]+b[0]*a[8];
  res[9]=b[9]*a[0]-b[4]*a[2]+b[2]*a[4]+b[0]*a[9];
  res[10]=b[10]*a[0]+b[4]*a[3]-b[3]*a[4]+b[0]*a[10];
  res[11]=b[11]*a[0]-b[8]*a[1]+b[6]*a[2]-b[5]*a[3]-b[3]*a[5]+b[2]*a[6]-b[1]*a[8]+b[0]*a[11];
  res[12]=b[12]*a[0]-b[9]*a[1]-b[7]*a[2]+b[5]*a[4]+b[4]*a[5]-b[2]*a[7]-b[1]*a[9]+b[0]*a[12];
  res[13]=b[13]*a[0]-b[10]*a[1]+b[7]*a[3]-b[6]*a[4]-b[4]*a[6]+b[3]*a[7]-b[1]*a[10]+b[0]*a[13];
  res[14]=b[14]*a[0]+b[10]*a[2]+b[9]*a[3]+b[8]*a[4]+b[4]*a[8]+b[3]*a[9]+b[2]*a[10]+b[0]*a[14];
  res[15]=b[15]*a[0]+b[14]*a[1]+b[13]*a[2]+b[12]*a[3]+b[11]*a[4]+b[10]*a[5]+b[9]*a[6]+b[8]*a[7]+b[7]*a[8]+b[6]*a[9]+b[5]*a[10]-b[4]*a[11]-b[3]*a[12]-b[2]*a[13]-b[1]*a[14]+b[0]*a[15];
  return r;
}

//***********************
// PGA3D.Vee : res = a & b 
// The regressive product. (JOIN)
//***********************
static inline multivector_t PGA3D_vee (const multivector_t amv,const multivector_t bmv) {
  
  multivector_t r;
  const float * a = &amv.mvec[0];
  const float * b = &bmv.mvec[0];
  float * res = &r.mvec[0];

  res[15]=1*(a[15]*b[15]);
  res[14]=-1*(a[14]*-1*b[15]+a[15]*b[14]*-1);
  res[13]=-1*(a[13]*-1*b[15]+a[15]*b[13]*-1);
  res[12]=-1*(a[12]*-1*b[15]+a[15]*b[12]*-1);
  res[11]=-1*(a[11]*-1*b[15]+a[15]*b[11]*-1);
  res[10]=1*(a[10]*b[15]+a[13]*-1*b[14]*-1-a[14]*-1*b[13]*-1+a[15]*b[10]);
  res[9]=1*(a[9]*b[15]+a[12]*-1*b[14]*-1-a[14]*-1*b[12]*-1+a[15]*b[9]);
  res[8]=1*(a[8]*b[15]+a[11]*-1*b[14]*-1-a[14]*-1*b[11]*-1+a[15]*b[8]);
  res[7]=1*(a[7]*b[15]+a[12]*-1*b[13]*-1-a[13]*-1*b[12]*-1+a[15]*b[7]);
  res[6]=1*(a[6]*b[15]-a[11]*-1*b[13]*-1+a[13]*-1*b[11]*-1+a[15]*b[6]);
  res[5]=1*(a[5]*b[15]+a[11]*-1*b[12]*-1-a[12]*-1*b[11]*-1+a[15]*b[5]);
  res[4]=1*(a[4]*b[15]-a[7]*b[14]*-1+a[9]*b[13]*-1-a[10]*b[12]*-1-a[12]*-1*b[10]+a[13]*-1*b[9]-a[14]*-1*b[7]+a[15]*b[4]);
  res[3]=1*(a[3]*b[15]-a[6]*b[14]*-1-a[8]*b[13]*-1+a[10]*b[11]*-1+a[11]*-1*b[10]-a[13]*-1*b[8]-a[14]*-1*b[6]+a[15]*b[3]);
  res[2]=1*(a[2]*b[15]-a[5]*b[14]*-1+a[8]*b[12]*-1-a[9]*b[11]*-1-a[11]*-1*b[9]+a[12]*-1*b[8]-a[14]*-1*b[5]+a[15]*b[2]);
  res[1]=1*(a[1]*b[15]+a[5]*b[13]*-1+a[6]*b[12]*-1+a[7]*b[11]*-1+a[11]*-1*b[7]+a[12]*-1*b[6]+a[13]*-1*b[5]+a[15]*b[1]);
  res[0]=1*(a[0]*b[15]+a[1]*b[14]*-1+a[2]*b[13]*-1+a[3]*b[12]*-1+a[4]*b[11]*-1+a[5]*b[10]+a[6]*b[9]+a[7]*b[8]+a[8]*b[7]+a[9]*b[6]+a[10]*b[5]-a[11]*-1*b[4]-a[12]*-1*b[3]-a[13]*-1*b[2]-a[14]*-1*b[1]+a[15]*b[0]);
  return r;
}

//***********************
// PGA3D.Dot : res = a | b 
// The inner product.
//***********************
static inline multivector_t PGA3D_dot (const multivector_t amv,const multivector_t bmv) {
  multivector_t r;
  const float * a = &amv.mvec[0];
  const float * b = &bmv.mvec[0];
  float * res = &r.mvec[0];
  res[0]=b[0]*a[0]+b[2]*a[2]+b[3]*a[3]+b[4]*a[4]-b[8]*a[8]-b[9]*a[9]-b[10]*a[10]-b[14]*a[14];
  res[1]=b[1]*a[0]+b[0]*a[1]-b[5]*a[2]-b[6]*a[3]-b[7]*a[4]+b[2]*a[5]+b[3]*a[6]+b[4]*a[7]+b[11]*a[8]+b[12]*a[9]+b[13]*a[10]+b[8]*a[11]+b[9]*a[12]+b[10]*a[13]+b[15]*a[14]-b[14]*a[15];
  res[2]=b[2]*a[0]+b[0]*a[2]-b[8]*a[3]+b[9]*a[4]+b[3]*a[8]-b[4]*a[9]-b[14]*a[10]-b[10]*a[14];
  res[3]=b[3]*a[0]+b[8]*a[2]+b[0]*a[3]-b[10]*a[4]-b[2]*a[8]-b[14]*a[9]+b[4]*a[10]-b[9]*a[14];
  res[4]=b[4]*a[0]-b[9]*a[2]+b[10]*a[3]+b[0]*a[4]-b[14]*a[8]+b[2]*a[9]-b[3]*a[10]-b[8]*a[14];
  res[5]=b[5]*a[0]-b[11]*a[3]+b[12]*a[4]+b[0]*a[5]-b[15]*a[10]-b[3]*a[11]+b[4]*a[12]-b[10]*a[15];
  res[6]=b[6]*a[0]+b[11]*a[2]-b[13]*a[4]+b[0]*a[6]-b[15]*a[9]+b[2]*a[11]-b[4]*a[13]-b[9]*a[15];
  res[7]=b[7]*a[0]-b[12]*a[2]+b[13]*a[3]+b[0]*a[7]-b[15]*a[8]-b[2]*a[12]+b[3]*a[13]-b[8]*a[15];
  res[8]=b[8]*a[0]+b[14]*a[4]+b[0]*a[8]+b[4]*a[14];
  res[9]=b[9]*a[0]+b[14]*a[3]+b[0]*a[9]+b[3]*a[14];
  res[10]=b[10]*a[0]+b[14]*a[2]+b[0]*a[10]+b[2]*a[14];
  res[11]=b[11]*a[0]+b[15]*a[4]+b[0]*a[11]-b[4]*a[15];
  res[12]=b[12]*a[0]+b[15]*a[3]+b[0]*a[12]-b[3]*a[15];
  res[13]=b[13]*a[0]+b[15]*a[2]+b[0]*a[13]-b[2]*a[15];
  res[14]=b[14]*a[0]+b[0]*a[14];
  res[15]=b[15]*a[0]+b[0]*a[15];
  return r;
}

//***********************
// PGA3D.Add : res = a + b 
// Multivector addition
//***********************
static inline multivector_t PGA3D_add (const multivector_t amv,const multivector_t bmv) {
  multivector_t r;
  const float * a = &amv.mvec[0];
  const float * b = &bmv.mvec[0];
  float * res = &r.mvec[0];
  res[0] = a[0]+b[0];
    res[1] = a[1]+b[1];
    res[2] = a[2]+b[2];
    res[3] = a[3]+b[3];
    res[4] = a[4]+b[4];
    res[5] = a[5]+b[5];
    res[6] = a[6]+b[6];
    res[7] = a[7]+b[7];
    res[8] = a[8]+b[8];
    res[9] = a[9]+b[9];
    res[10] = a[10]+b[10];
    res[11] = a[11]+b[11];
    res[12] = a[12]+b[12];
    res[13] = a[13]+b[13];
    res[14] = a[14]+b[14];
    res[15] = a[15]+b[15];
  return r;
}

//***********************
// PGA3D.Sub : res = a - b 
// Multivector subtraction
//***********************
static inline multivector_t PGA3D_sub (const multivector_t amv,const multivector_t bmv) {
  multivector_t r;
  const float * a = &amv.mvec[0];
  const float * b = &bmv.mvec[0];
  float * res = &r.mvec[0];
      res[0] = a[0]-b[0];
    res[1] = a[1]-b[1];
    res[2] = a[2]-b[2];
    res[3] = a[3]-b[3];
    res[4] = a[4]-b[4];
    res[5] = a[5]-b[5];
    res[6] = a[6]-b[6];
    res[7] = a[7]-b[7];
    res[8] = a[8]-b[8];
    res[9] = a[9]-b[9];
    res[10] = a[10]-b[10];
    res[11] = a[11]-b[11];
    res[12] = a[12]-b[12];
    res[13] = a[13]-b[13];
    res[14] = a[14]-b[14];
    res[15] = a[15]-b[15];
  return r;
}

//***********************
// PGA3D.smul : res = a * b 
// scalar/multivector multiplication
//***********************
static inline multivector_t PGA3D_smul (const float a,const multivector_t bmv) {
  multivector_t r;
  const float * b = &bmv.mvec[0];
  float * res = &r.mvec[0];
      res[0] = a*b[0];
    res[1] = a*b[1];
    res[2] = a*b[2];
    res[3] = a*b[3];
    res[4] = a*b[4];
    res[5] = a*b[5];
    res[6] = a*b[6];
    res[7] = a*b[7];
    res[8] = a*b[8];
    res[9] = a*b[9];
    res[10] = a*b[10];
    res[11] = a*b[11];
    res[12] = a*b[12];
    res[13] = a*b[13];
    res[14] = a*b[14];
    res[15] = a*b[15];
  return r;
}

//***********************
// PGA3D.muls : res = a * b 
// multivector/scalar multiplication
//***********************
static inline multivector_t PGA3D_muls (const multivector_t amv, const float b) {
  multivector_t r;
  const float * a = &amv.mvec[0];
  float * res = &r.mvec[0];
      res[0] = a[0]*b;
    res[1] = a[1]*b;
    res[2] = a[2]*b;
    res[3] = a[3]*b;
    res[4] = a[4]*b;
    res[5] = a[5]*b;
    res[6] = a[6]*b;
    res[7] = a[7]*b;
    res[8] = a[8]*b;
    res[9] = a[9]*b;
    res[10] = a[10]*b;
    res[11] = a[11]*b;
    res[12] = a[12]*b;
    res[13] = a[13]*b;
    res[14] = a[14]*b;
    res[15] = a[15]*b;
  return r;
}

//***********************
// PGA3D.sadd : res = a + b 
// scalar/multivector addition
//***********************
static inline multivector_t PGA3D_sadd (const float a,const multivector_t bmv) {
  multivector_t r;
  const float * b = &bmv.mvec[0];
  float * res = &r.mvec[0];
    res[0] = a+b[0];
      res[1] = b[1];
    res[2] = b[2];
    res[3] = b[3];
    res[4] = b[4];
    res[5] = b[5];
    res[6] = b[6];
    res[7] = b[7];
    res[8] = b[8];
    res[9] = b[9];
    res[10] = b[10];
    res[11] = b[11];
    res[12] = b[12];
    res[13] = b[13];
    res[14] = b[14];
    res[15] = b[15];
  return r;
}

//***********************
// PGA3D.adds : res = a + b 
// multivector/scalar addition
//***********************
static inline multivector_t PGA3D_adds (const multivector_t amv, const float b) {
  multivector_t r;
  const float * a = &amv.mvec[0];
  float * res = &r.mvec[0];
    res[0] = a[0]+b;
      res[1] = a[1];
    res[2] = a[2];
    res[3] = a[3];
    res[4] = a[4];
    res[5] = a[5];
    res[6] = a[6];
    res[7] = a[7];
    res[8] = a[8];
    res[9] = a[9];
    res[10] = a[10];
    res[11] = a[11];
    res[12] = a[12];
    res[13] = a[13];
    res[14] = a[14];
    res[15] = a[15];
  return r;
}

//***********************
// PGA3D.ssub : res = a - b 
// scalar/multivector subtraction
//***********************
static inline multivector_t PGA3D_ssub (const float a,const multivector_t bmv) {
  multivector_t r;
  const float * b = &bmv.mvec[0];
  float * res = &r.mvec[0];
    res[0] = a-b[0];
      res[1] = -b[1];
    res[2] = -b[2];
    res[3] = -b[3];
    res[4] = -b[4];
    res[5] = -b[5];
    res[6] = -b[6];
    res[7] = -b[7];
    res[8] = -b[8];
    res[9] = -b[9];
    res[10] = -b[10];
    res[11] = -b[11];
    res[12] = -b[12];
    res[13] = -b[13];
    res[14] = -b[14];
    res[15] = -b[15];
  return r;
}

//***********************
// PGA3D.subs : res = a - b 
// multivector/scalar subtraction
//***********************
static inline multivector_t PGA3D_subs (const multivector_t amv, const float b) {
  multivector_t r;
  const float * a = &amv.mvec[0];
  float * res = &r.mvec[0];
    res[0] = a[0]-b;
      res[1] = a[1];
    res[2] = a[2];
    res[3] = a[3];
    res[4] = a[4];
    res[5] = a[5];
    res[6] = a[6];
    res[7] = a[7];
    res[8] = a[8];
    res[9] = a[9];
    res[10] = a[10];
    res[11] = a[11];
    res[12] = a[12];
    res[13] = a[13];
    res[14] = a[14];
    res[15] = a[15];
  return r;
}


static inline float PGA3D_norm(multivector_t mv) { 
    return sqrt(fabsf(PGA3D_mul(mv, PGA3D_conjugate(mv)).mvec[0])); 
}
static inline float PGA3D_inorm(multivector_t mv) { return PGA3D_norm(PGA3D_dual(mv)); }
static inline multivector_t PGA3D_normalized(multivector_t mv) { return PGA3D_muls(mv, 1.0f/PGA3D_norm(mv)); }


// A rotor (Euclidean line) and translator (Ideal line)
static inline multivector_t PGA3D_rotor(float angle, multivector_t line) { return PGA3D_sadd(cosf(angle/2.0f), PGA3D_smul(sinf(angle/2.0f),PGA3D_normalized(line))); }

// PGA is plane based. Vectors are planes. (think linear functionals)
static const multivector_t e0(void) {return (multivector_t){{0,1.0f,0,0,0,0,0,0,0,0,0,0,0,0,0,0}};}
static const multivector_t e1(void) {return (multivector_t){{0,0,1.0f,0,0,0,0,0,0,0,0,0,0,0,0,0}};}
static const multivector_t e2(void) {return (multivector_t){{0,0,0,1.0f,0,0,0,0,0,0,0,0,0,0,0,0}};}
static const multivector_t e3(void) {return (multivector_t){{0,0,0,0,1.0f,0,0,0,0,0,0,0,0,0,0,0}};}


// PGA points are trivectors.
static const multivector_t e123(void) { return PGA3D_wedge(e1(), PGA3D_wedge(e2(), e3())); }
static const multivector_t e032(void) { return PGA3D_wedge(e0(), PGA3D_wedge(e3(), e2())); }
static const multivector_t e013(void) { return PGA3D_wedge(e0(), PGA3D_wedge(e1(), e3())); }
static const multivector_t e021(void) { return PGA3D_wedge(e0(), PGA3D_wedge(e2(), e1())); }

// A point is just a homogeneous point, euclidean coordinates plus the origin
static multivector_t point(float x, float y, float z) { return PGA3D_add(e123(), PGA3D_add(PGA3D_smul(x,e032()), PGA3D_add(PGA3D_smul(y,e013()), PGA3D_smul(z,e021())))); }

/*
// A plane is defined using its homogenous equation ax + by + cz + d = 0
static PGA3D plane(float a,float b,float c,float d) { return a*e1 + b*e2 + c*e3 + d*e0; }
static PGA3D translator(float dist, PGA3D line) { return 1.0f + dist/2.0f*line; }
// for our toy problem (generate points on the surface of a torus)
// we start with a function that generates motors.
// circle(t) with t going from 0 to 1.
static PGA3D circle(float t, float radius, PGA3D line) {
  return rotor(t*2.0f*PI,line) * translator(radius,e1*e0);
}

// a torus is now the product of two circles.
static PGA3D torus(float s, float t, float r1, PGA3D l1, float r2, PGA3D l2) {
  return circle(s,r2,l2)*circle(t,r1,l1);
}

// and to sample its points we simply sandwich the origin ..
static PGA3D point_on_torus(float s, float t) {
  PGA3D to = torus(s,t,0.25f,e1*e2,0.6f,e1*e3);
  return to * e123 * ~to;
}


int main (int argc, char **argv) {
  
  // Elements of the even subalgebra (scalar + bivector + pss) of unit length are motors
  PGA3D rot = rotor( PI/2.0f, e1 * e2 );
  
  // The outer product ^ is the MEET. Here we intersect the yz (x=0) and xz (y=0) planes.
  PGA3D ax_z = e1 ^ e2;
  
  // line and plane meet in point. We intersect the line along the z-axis (x=0,y=0) with the xy (z=0) plane.
  PGA3D orig = ax_z ^ e3;
  
  // We can also easily create points and join them into a line using the regressive (vee, &) product.
  PGA3D px = point(1.0,0.0,0.0);
  PGA3D line = orig & px;
  
  // Lets also create the plane with equation 2x + z - 3 = 0
  PGA3D p = plane(2,0,1,-3);
  
  // rotations work on all elements
  PGA3D rotated_plane = rot * p * ~rot;
  PGA3D rotated_line  = rot * line * ~rot;
  PGA3D rotated_point = rot * px * ~rot;
  
  // See the 3D PGA Cheat sheet for a huge collection of useful formulas
  PGA3D point_on_plane = (p | px) * p;
  
  // Some output.
  printf("a point       : "); px.log();
  printf("a line        : "); line.log();
  printf("a plane       : "); p.log();
  printf("a rotor       : "); rot.log();
  printf("rotated line  : "); rotated_line.log();
  printf("rotated point : "); rotated_point.log();
  printf("rotated plane : "); rotated_plane.log();
  printf("point on plane: "); point_on_plane.normalized().log();
  printf("point on torus: "); point_on_torus(0.0f,0.0f).log();
  (e0-1.0f).log();
  (1.0f-e0).log();

  return 0;
}
*/

#endif /* PGA3D_H */
