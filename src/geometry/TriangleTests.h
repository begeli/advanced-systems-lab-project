#pragma once

/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * Source: https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/opttritri.txt
 * Quite heavily modified to make it compatible with our data formats.
*/

#include <math.h>
#include "LinAlg.h"
#include "../include/bvhstats.h"

/* if USE_EPSILON_TEST is true then we do a check:
         if |dv|<EPSILON then dv=0.0;
   else no check is done (which is less robust)
*/
#define USE_EPSILON_TEST TRUE
#define EPSILON 0.000001


/* some macros */
//#define CROSS(dest, v1, v2){                     \
//              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
//              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
//              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];}
//
//#define DOT(v1, v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
//
//#define SUB(dest, v1, v2){         \
//            dest[0]=v1[0]-v2[0]; \
//            dest[1]=v1[1]-v2[1]; \
//            dest[2]=v1[2]-v2[2];}

/* sort so that a<=b */
//#define SORT(a, b)       \
//             if(a>b)    \
//             {          \
//               float c; \
//               c=a;     \
//               a=b;     \
//               b=c;     \
//             }
//

/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */
//#define EDGE_EDGE_TEST(V0, U0, U1)                      \
//  Bx = U0[i0]-U1[i0];                                   \
//  By=U0[i1]-U1[i1];                                   \
//  Cx=V0[i0]-U0[i0];                                   \
//  Cy=V0[i1]-U0[i1];                                   \
//  f=Ay*Bx-Ax*By;                                      \
//  d=By*Cx-Bx*Cy;                                      \
//  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
//  {                                                   \
//    e=Ax*Cy-Ay*Cx;                                    \
//    if(f>0)                                           \
//    {                                                 \
//      if(e>=0 && e<=f) return 1;                      \
//    }                                                 \
//    else                                              \
//    {                                                 \
//      if(e<=0 && e>=f) return 1;                      \
//    }                                                 \
//  }

//#define EDGE_AGAINST_TRI_EDGES(V0, V1, U0, U1, U2) \
//{                                              \
//  float Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
//  Ax=V1[i0]-V0[i0];                            \
//  Ay=V1[i1]-V0[i1];                            \
//  /* test edge U0,U1 against V0,V1 */          \
//  EDGE_EDGE_TEST(V0,U0,U1);                    \
//  /* test edge U1,U2 against V0,V1 */          \
//  EDGE_EDGE_TEST(V0,U1,U2);                    \
//  /* test edge U2,U1 against V0,V1 */          \
//  EDGE_EDGE_TEST(V0,U2,U0);                    \
//}
//
//#define POINT_IN_TRI(V0, U0, U1, U2)           \
//{                                           \
//  float a,b,c,d0,d1,d2;                     \
//  /* is T1 completly inside T2? */          \
//  /* check if V0 is inside tri(U0,U1,U2) */ \
//  a=U1[i1]-U0[i1];                          \
//  b=-(U1[i0]-U0[i0]);                       \
//  c=-a*U0[i0]-b*U0[i1];                     \
//  d0=a*V0[i0]+b*V0[i1]+c;                   \
//                                            \
//  a=U2[i1]-U1[i1];                          \
//  b=-(U2[i0]-U1[i0]);                       \
//  c=-a*U1[i0]-b*U1[i1];                     \
//  d1=a*V0[i0]+b*V0[i1]+c;                   \
//                                            \
//  a=U0[i1]-U2[i1];                          \
//  b=-(U0[i0]-U2[i0]);                       \
//  c=-a*U2[i0]-b*U2[i1];                     \
//  d2=a*V0[i0]+b*V0[i1]+c;                   \
//  if(d0*d1>0.0)                             \
//  {                                         \
//    if(d0*d2>0.0) return 1;                 \
//  }                                         \
//}

__always_inline float gc(struct Vec3f p, short i) {
  assume(0 <= i && i <= 2);
  switch (i) {
    case 0: return p.x;
    case 1: return p.y;
    case 2: return p.z;
    default:;
  }
  return FLT_MIN;
}

__always_inline int coplanar_tri_tri(struct Vec3f N, struct Vec3f V0, struct Vec3f V1, struct Vec3f V2,
                     struct Vec3f U0, struct Vec3f U1, struct Vec3f U2) {
  short i0, i1;
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 4;
#endif
  /* first project onto an axis-aligned plane, that maximizes the area */
  /* of the triangles, compute indices: i0,i1. */
  struct Vec3f A = vabs3D(N);
  if (A.x > A.y) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 1;
#endif
    if (A.x > A.z) {
      i0 = 1;      /* A[0] is greatest */
      i1 = 2;
    } else {
      i0 = 0;      /* A[2] is greatest */
      i1 = 1;
    }
  } else   /* A[0]<=A[1] */
  {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 1;
#endif
    if (A.z > A.y) {
      i0 = 0;      /* A[2] is greatest */
      i1 = 1;
    } else {
      i0 = 0;      /* A[1] is greatest */
      i1 = 2;
    }
  }

  /* test all edges of triangle 1 against the edges of triangle 2 */
  //EDGE_AGAINST_TRI_EDGES(V0, V1, U0, U1, U2);
  {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 18;
#endif
    float Ax, Ay, Bx, By, Cx, Cy, e, d, f;
    //Ax = V1[i0] - V0[i0];
    Ax = gc(V1, i0) - gc(V0, i0);
    //Ay = V1[i1] - V0[i1];
    Ay = gc(V1, i1) - gc(V0, i1);
    //Bx = U0[i0] - U1[i0];
    Bx = gc(U0, i0) - gc(U1, i0);
    //By = U0[i1] - U1[i1];
    By = gc(U0, i1) - gc(U1, i1);
    //Cx = V0[i0] - U0[i0];
    Cx = gc(V0, i0) - gc(U0, i0);
    //Cy = V0[i1] - U0[i1];
    Cy = gc(V0, i1) - gc(U0, i1);

    f = Ay * Bx - Ax * By;
    d = By * Cx - Bx * Cy;
    if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 4;
#endif
      e = Ax * Cy - Ay * Cx;
      if (f > 0) {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e >= 0 && e <= f)
          return 1;
      }
      else {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e <= 0 && e >= f)
          return 1;
      }
    }
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 16;
#endif
    //Bx = U1[i0] - U2[i0];
    Bx = gc(U1, i0) - gc(U2, i0);
    //By = U1[i1] - U2[i1];
    By = gc(U1, i1) - gc(U2, i1);
    //Cx = V0[i0] - U1[i0];
    Cx = gc(V0, i0) - gc(U1, i0);
    //Cy = V0[i1] - U1[i1];
    Cy = gc(V0, i1) - gc(U1, i1);

    f = Ay * Bx - Ax * By;
    d = By * Cx - Bx * Cy;
    if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 4;
#endif
      e = Ax * Cy - Ay * Cx;
      if (f > 0) {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e >= 0 && e <= f)
          return 1;
      } else {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e <= 0 && e >= f)
          return 1;
      }
    }
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 16;
#endif
    //Bx = U2[i0] - U0[i0];
    Bx = gc(U2, i0) - gc(U0, i0);
    //By = U2[i1] - U0[i1];
    By = gc(U2, i1) - gc(U0, i1);
    //Cx = V0[i0] - U2[i0];
    Cx = gc(V0, i0) - gc(U2, i0);
    //Cy = V0[i1] - U2[i1];
    Cy = gc(V0, i1) - gc(U2, i1);

    f = Ay * Bx - Ax * By;
    d = By * Cx - Bx * Cy;
    if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 4;
#endif
      e = Ax * Cy - Ay * Cx;
      if (f > 0) {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e >= 0 && e <= f)
          return 1;
      } else {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e <= 0 && e >= f)
          return 1;
      }
    }
  }

  //EDGE_AGAINST_TRI_EDGES(V1, V2, U0, U1, U2);
  {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 18;
#endif
    float Ax, Ay, Bx, By, Cx, Cy, e, d, f;
    Ax = gc(V2, i0) - gc(V1, i0);
    Ay = gc(V2, i1) - gc(V1, i1);
    Bx = gc(U0, i0) - gc(U1, i0);
    By = gc(U0, i1) - gc(U1, i1);
    Cx = gc(V1, i0) - gc(U0, i0);
    Cy = gc(V1, i1) - gc(U0, i1);
    f = Ay * Bx - Ax * By;
    d = By * Cx - Bx * Cy;
    if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 4;
#endif
      e = Ax * Cy - Ay * Cx;
      if (f > 0) {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e >= 0 && e <= f)
          return 1;
      } else {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e <= 0 && e >= f)
          return 1;
      }
    };
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 16;
#endif
    Bx = gc(U1, i0) - gc(U2, i0);
    By = gc(U1, i1) - gc(U2, i1);
    Cx = gc(V1, i0) - gc(U1, i0);
    Cy = gc(V1, i1) - gc(U1, i1);
    f = Ay * Bx - Ax * By;
    d = By * Cx - Bx * Cy;
    if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 4;
#endif
      e = Ax * Cy - Ay * Cx;
      if (f > 0) {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e >= 0 && e <= f)
          return 1;
      } else {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e <= 0 && e >= f)
          return 1;
      }
    };
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 16;
#endif
    Bx = gc(U2, i0) - gc(U0, i0);
    By = gc(U2, i1) - gc(U0, i1);
    Cx = gc(V1, i0) - gc(U2, i0);
    Cy = gc(V1, i1) - gc(U2, i1);
    f = Ay * Bx - Ax * By;
    d = By * Cx - Bx * Cy;
    if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 4;
#endif
      e = Ax * Cy - Ay * Cx;
      if (f > 0) {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e >= 0 && e <= f)
          return 1;
      } else {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e <= 0 && e >= f)
          return 1;
      }
    };
  }

  //EDGE_AGAINST_TRI_EDGES(V2, V0, U0, U1, U2);
  {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 18;
#endif
    float Ax, Ay, Bx, By, Cx, Cy, e, d, f;
    Ax = gc(V0, i0) - gc(V2, i0);
    Ay = gc(V0, i1) - gc(V2, i1);
    Bx = gc(U0, i0) - gc(U1, i0);
    By = gc(U0, i1) - gc(U1, i1);
    Cx = gc(V2, i0) - gc(U0, i0);
    Cy = gc(V2, i1) - gc(U0, i1);
    f = Ay * Bx - Ax * By;
    d = By * Cx - Bx * Cy;
    if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 4;
#endif
      e = Ax * Cy - Ay * Cx;
      if (f > 0) {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e >= 0 && e <= f)
          return 1;
      } else {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e <= 0 && e >= f)
          return 1;
      }
    };
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 16;
#endif
    Bx = gc(U1, i0) - gc(U2, i0);
    By = gc(U1, i1) - gc(U2, i1);
    Cx = gc(V2, i0) - gc(U1, i0);
    Cy = gc(V2, i1) - gc(U1, i1);
    f = Ay * Bx - Ax * By;
    d = By * Cx - Bx * Cy;
    if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 4;
#endif
      e = Ax * Cy - Ay * Cx;
      if (f > 0) {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e >= 0 && e <= f)
          return 1;
      } else {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e <= 0 && e >= f)
          return 1;
      }
    };
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 16;
#endif
    Bx = gc(U2, i0) - gc(U0, i0);
    By = gc(U2, i1) - gc(U0, i1);
    Cx = gc(V2, i0) - gc(U2, i0);
    Cy = gc(V2, i1) - gc(U2, i1);
    f = Ay * Bx - Ax * By;
    d = By * Cx - Bx * Cy;
    if ((f > 0 && d >= 0 && d <= f) || (f < 0 && d <= 0 && d >= f)) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 4;
#endif
      e = Ax * Cy - Ay * Cx;
      if (f > 0) {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e >= 0 && e <= f)
          return 1;
      } else {
#ifdef BVH_COUNT_FLOPS
        bvh_flopcount += 2;
#endif
        if (e <= 0 && e >= f)
          return 1;
      }
    };
  }

  /* finally, test if tri1 is totally contained in tri2 or vice versa */
  //POINT_IN_TRI(V0, U0, U1, U2);
  {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 35;
#endif
    float a, b, c, d0, d1, d2;
    a = gc(U1, i1) - gc(U0, i1);
    b = -(gc(U1, i0) - gc(U0, i0));
    c = -a * gc(U0, i0) - b * gc(U0, i1);
    d0 = a * gc(V0, i0) + b * gc(V0, i1) + c;
    a = gc(U2, i1) - gc(U1, i1);
    b = -(gc(U2, i0) - gc(U1, i0));
    c = -a * gc(U1, i0) - b * gc(U1, i1);
    d1 = a * gc(V0, i0) + b * gc(V0, i1) + c;
    a = gc(U0, i1) - gc(U2, i1);
    b = -(gc(U0, i0) - gc(U2, i0));
    c = -a * gc(U2, i0) - b * gc(U2, i1);
    d2 = a * gc(V0, i0) + b * gc(V0, i1) + c;
    if (d0 * d1 > 0.0) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 2;
#endif
      if (d0 * d2 > 0.0)return 1;
    }
  }

  //POINT_IN_TRI(U0, V0, V1, V2);
  {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 35;
#endif
    float a, b, c, d0, d1, d2;
    a = gc(V1, i1) - gc(V0, i1);
    b = -(gc(V1, i0) - gc(V0, i0));
    c = -a * gc(V0, i0) - b * gc(V0, i1);
    d0 = a * gc(U0, i0) + b * gc(U0, i1) + c;
    a = gc(V2, i1) - gc(V1, i1);
    b = -(gc(V2, i0) - gc(V1, i0));
    c = -a * gc(V1, i0) - b * gc(V1, i1);
    d1 = a * gc(U0, i0) + b * gc(U0, i1) + c;
    a = gc(V0, i1) - gc(V2, i1);
    b = -(gc(V0, i0) - gc(V2, i0));
    c = -a * gc(V2, i0) - b * gc(V2, i1);
    d2 = a * gc(U0, i0) + b * gc(U0, i1) + c;
    if (d0 * d1 > 0.0) {
#ifdef BVH_COUNT_FLOPS
      bvh_flopcount += 2;
#endif
      if (d0 * d2 > 0.0)return 1;
    }
  }

  return 0;
}

//#define NEWCOMPUTE_INTERVALS(VV0, VV1, VV2, D0, D1, D2, D0D1, D0D2, A, B, C, X0, X1) \
//{ \
//        if(D0D1>0.0f) \
//        { \
//                /* here we know that D0D2<=0.0 */ \
//            /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
//                A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
//        } \
//        else if(D0D2>0.0f)\
//        { \
//                /* here we know that d0d1<=0.0 */ \
//            A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
//        } \
//        else if(D1*D2>0.0f || D0!=0.0f) \
//        { \
//                /* here we know that d0d1<=0.0 or that D0!=0.0 */ \
//                A=VV0; B=(VV1-VV0)*D0; C=(VV2-VV0)*D0; X0=D0-D1; X1=D0-D2; \
//        } \
//        else if(D1!=0.0f) \
//        { \
//                A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
//        } \
//        else if(D2!=0.0f) \
//        { \
//                A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
//        } \
//        else \
//        { \
//                /* triangles are coplanar */ \
//                return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2); \
//        } \
//}

__always_inline int NoDivTriTriIsect(struct Vec3f V0, struct Vec3f V1, struct Vec3f V2,
                     struct Vec3f U0, struct Vec3f U1, struct Vec3f U2) {
  //struct Vec3f E1, E2;
  //struct Vec3f N1, N2;
  float d1, d2;
  float du0, du1, du2, dv0, dv1, dv2;
  float isect1[2], isect2[2];
  float du0du1, du0du2, dv0dv1, dv0dv2;
  float vp0, vp1, vp2;
  float up0, up1, up2;
  float bb, cc, max;
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 41;
#endif

  /* compute plane equation of triangle(V0,V1,V2) */
  struct Vec3f E1 = vsub(V1, V0);
  struct Vec3f E2 = vsub(V2, V0);
  struct Vec3f N1 = vcross3D(E1, E2);
  d1 = -vdot3D(N1, V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0 = vdot3D(N1, U0) + d1;
  du1 = vdot3D(N1, U1) + d1;
  du2 = vdot3D(N1, U2) + d1;

  /* coplanarity robustness check */
#if USE_EPSILON_TEST == TRUE
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 6;
#endif
  if (fabsf(du0) < EPSILON) du0 = 0.0;
  if (fabsf(du1) < EPSILON) du1 = 0.0;
  if (fabsf(du2) < EPSILON) du2 = 0.0;
#endif
  du0du1 = du0 * du1;
  du0du2 = du0 * du2;

  if (du0du1 > 0.0f && du0du2 > 0.0f) {/* same sign on all of them + not equal 0 ? */
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 2;
#endif
    return 0;                    /* no intersection occurs */
  }
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 41;
#endif

  /* compute plane of triangle (U0,U1,U2) */
  E1 = vsub(U1, U0);
  E2 = vsub(U2, U0);
  struct Vec3f N2 = vcross3D(E1, E2);
  d2 = -vdot3D(N2, U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0 = vdot3D(N2, V0) + d2;
  dv1 = vdot3D(N2, V1) + d2;
  dv2 = vdot3D(N2, V2) + d2;

#if USE_EPSILON_TEST == TRUE
#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 6;
#endif
  if (fabsf(dv0) < EPSILON) dv0 = 0.0f;
  if (fabsf(dv1) < EPSILON) dv1 = 0.0f;
  if (fabsf(dv2) < EPSILON) dv2 = 0.0f;
#endif

  dv0dv1 = dv0 * dv1;
  dv0dv2 = dv0 * dv2;

  if (dv0dv1 > 0.0f && dv0dv2 > 0.0f) { /* same sign on all of them + not equal 0 ? */
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 2;
#endif
    return 0;                    /* no intersection occurs */
  }

#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 15;
#endif

  /* compute direction of intersection line */
  struct Vec3f D = vcross3D(N1, N2);

  /* compute and index to the largest component of D */
  max = fabsf(D.x);
  bb = fabsf(D.y);
  cc = fabsf(D.z);

  /* this is the simplified projection onto L*/
  vp0 = V0.x;
  vp1 = V1.x;
  vp2 = V2.x;

  up0 = U0.x;
  up1 = U1.x;
  up2 = U2.x;

  if (bb > max) {
    max = bb;

    /* this is the simplified projection onto L*/
    vp0 = V0.y;
    vp1 = V1.y;
    vp2 = V2.y;

    up0 = U0.y;
    up1 = U1.y;
    up2 = U2.y;
  }
  if (cc > max) {
    max = cc;

    /* this is the simplified projection onto L*/
    vp0 = V0.z;
    vp1 = V1.z;
    vp2 = V2.z;

    up0 = U0.z;
    up1 = U1.z;
    up2 = U2.z;
  }

  /* compute interval for triangle 1 */
  float a, b, c, x0, x1;
  //NEWCOMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,a,b,c,x0,x1);
  if (dv0dv1 > 0.0f) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 6;
#endif
    a = vp2;
    b = (vp0 - vp2) * dv2;
    c = (vp1 - vp2) * dv2;
    x0 = dv2 - dv0;
    x1 = dv2 - dv1;
  } else if (dv0dv2 > 0.0f) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 7;
#endif
    a = vp1;
    b = (vp0 - vp1) * dv1;
    c = (vp2 - vp1) * dv1;
    x0 = dv1 - dv0;
    x1 = dv1 - dv2;
  } else if (dv1 * dv2 > 0.0f || dv0 != 0.0f) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 10;
#endif
    a = vp0;
    b = (vp1 - vp0) * dv0;
    c = (vp2 - vp0) * dv0;
    x0 = dv0 - dv1;
    x1 = dv0 - dv2;
  } else if (dv1 != 0.0f) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 11;
#endif
    a = vp1;
    b = (vp0 - vp1) * dv1;
    c = (vp2 - vp1) * dv1;
    x0 = dv1 - dv0;
    x1 = dv1 - dv2;
  } else if (dv2 != 0.0f) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 12;
#endif
    a = vp2;
    b = (vp0 - vp2) * dv2;
    c = (vp1 - vp2) * dv2;
    x0 = dv2 - dv0;
    x1 = dv2 - dv1;
  } else {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 6;
#endif
    return coplanar_tri_tri(N1, V0, V1, V2, U0, U1, U2);
  }

#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 1;
#endif
  /* compute interval for triangle 2 */
  float d, e, f, y0, y1;
  // NEWCOMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,d,e,f,y0,y1);
  if (du0du1 > 0.0f) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 6;
#endif
    d = up2;
    e = (up0 - up2) * du2;
    f = (up1 - up2) * du2;
    y0 = du2 - du0;
    y1 = du2 - du1;
  } else if (du0du2 > 0.0f) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 7;
#endif
    d = up1;
    e = (up0 - up1) * du1;
    f = (up2 - up1) * du1;
    y0 = du1 - du0;
    y1 = du1 - du2;
  } else if (du1 * du2 > 0.0f || du0 != 0.0f) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 10;
#endif
    d = up0;
    e = (up1 - up0) * du0;
    f = (up2 - up0) * du0;
    y0 = du0 - du1;
    y1 = du0 - du2;
  } else if (du1 != 0.0f) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 11;
#endif
    d = up1;
    e = (up0 - up1) * du1;
    f = (up2 - up1) * du1;
    y0 = du1 - du0;
    y1 = du1 - du2;
  } else if (du2 != 0.0f) {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 12;
#endif
    d = up2;
    e = (up0 - up2) * du2;
    f = (up1 - up2) * du2;
    y0 = du2 - du0;
    y1 = du2 - du1;
  } else {
#ifdef BVH_COUNT_FLOPS
    bvh_flopcount += 6;
#endif
    return coplanar_tri_tri(N1, V0, V1, V2, U0, U1, U2);
  }

#ifdef BVH_COUNT_FLOPS
  bvh_flopcount += 21;
#endif
  float xx, yy, xxyy, tmp;
  xx = x0 * x1;
  yy = y0 * y1;
  xxyy = xx * yy;

  tmp = a * xxyy;
  isect1[0] = tmp + b * x1 * yy;
  isect1[1] = tmp + c * x0 * yy;

  tmp = d * xxyy;
  isect2[0] = tmp + e * xx * y1;
  isect2[1] = tmp + f * xx * y0;

  //SORT(isect1[0],isect1[1]);
  if (isect1[0] > isect1[1]) {
    float temp;
    temp = isect1[0];
    isect1[0] = isect1[1];
    isect1[1] = temp;
  }
  //SORT(isect2[0],isect2[1]);
  if (isect2[0] > isect2[1]) {
    float temp;
    temp = isect2[0];
    isect2[0] = isect2[1];
    isect2[1] = temp;
  }

  if (isect1[1] < isect2[0] || isect2[1] < isect1[0]) return 0;
  return 1;
}
