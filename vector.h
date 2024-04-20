#pragma once
#include <math.h>
#include <stdarg.h>
#include <stdio.h>

typedef struct vector4_t {
  double x;
  double y;
  double z;
  double w;
} vector4_t;

typedef struct vector3_t {
  double x;
  double y;
  double z;
} vector3_t;

typedef struct vector2_t {
  double x;
  double y;
} vector2_t;

typedef struct matrix4x4_t {
  double entries[4][4];
} matrix4x4_t;

#define DOT(v, u) (v.x * u.x + v.y * u.y + v.z * u.z)
#define CROSS(v, u)                                                            \
  (vector3_t) {                                                                \
    (v.y * u.z - v.z * u.y), -(v.x * u.z - v.z * u.x), (v.x * u.y - v.y * u.x) \
  }
#define NORM(u)                                                                \
  (vector3_t) {                                                                \
    u.x / (sqrtf(u.x * u.x + u.y * u.y + u.z * u.z)),                          \
        u.y / (sqrtf(u.x * u.x + u.y * u.y + u.z * u.z)),                      \
        u.z / (sqrtf(u.x * u.x + u.y * u.y + u.z * u.z))                       \
  }

#define ADD(v, u)                                                              \
  (vector3_t) { v.x + u.x, v.y + u.y, v.z + u.z }

#define SUB(v, u)                                                              \
  (vector3_t) { v.x - u.x, v.y - u.y, v.z - u.z }

#define LEN(u) (sqrtf(u.x * u.x + u.y * u.y + u.z * u.z))
#define MUL(u, s)                                                              \
  (vector3_t) { u.x *s, u.y *s, u.z *s }

static inline void m4x4print(matrix4x4_t A);
static inline vector3_t v3m4x4mul(vector3_t v, matrix4x4_t A);
static inline matrix4x4_t m4x4mul(int argcount, ...);

static inline vector2_t v3tov2(vector3_t v) {
  return (vector2_t){v.x / v.z, v.y / v.z};
}

static inline vector4_t v3tov4(vector3_t v) {
  return (vector4_t){v.x, v.y, v.z, 1};
}

static inline matrix4x4_t m4x4create(double m[4][4]) {
  matrix4x4_t ret;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      ret.entries[i][j] = m[i][j];
    }
  }
  return ret;
}

static inline vector4_t v4m4x4mul(vector4_t v, matrix4x4_t A) {
  vector4_t ret;
  ret.x = A.entries[0][0] * v.x + A.entries[1][0] * v.y +
          A.entries[2][0] * v.z + A.entries[3][0] * v.w;
  ret.y = A.entries[0][1] * v.x + A.entries[1][1] * v.y +
          A.entries[2][1] * v.z + A.entries[3][1] * v.w;
  ret.z = A.entries[0][2] * v.x + A.entries[1][2] * v.y +
          A.entries[2][2] * v.z + A.entries[3][2] * v.w;
  ret.w = A.entries[0][3] * v.x + A.entries[1][3] * v.y +
          A.entries[2][3] * v.z + A.entries[3][3] * v.w;

  return ret;
}

static inline vector4_t m4x4v4mul(matrix4x4_t A, vector4_t v) {
  vector4_t ret;
  ret.x = v.x * A.entries[0][0] + v.y * A.entries[0][1] +
          v.z * A.entries[0][2] + v.w * A.entries[0][3];
  ret.y = v.x * A.entries[1][0] + v.y * A.entries[1][1] +
          v.z * A.entries[1][2] + v.w * A.entries[1][3];
  ret.z = v.x * A.entries[2][0] + v.y * A.entries[2][1] +
          v.z * A.entries[2][2] + v.w * A.entries[2][3];
  ret.w = v.x * A.entries[3][0] + v.y * A.entries[3][1] +
          v.z * A.entries[3][2] + v.w * A.entries[3][3];
  return ret;
}

static inline matrix4x4_t _m4x4mul_(matrix4x4_t A, matrix4x4_t B) {
  matrix4x4_t ret;
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      ret.entries[i][j] = 0;

  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++) {
      for (int k = 0; k < 4; k++)
        ret.entries[i][j] += (A.entries[i][k] * B.entries[k][j]);
    }
  return ret;
}

static inline matrix4x4_t m4x4mul(int argcount, ...) {

  va_list argptr;
  va_start(argptr, argcount);

  matrix4x4_t ret = va_arg(argptr, matrix4x4_t);

  for (int counter = 1; counter < argcount; counter++) {
    ret = _m4x4mul_(ret, va_arg(argptr, matrix4x4_t));
  }

  va_end(argptr);

  return ret;
}

static inline void v2print(vector2_t v) { printf("<%lf, %lf>\n", v.x, v.y); }

static inline void v4print(vector4_t v) {
  printf("<%lf, %lf, %lf, %lf>\n", v.x, v.y, v.z, v.w);
}

static inline void v3print(vector3_t v) {
  printf("<%lf, %lf, %lf>\n", v.x, v.y, v.z);
}

static inline void m4x4print(matrix4x4_t A) {
  printf("[");
  for (int i = 0; i < 4; i++) {
    if (i != 0)
      printf(" ");
    for (int j = 0; j < 4; j++) {
      printf(" %lf ", A.entries[i][j]);
    }
    if (i == 3)
      printf("]");
    printf("\n");
  }
}

static inline vector2_t v4tov2(vector4_t proj) {
  vector3_t p;
  vector2_t ret;
  p.x = proj.x;
  p.y = proj.y;
  p.z = proj.z;
  if (proj.w != 0) {
    p.x /= proj.w;
    p.y /= proj.w;
    p.z /= proj.w;
  }
  ret.x = p.x;
  ret.y = p.y;
  /*
  if (p.z != 0) {
    ret.x /= p.z;
    ret.y /= p.z;
  }
  */
  return ret;
}

static inline vector3_t v4tov3(vector4_t v) {
  return (vector3_t){v.x, v.y, v.z};
}

static inline vector3_t v3m4x4mul(vector3_t v, matrix4x4_t A) {
  return v4tov3(v4m4x4mul(v3tov4(v), A));
}

static inline matrix4x4_t createRotationMatrix(double x, double y, double z) {
  double yrotation[4][4] = {
      {cos(y), 0, sin(y), 0},
      {0, 1, 0, 0},
      {-sin(y), 0, cos(y), 0},
      {0, 0, 0, 0},
  };
  matrix4x4_t yrot = m4x4create(yrotation);

  double xrotation[4][4] = {
      {1, 0, 0, 0},
      {0, cos(x), -sin(x), 0},
      {0, sin(x), cos(x), 0},
      {0, 0, 0, 1},
  };
  matrix4x4_t xrot = m4x4create(xrotation);

  double zrotation[4][4] = {
      {cos(z), -sin(z), 0, 0},
      {sin(z), cos(z), 0, 0},
      {0, 0, 1, 0},
      {0, 0, 0, 1},
  };
  matrix4x4_t zrot = m4x4create(zrotation);

  return m4x4mul(3, xrot, yrot, zrot);
}

static inline matrix4x4_t createTranslationMatrix(double x, double y,
                                                  double z) {

  double entries2[4][4] = {
      {1, 0, 0, 0},
      {0, 1, 0, 0},
      {0, 0, 1, 0},
      {x, y, z, 1},
  };
  return m4x4create(entries2);
}

static inline matrix4x4_t createProjectionMatrix(double near, double far,
                                                 double a, double fov) {
  double entries5[4][4] = {
      {a / tan(fov / 2), 0, 0, 0},
      {0, 1 / tan(fov / 2), 0, 0},
      {0, 0, (far) / (far - near), 1},
      {0, 0, (far * near) / (far - near), 0},
  };
  return m4x4create(entries5);
}

static inline double _det2(double A[2][2]) {
  return A[0][0] * A[1][1] - A[0][1] * A[1][0];
}

static inline double _det3(double A[3][3]) {
  double m[2][2];
  double sum = 0;

  // a00 column
  m[0][0] = A[1][1];
  m[1][0] = A[2][1];
  m[0][1] = A[1][2];
  m[1][1] = A[2][2];
  sum += A[0][0] * _det2(m);

  // a01 column
  m[0][0] = A[1][0];
  m[1][0] = A[2][0];
  m[0][1] = A[1][2];
  m[1][1] = A[2][2];
  sum -= A[0][1] * _det2(m);

  // a02 column
  m[0][0] = A[1][0];
  m[1][0] = A[2][0];
  m[0][1] = A[1][1];
  m[1][1] = A[2][1];
  sum += A[0][2] * _det2(m);

  return sum;
}

static inline double _det4(double A[4][4]) {
  double m[3][3];
  double sum = 0;

  // a00 column
  m[0][0] = A[1][1];
  m[1][0] = A[2][1];
  m[2][0] = A[3][1];
  m[0][1] = A[1][2];
  m[1][1] = A[2][2];
  m[2][1] = A[3][2];
  m[0][2] = A[1][3];
  m[1][2] = A[2][3];
  m[2][2] = A[3][3];
  sum += A[0][0] * _det3(m);

  // a01 column
  m[0][0] = A[1][0];
  m[1][0] = A[2][0];
  m[2][0] = A[3][0];
  m[0][1] = A[1][2];
  m[1][1] = A[2][2];
  m[2][1] = A[3][2];
  m[0][2] = A[1][3];
  m[1][2] = A[2][3];
  m[2][2] = A[3][3];
  sum -= A[0][1] * _det3(m);

  // a02 column
  m[0][0] = A[1][0];
  m[1][0] = A[2][0];
  m[2][0] = A[3][0];
  m[0][1] = A[1][1];
  m[1][1] = A[2][1];
  m[2][1] = A[3][1];
  m[0][2] = A[1][3];
  m[1][2] = A[2][3];
  m[2][2] = A[3][3];
  sum += A[0][2] * _det3(m);

  // a03 column
  m[0][0] = A[1][0];
  m[1][0] = A[2][0];
  m[2][0] = A[3][0];
  m[0][1] = A[1][1];
  m[1][1] = A[2][1];
  m[2][1] = A[3][1];
  m[0][2] = A[1][2];
  m[1][2] = A[2][2];
  m[2][2] = A[3][2];
  sum -= A[0][3] * _det3(m);

  return sum;
}

static inline double m4x4det(matrix4x4_t A) { return _det4(A.entries); }

static inline matrix4x4_t m4x4cofactor(matrix4x4_t A) {
  double cofactor[4][4];
  double m[3][3];

  int factor_r = 1;
  int factor_c = 1;
  int m_r = 0;
  int m_c = 0;
  int alt = 1;

  for (int factor_r = 0; factor_r < 4; factor_r++) {
    for (int factor_c = 0; factor_c < 4; factor_c++) {
      for (int r = 0; r < 4; r++) {
        for (int c = 0; c < 4; c++) {
          if (r == factor_r) {
            m_r--;
            break;
          } else if (c != factor_c) {
            m[m_r][m_c++] = A.entries[r][c];
          }
        }

        /*
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++)
            printf("%lf ", m[i][j]);
          printf("\n");
        }
        printf("\n");
        printf("\n");
       */
        m_r++;
        m_c = 0;
      }
      m_r = 0;
      cofactor[factor_r][factor_c] = _det3(m) * alt;
      alt *= -1;
    }
    alt *= -1;
  }
  return m4x4create(cofactor);
}

static inline matrix4x4_t m4x4transpose(matrix4x4_t A) {
  double copy[4][4] = {
    {A.entries[0][0], A.entries[0][1], A.entries[0][2], A.entries[0][3]},
    {A.entries[1][0], A.entries[1][1], A.entries[1][2], A.entries[1][3]},
    {A.entries[2][0], A.entries[2][1], A.entries[2][2], A.entries[2][3]},
    {A.entries[3][0], A.entries[3][1], A.entries[3][2], A.entries[3][3]},
  };
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      A.entries[i][j] = copy[j][i];

  return A;
}

static inline matrix4x4_t m4x4adjoint(matrix4x4_t A) {
  A = m4x4cofactor(A);
  A = m4x4transpose(A);
  return A;
}

static inline matrix4x4_t m4x4scale(matrix4x4_t A, double s) {
  for (int i = 0; i < 4; i++)
    for (int j = 0; j < 4; j++)
      A.entries[i][j] *= s;
  return A;
}

static inline matrix4x4_t m4x4invert(matrix4x4_t A) {
  double det = m4x4det(A);
  matrix4x4_t adjoint = m4x4adjoint(A);
  adjoint = m4x4scale(adjoint, 1/det);
  return adjoint;
}
