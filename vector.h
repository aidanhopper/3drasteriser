#pragma once
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

#define DOT(v, u) v.x *u.x + v.y *u.y + v.z *u.z
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
  return v;
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
  if (p.z != 0) {
    ret.x /= p.z;
    ret.y /= p.z;
  }
  return ret;
}
