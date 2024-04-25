#pragma once

#include "vector.h"
#include <stdlib.h>

typedef struct v4queue_t {
  vector4_t *list;
  int head;
  int tail;
  int capacity;
  int count;
} v4queue_t;

static inline v4queue_t *v4qcreate(int capacity) {
  v4queue_t *ret = (v4queue_t *)malloc(sizeof(v4queue_t));
  ret->list = (vector4_t *)malloc(sizeof(vector4_t) * capacity);
  ret->capacity = capacity;
  ret->head = 0;
  ret->tail = 0;
  ret->count = 0;
  return ret;
}

static inline void v4qfree(v4queue_t *q) {
  free(q->list);
  free(q);
}

static inline int v4qisempty(v4queue_t *q) {
  if (q->count == 0)
    return 1;
  return 0;
}

static inline int v4qisfull(v4queue_t *q) {
  if (q->count == q->capacity)
    return 1;
  return 0;
}

static inline int v4enq(v4queue_t *q, vector4_t item) {
  if (v4qisfull(q))
    return 1;
  q->list[q->head++] = item;
  q->head = (q->head + 1) % q->capacity;
  q->count++;
  // need to add overflow protection
  //
  // and realloc possibly
  return 0;
}

static inline int v4deq(v4queue_t *q, vector4_t *ret) {
  if (v4qisempty(q))
    return 1;
  *ret = q->list[q->tail];
  q->count--;
  q->tail = (q->tail + 1) % q->capacity;
  return 0;
}

static inline vector3_t v3intersectplane(vector3_t planepoint, vector3_t planenormal,
                                   vector3_t linestart, vector3_t lineend) {
  planenormal = NORM(planenormal);
  double plane_d = -DOT(planenormal, planepoint); 
  double ad = DOT(linestart, planenormal);
  double bd = DOT(lineend, planenormal);
  double t = (-plane_d - ad)/(bd-ad);
  vector3_t linestarttoend = SUB(lineend, linestart);
  vector3_t intersect = MUL(linestarttoend, t);
  return ADD(linestart, intersect);
}
