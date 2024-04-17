#pragma once

#include "vector.h"
#include <stdlib.h>

typedef struct v4queue_t {
  vector4_t *list;
  int head;
  int tail;
  int capacity;
  char empty;
} v4queue_t;

static inline v4queue_t *v4qcreate(int capacity) {
  v4queue_t *ret = (v4queue_t *)malloc(sizeof(v4queue_t));
  ret->list = (vector4_t *)malloc(sizeof(vector4_t) * capacity);
  ret->capacity = capacity;
  ret->head = 0;
  ret->tail = 0;
  ret->empty = 1;
  return ret;
}

static inline void v4qfree(v4queue_t *q) {
  free(q->list);
  free(q);
}

static inline int v4qisempty(v4queue_t *q) { 
  if (q->empty)
    return 1; 
  return 0;
}

static inline int v4qisfull(v4queue_t *q) { 
  if (q->tail == q->head && !q->empty)
    return 1;
  return 0; 
}

static inline int v4enq(v4queue_t *q, vector4_t item) {
  if (v4qisfull(q))
    return 1;
  q->list[q->head++] = item;
  q->head = (q->head + 1) % q->capacity;
  q->empty = 0;
  // need to add overflow protection
  //
  // and realloc possibly
  return 0;
}

static inline int v4deq(v4queue_t *q, vector4_t *ret) {
  if (v4qisempty(q))
    return 1;
  *ret = q->list[q->tail];
  if (q->tail == q->head)
    q->empty = 1;
  q->tail = (q->tail + 1) % q->capacity;
  return 0;
}
