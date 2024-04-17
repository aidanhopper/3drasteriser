#pragma once

#include "vector.h"
#include <stdlib.h>

typedef struct v4queue_t {
  vector4_t *list;
  int head;
  int tail;
  int capacity;
} v4queue_t;

static inline v4queue_t *v4qcreate(int capacity) {
  v4queue_t *ret = (v4queue_t *)malloc(sizeof(v4queue_t));
  ret->list = (vector4_t *)malloc(sizeof(vector4_t) * ret->capacity);
  ret->capacity = capacity;
  ret->head = 0;
  ret->tail = 0;
  return ret;
}

static inline void v4qfree(v4queue_t *q) {
  free(q->list);
  free(q);
}

static inline void v4enqueue(v4queue_t *q, vector4_t item) {
  q->list[q->head] = item;
  q->head = (q->head + 1) % q->capacity;
  // need to add overflow protection
  // and realloc possibly
}

static inline vector4_t v4dequeue(v4queue_t *q) {
  vector4_t ret = q->list[q->tail];
  q->tail = (q->tail + 1) % q->capacity;
  return ret;
}

static inline int v4qisempty() { return 0; }

static inline int v4qisfull() { return 0; }
