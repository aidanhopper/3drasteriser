#include "queue.h"
#include "vector.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_events.h>
#include <SDL2/SDL_surface.h>
#include <math.h>

#define SCREEN_WIDTH 1920
#define SCREEN_HEIGHT 1080
#define PI 3.1415

#define MIN(x, y) x > y ? y : x
#define MAX(x, y) x < y ? y : x

#define DEBUG
// #define WIREFRAMEONLY

typedef struct v3mesh_t {
  vector3_t *vertices;
  int **triangles;
  int vlen;
  int tlen;
  int color;
} v3mesh_t;

SDL_Window *window = NULL;
SDL_Renderer *renderer = NULL;
SDL_Texture *texture = NULL;
SDL_Surface *surface = NULL;
float theta = PI / 4;
float phi = 0;
float alpha = 0;
float sensitivity = 1;
vector3_t camera = {0, 0, 0};
vector3_t forward = {0, 0, -1};
vector3_t up = {0, 1, 0};
vector3_t right = {1, 0, 0};
unsigned int deltatime;
float near = 1;
float far = 100;

float aspectratio = (float)SCREEN_HEIGHT / (float)(SCREEN_WIDTH);

int init();
int cleanup();
void pixel(int x, int y, int color);
void present();
void line(int x0, int y0, int x1, int y1, int color);
void orthographicProjection();
void matrixMultiply(int n, float matrix0[][n], float matrix1[][n],
                    float ret[n][n]);
void handleMouse();
void triangle(int x0, int y0, int x1, int y1, int x2, int y2, int color);
int pixelList(int *list, int x0, int y0, int x1, int y1);
void clear(int color);
void v2triangle(vector2_t u, vector2_t v, vector2_t w, int color);
void matrixstr(char *str, int n, float matrix[n][n]);
vector2_t v2NormalizedToScreen(vector2_t v);
void v2line(vector2_t u, vector2_t v, int color);
void normv2line(vector2_t u, vector2_t v, int color);
vector3_t project(vector3_t v);
void normv2triangle(vector2_t u, vector2_t v, vector2_t w, int color);
void v3meshdraw(v3mesh_t mesh);
void v3meshfree(v3mesh_t mesh);
v3mesh_t v3meshfromlist(int tlen, int vlen, int tris[tlen][3],
                        vector3_t verts[vlen], int color);
v3mesh_t loadobjfile(char *path);
void handleKeyboard();

int main() {

  if (init())
    return 1;

  float r = 1;

  int x = 0;
  int y = 0;
  int z = 0;

  int quit = 0;
  SDL_Event e;

  v3mesh_t mesh = loadobjfile("teapot.obj");
  mesh.color = 0xFFFFFF;

  unsigned int starttime = SDL_GetTicks();
  unsigned int currenttime = SDL_GetTicks();
  deltatime = SDL_GetTicks() - currenttime;
  float elapsedtime = (float)(currenttime - starttime) / 1000;

  while (!quit) {

    while (SDL_PollEvent(&e) != 0) {
      if (e.type == SDL_QUIT)
        quit = 1;
    }
    deltatime = SDL_GetTicks() - currenttime;
    currenttime = SDL_GetTicks();
    elapsedtime = (float)(currenttime - starttime) / 1000;

    handleKeyboard();

    handleMouse();

    clear(0x0);

    v3meshdraw(mesh);

    present();

    // SDL_Delay(1000 / 60);
  }

  v3meshfree(mesh);
  cleanup();

  SDL_Quit();

  return 0;
}

void handleKeyboard() {

  const unsigned char *keystates = SDL_GetKeyboardState(NULL);

  float speed = 0.01;
  vector3_t velocity = {0, 0, 0};

  if (keystates[SDL_SCANCODE_UP]) {
    velocity.y += -speed * deltatime;
  }
  if (keystates[SDL_SCANCODE_DOWN]) {
    velocity.y += speed * deltatime;
  }
  if (keystates[SDL_SCANCODE_A]) {
    vector3_t velocityscaled = MUL(right, speed * deltatime);
    velocity = ADD(velocity, velocityscaled);
  }
  if (keystates[SDL_SCANCODE_D]) {
    vector3_t velocityscaled = MUL(right, -speed * deltatime);
    velocity = ADD(velocity, velocityscaled);
  }
  if (keystates[SDL_SCANCODE_W]) {
    vector3_t velocityscaled = MUL(forward, -speed * deltatime);
    velocity = ADD(velocity, velocityscaled);
    // velocity.z -= 1;
  }
  if (keystates[SDL_SCANCODE_S]) {
    vector3_t velocityscaled = MUL(forward, speed * deltatime);
    velocity = ADD(velocity, velocityscaled);
    // velocity.z += 1;
  }

  camera = ADD(camera, velocity);

  if (keystates[SDL_SCANCODE_LEFT])
    phi += 0.002 * deltatime;
  if (keystates[SDL_SCANCODE_RIGHT])
    phi += -0.002 * deltatime;
}

v3mesh_t loadobjfile(char *path) {
  v3mesh_t ret;
  int tlen = 0;
  int vlen = 0;
  size_t filelen = 0;
  int read;

  char *line;
  FILE *f = fopen(path, "r");
  if (f == NULL)
    return ret;

  while ((read = getline(&line, &filelen, f)) != -1) {
    if (line[0] == 'v')
      vlen++;
    else if (line[0] == 'f')
      tlen++;
  }

  fclose(f);
  if (line)
    free(line);

  f = fopen(path, "r");
  filelen = 0;
  read = 0;

  int triangles[tlen][3];
  vector3_t vertices[vlen];
  int v = 0;
  int t = 0;
  while ((read = getline(&line, &filelen, f)) != -1) {
    char *tok = strtok(line, " ");
    if (tok[0] == 'v') {
      float x, y, z;
      tok = strtok(NULL, " ");
      sscanf(tok, "%f", &x);
      tok = strtok(NULL, " ");
      sscanf(tok, "%f", &y);
      tok = strtok(NULL, " ");
      sscanf(tok, "%f", &z);
      vertices[v++] = (vector3_t){x, y, z};
    } else if (tok[0] == 'f') {
      int p0, p1, p2;
      tok = strtok(NULL, " ");
      sscanf(tok, "%d", &p0);
      tok = strtok(NULL, " ");
      sscanf(tok, "%d", &p1);
      tok = strtok(NULL, " ");
      sscanf(tok, "%d", &p2);
      triangles[t][0] = p0 - 1;
      triangles[t][1] = p1 - 1;
      triangles[t][2] = p2 - 1;
      t++;
    }
  }

  fclose(f);
  if (line)
    free(line);

  ret = v3meshfromlist(tlen, vlen, triangles, vertices, 0xFFFFFF);

  return ret;
}

v3mesh_t v3meshfromlist(int tlen, int vlen, int tris[tlen][3],
                        vector3_t verts[vlen], int color) {

  int **triangles = malloc(sizeof(int *) * tlen);
  for (int i = 0; i < tlen; i++)
    triangles[i] = malloc(sizeof(int) * 3);

  for (int i = 0; i < tlen; i++)
    for (int j = 0; j < 3; j++)
      triangles[i][j] = tris[i][j];

  vector3_t *vertices = malloc(sizeof(vector3_t) * vlen);
  for (int i = 0; i < vlen; i++)
    vertices[i] = verts[i];

  v3mesh_t ret;
  ret.triangles = triangles;
  ret.tlen = tlen;
  ret.vertices = vertices;
  ret.vlen = vlen;
  ret.color = color;
  return ret;
}

void hextohsv(int color, float *hue, float *saturation, float *value) {
  float r = (float)((color & 0xFF0000) >> 16) / 255;
  float g = (float)((color & 0x00FF00) >> 8) / 255;
  float b = (float)((color & 0x0000FF)) / 255;

  float cmax = MAX(MAX(r, g), b);
  float cmin = MIN(MIN(r, g), b);
  float range = cmax - cmin;

  *hue = 0;
  if (range != 0) {
    if (cmax == r)
      *hue = 60 * ((g - b) / range);
    else if (cmax == g)
      *hue = 60 * ((b - r) / range + 2);
    else if (cmax == b)
      *hue = 60 * ((r - g) / range + 4);
  }
  if (*hue < 0)
    *hue += 360;

  *saturation = 0;
  if (cmax != 0)
    *saturation = (range / cmax);

  *value = cmax;
}

int hsvtohex(float hue, float saturation, float value) {
  float H = hue, S = saturation, V = value, P, Q, T, fract;
  float r, g, b;

  (H == 360.) ? (H = 0.) : (H /= 60.);
  fract = H - floor(H);

  P = V * (1. - S);
  Q = V * (1. - S * fract);
  T = V * (1. - S * (1. - fract));

  if (0. <= H && H < 1.) {
    // RGB = (rgb){.r = V, .g = T, .b = P};
    r = V;
    g = T;
    b = P;
  } else if (1. <= H && H < 2.) {
    // RGB = (rgb){.r = Q, .g = V, .b = P};
    r = Q;
    g = V;
    b = P;
  } else if (2. <= H && H < 3.) {
    // RGB = (rgb){.r = P, .g = V, .b = T};
    r = P;
    g = V;
    b = T;
  } else if (3. <= H && H < 4.) {
    // RGB = (rgb){.r = P, .g = Q, .b = V};
    r = P;
    g = Q;
    b = V;
  } else if (4. <= H && H < 5.) {
    // RGB = (rgb){.r = T, .g = P, .b = V};
    r = T;
    g = P;
    b = V;
  } else if (5. <= H && H < 6.) {
    // RGB = (rgb){.r = V, .g = P, .b = Q};
    r = V;
    g = P;
    b = Q;
  } else {
    // RGB = (rgb){.r = 0., .g = 0., .b = 0.};
    r = 0;
    g = 0;
    b = 0;
  }

  r *= 255;
  g *= 255;
  b *= 255;

  int color = (int)r << 16 | (int)g << 8 | (int)b;

  return color;
};

int v4cmpz(const void *u, const void *v) {
  vector4_t *v4u = (vector4_t *)u;
  vector4_t *v4v = (vector4_t *)v;

  float z00 = v4u[0].z;
  float z01 = v4u[1].z;
  float z02 = v4u[2].z;

  float z10 = v4v[0].z;
  float z11 = v4v[1].z;
  float z12 = v4v[2].z;

  float z0mid = (z00 + z01 + z02) / 3;
  float z1mid = (z10 + z11 + z12) / 3;

  if (z0mid < z1mid)
    return -1;
  if (z0mid > z1mid)
    return 1;
  return 0;
}

vector3_t getLeftPlanePoint(vector4_t p) { return (vector3_t){-p.w, 0, p.z}; }
vector3_t getLeftPlaneNorm(vector4_t p) {
  vector3_t absoluteup = {0, 1, 0};
  vector3_t leftp = CROSS(absoluteup, getLeftPlanePoint(p));
  return NORM(leftp);
}

vector3_t getTopPlanePoint(vector4_t p) { return (vector3_t){0, p.w, p.z}; }
vector3_t getTopPlaneNorm(vector4_t p) {
  vector3_t absoluteright = {1, 0, 0};
  vector3_t topp = CROSS(getTopPlanePoint(p), absoluteright);
  return NORM(topp);
}

vector3_t getBottomPlanePoint(vector4_t p) { return (vector3_t){0, -p.w, p.z}; }
vector3_t getBottomPlaneNorm(vector4_t p) {
  vector3_t absoluteright = {1, 0, 0};
  vector3_t bottomp = CROSS(absoluteright, getBottomPlanePoint(p));
  return NORM(bottomp);
}

vector3_t getNearPlanePoint() { return (vector3_t){0, 0, -near}; }
vector3_t getNearPlaneNorm() { return NORM(getNearPlanePoint()); }

vector3_t getFarPlaneNorm() {
  vector3_t farp = {0, 0, -far};
  return NORM(farp);
}

vector4_t intersection(vector4_t q1, vector4_t q2, vector3_t p, vector3_t n) {
  /*
  */

  vector2_t v2q1 = v4tov2(q1);
  vector2_t v2q2 = v4tov2(q2);

  n = getLeftPlaneNorm(q1);

  vector4_t i1;
  i1.w = q2.w;
  i1.x = -i1.w;
  i1.y = q2.y;
  i1.z = q2.z;

  //v4print(q1);
  //v4print(q2);
  //v4print(i1);
  //v2print(v4tov2(i1));

  float d1 = DOT(n, SUB(q1, p));
  //printf("d1: %f\n", d1);
  float d2 = DOT(n, SUB(q2, p));
  //printf("d2: %f\n", d2);
  //printf("d1-d2: %f\n", d1-d2);
  float t = d1 / (d1 - d2);
  //printf("t %f\n", t);
  vector3_t dir = SUB(q2, q1);
  vector4_t i2 = v3tov4(ADD(q1, MUL(dir, t)));
  i2.w = q2.w;
  //i2.w = q2.w;
  //i2.x = -q2.w;

  printf("q1v4 = ");
  v4print(q1);
  printf("q2v4 = ");
  v4print(q2);
  printf("i = ");
  v4print(i2);
  printf("div = ");
  v2print(v4tov2(i2));

  /*
   */

  printf("\n");

  return i2;
}

int clip(vector4_t outputPoints[100]) {
#define TEST(n, q, p) (DOT(n, SUB(q, p)))
#define OUT(w, q) (q <= w)
  vector4_t q0 = outputPoints[0];
  vector4_t q1 = outputPoints[1];
  vector4_t q2 = outputPoints[2];

  vector3_t p = getBottomPlanePoint(q0);
  vector3_t n = getBottomPlaneNorm(q0);

  p = getLeftPlanePoint(q0);
  n = getLeftPlaneNorm(q0);

  vector4_t in[10];
  vector4_t out[10];

  int outlen = 0;
  int inlen = 0;

  if (OUT(-q0.w, q0.x)) {
    out[outlen++] = q0;
  } else {
    in[inlen++] = q0;
  }

  if (OUT(-q1.w, q1.x)) {
    out[outlen++] = q1;
  } else {
    in[inlen++] = q1;
  }

  if (OUT(-q2.w, q2.x)) {
    out[outlen++] = q2;
  } else {
    in[inlen++] = q2;
  }

  if (inlen == 0) {
    return 0;
  }

  normv2line(v4tov2(q0), v4tov2(q1), 0x00FFFF);
  normv2line(v4tov2(q0), v4tov2(q2), 0x00FFFF);
  normv2line(v4tov2(q1), v4tov2(q2), 0x00FFFF);
  if (inlen == 1) {
    vector4_t I0 = intersection(in[0], out[0], p, n);
    vector4_t I1 = intersection(in[0], out[1], p, n);

    // v2print(v4tov2(I0));


    //normv2triangle(v4tov2(in[0]), v4tov2(I0), v4tov2(I1), 0x00FFFF);
    //normv2line(v4tov2(in[0]), v4tov2(I1), 0xFF0000);
    //normv2line(v4tov2(in[0]), v4tov2(I0), 0xFF0000);
    //normv2line(v4tov2(I0), v4tov2(I1), 0xFF0000);
  }

  if (inlen == 2) {
    vector4_t I0 = intersection(in[0], out[0], p, n);
    vector4_t I1 = intersection(in[1], out[0], p, n);

    //normv2triangle(v4tov2(in[0]), v4tov2(I0), v4tov2(in[1]), 0x00FFFF);
    //normv2triangle(v4tov2(in[1]), v4tov2(I0), v4tov2(I1), 0x00FFFF);

    //normv2line(v4tov2(in[0]), v4tov2(I0), 0xFF0000);
    //normv2line(v4tov2(in[0]), v4tov2(in[1]), 0xFF0000);
    //normv2line(v4tov2(I0), v4tov2(in[1]), 0xFF0000);

    //normv2line(v4tov2(in[1]), v4tov2(I0), 0xFF0000);
    //normv2line(v4tov2(in[1]), v4tov2(I0), 0xFF0000);
    //normv2line(v4tov2(I0), v4tov2(I1), 0xFF0000);
  }

  return 0;
  // printf("%d\n", outlen);

  /*
  float epsilon = 0;

  int outputPointsLen = 3;

  vector4_t tmp = outputPoints[0];

  vector3_t planePoints[] = {getLeftPlanePoint(tmp), getRightPlanePoint(tmp),
                             getTopPlanePoint(tmp), getBottomPlanePoint(tmp),
                             getNearPlanePoint()};

  vector3_t planeNormals[] = {
      getLeftPlaneNorm(tmp),   getRightPlaneNorm(tmp), getTopPlaneNorm(tmp),
      getBottomPlaneNorm(tmp), getNearPlaneNorm(),
  };

  // for each clip edge in clip polygon
  for (int i = 3; i < 4; i++) {

    // get clip plane normal
    vector3_t planeNormal = planeNormals[i];
    vector3_t planePoint = planePoints[i];

    // copy output points to input points
    vector4_t inputPoints[100];
    int inputPointsLen = 0;
    for (int j = 0; j < outputPointsLen; j++)
      inputPoints[inputPointsLen++] = outputPoints[j];

    // clear output points
    outputPointsLen = 0;

    // clip every point in input list
    for (int j = 0; j < inputPointsLen; j++) {

      // get points
      vector4_t currentPoint = inputPoints[j];
      vector4_t prevPoint = inputPoints[(j - 1) % inputPointsLen];

      // normv2line(v4tov2(intersectionPoint), v4tov2(prevPoint), 0xFF0000);
      vector4_t intersectionPoint =
          intersection(currentPoint, prevPoint, planePoint, planeNormal);

      // if current point is inside
      if (TEST(planeNormal, currentPoint, planePoint) > 0) {

        // prev point is outside
        if (TEST(planeNormal, prevPoint, planePoint) < 0) {
          outputPoints[outputPointsLen++] = intersectionPoint;
        }
        outputPoints[outputPointsLen++] = currentPoint;

      }

      else if (TEST(planeNormal, prevPoint, planePoint) > 0) {
        outputPoints[outputPointsLen++] = intersectionPoint;
      }
    }
  }
  return outputPointsLen;
  */
}

void v3meshdraw(v3mesh_t mesh) {
  matrix4x4_t camerarot = createRotationMatrix(0, phi, 0);

  // I need figure out the exact operations I should be performing
  // because my ordering might be wrong
  matrix4x4_t rot1 = createRotationMatrix(0, phi, 0);

  // matrix4x4_t rot2 = createRotationMatrix(0, phi/2, 0);
  matrix4x4_t translation =
      createTranslationMatrix(-camera.x, -camera.y, -camera.z + 10);
  // matrix4x4_t view = createViewMatrix();
  matrix4x4_t projection =
      createProjectionMatrix(near, far, aspectratio, theta);

  forward = (vector3_t){0, 0, -1};
  right = (vector3_t){1, 0, 0};

  vector3_t Y = {0, 1, 0};
  vector3_t center = {0, 0, -1};
  // target= v3m4x4mul(target, camerarot);

  forward = v3m4x4mul(forward, camerarot);
  right = v3m4x4mul(right, camerarot);

  right = CROSS(forward, Y);
  right = NORM(right);
  up = CROSS(right, forward);

  vector3_t U = CROSS(forward, up);
  vector3_t V = CROSS(U, forward);
  vector3_t W = MUL(forward, -1);

  float lookAtEntries[4][4] = {
      {U.x, V.x, W.x, 0},
      {U.y, V.y, W.y, 0},
      {U.z, V.z, W.z, 0},
      //{-camera.x, -camera.y, -camera.z, 1},
      {0, 0, 0, 1},

  };
  matrix4x4_t lookat = m4x4create(lookAtEntries);
  matrix4x4_t pointat = m4x4invert(lookat);

  matrix4x4_t dir = createRotationMatrix(0, -phi, 0);

  forward = (vector3_t){0, 0, -1};
  forward = v3m4x4mul(forward, dir);

  right = (vector3_t){1, 0, 0};
  right = v3m4x4mul(right, dir);

  matrix4x4_t view = m4x4invert(m4x4mul(2, lookat, translation));

  vector4_t draworder[mesh.tlen][3];
  int drawcount = 0;

  for (int i = 0; i < mesh.tlen; i++) {
    vector4_t p0 = v3tov4(mesh.vertices[mesh.triangles[i][0]]);
    vector4_t p1 = v3tov4(mesh.vertices[mesh.triangles[i][1]]);
    vector4_t p2 = v3tov4(mesh.vertices[mesh.triangles[i][2]]);

    // vector4_t p0projected = v4m4x4mul(p0, transform);
    // vector4_t p1projected = v4m4x4mul(p1, transform);
    // vector4_t p2projected = v4m4x4mul(p2, transform);

    vector4_t p0viewed = v4m4x4mul(p0, view);
    vector4_t p1viewed = v4m4x4mul(p1, view);
    vector4_t p2viewed = v4m4x4mul(p2, view);
    // v4print(p0viewed);

    vector3_t surfacenorm = SURFACENORM(p0viewed, p1viewed, p2viewed);

    if (DOT(surfacenorm, p0viewed) < 0) {

      draworder[drawcount][0] = p0viewed;
      draworder[drawcount][1] = p1viewed;
      draworder[drawcount][2] = p2viewed;
      drawcount++;
    }
  }

  // printf("\n");

  // sort points from largest z to smallest z
  qsort(draworder, drawcount, sizeof(vector4_t) * 3, v4cmpz);

  // uses painters algorithm for drawing
  for (int i = 0; i < drawcount; i++) {
    vector4_t p0projected = v4m4x4mul(draworder[i][0], projection);
    vector4_t p1projected = v4m4x4mul(draworder[i][1], projection);
    vector4_t p2projected = v4m4x4mul(draworder[i][2], projection);

    // does crude clipping
    if (-p0projected.w < p0projected.x && p0projected.x < p0projected.w &&
        -p0projected.w < p0projected.y && p0projected.y < p0projected.w &&
        -p1projected.w < p1projected.x && p1projected.x < p1projected.w &&
        -p1projected.w < p1projected.y && p1projected.y < p1projected.w &&
        -p2projected.w < p2projected.x && p2projected.x < p2projected.w &&
        -p2projected.w < p2projected.y && p2projected.y < p2projected.w) {

      // calculate normal for surface
      vector3_t surfacenorm =
          SURFACENORM(p0projected, p1projected, p2projected);

      // color based on normal
      vector3_t light = {0, 0, 1};
      vector3_t lightnorm = NORM(light);
      float luminence = DOT(lightnorm, surfacenorm) / LEN(light);

      float hue, saturation, value;
      hextohsv(mesh.color, &hue, &saturation, &value);

      value *= MAX(luminence, 0.4);

      int color = hsvtohex(hue, saturation, value);

      vector2_t p0image = v4tov2(p0projected);
      vector2_t p1image = v4tov2(p1projected);
      vector2_t p2image = v4tov2(p2projected);
      normv2triangle(p0image, p1image, p2image, color);
    }

    else {


      vector4_t outputPoints[100];
      outputPoints[0] = p0projected;
      outputPoints[1] = p1projected;
      outputPoints[2] = p2projected;
      clip(outputPoints);


    }
  }
}

void v3meshfree(v3mesh_t mesh) {
  for (int i = 0; i < mesh.tlen; i++)
    free(mesh.triangles[i]);
  free(mesh.triangles);
  free(mesh.vertices);
}

vector3_t project(vector3_t v) { return v; }

void handleMouse() {
  int mouseX;
  int mouseY;
  SDL_GetMouseState(&mouseX, &mouseY);
  int offsetX = mouseX - SCREEN_WIDTH / 2;
  int offsetY = mouseY - SCREEN_HEIGHT / 2;
  // SDL_WarpMouseInWindow(window, SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);
  //  theta = fmod(theta + offsetX * sensitivity / 100, (float)2 * PI);
  /*
  if (offsetX > 2)
    offsetX = 2;
  if (offsetX < -2)
    offsetX = -2;
  phi += ((float)offsetX / 100) * sensitivity;

  if (offsetY > 2)
    offsetY = 2;
  if (offsetY < -2)
    offsetY = -2;
  alpha += ((float)offsetY / 100) * sensitivity;

  */
  //  if ((fabs(phi + offsetY * sensitivity / 100) <= 2 * PI))
  //    phi = phi + offsetY * sensitivity / 100;
}

int init() {
  if (SDL_Init(SDL_INIT_VIDEO) < 0)
    return 1;

  window = SDL_CreateWindow("Framebuffer Example", SDL_WINDOWPOS_UNDEFINED,
                            SDL_WINDOWPOS_UNDEFINED, SCREEN_WIDTH,
                            SCREEN_HEIGHT, SDL_WINDOW_SHOWN);
  if (window == NULL)
    return 1;

  renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);
  if (renderer == NULL)
    return 1;

  surface =
      SDL_CreateRGBSurface(0, SCREEN_WIDTH, SCREEN_HEIGHT, 32, 0, 0, 0, 0);
  if (surface == NULL)
    return 1;

  return 0;
}

int cleanup() {
  SDL_FreeSurface(surface);
  SDL_DestroyRenderer(renderer);
  SDL_DestroyWindow(window);
  return 0;
}

void pixel(int x, int y, int color) {
  if (x >= SCREEN_WIDTH || y >= SCREEN_HEIGHT || x < 0 || y < 0)
    return;
  int *pixels = (int *)surface->pixels;
  pixels[y * surface->w + x] = color;
}

void normv2line(vector2_t u, vector2_t v, int color) {
  vector2_t u_ = v2NormalizedToScreen(u);
  vector2_t v_ = v2NormalizedToScreen(v);
  v2line(u_, v_, color);
}

void normv2triangle(vector2_t u, vector2_t v, vector2_t w, int color) {
  u = v2NormalizedToScreen(u);
  v = v2NormalizedToScreen(v);
  w = v2NormalizedToScreen(w);
  v2triangle(u, v, w, color);
}

void v2line(vector2_t u, vector2_t v, int color) {
  line(u.x, u.y, v.x, v.y, color);
}

void line(int x0, int y0, int x1, int y1, int color) {

  x0 = MAX(x0, 0);
  y0 = MAX(y0, 0);
  x1 = MAX(x1, 0);
  y1 = MAX(y1, 0);

  x0 = MIN(x0, SCREEN_WIDTH);
  y0 = MIN(y0, SCREEN_HEIGHT);
  x1 = MIN(x1, SCREEN_WIDTH);
  y1 = MIN(y1, SCREEN_HEIGHT);

  int dx = abs(x1 - x0);
  int sx = x0 < x1 ? 1 : -1;
  int dy = -abs(y1 - y0);
  int sy = y0 < y1 ? 1 : -1;
  int error = dx + dy;

  while (1) {
    pixel(x0, y0, color);
    if (x0 == x1 && y0 == y1)
      break;
    int e2 = 2 * error;
    if (e2 >= dy) {
      if (x0 == x1)
        break;
      error = error + dy;
      x0 = x0 + sx;
    }
    if (e2 <= dx) {
      if (y0 == y1)
        break;
      error = error + dx;
      y0 = y0 + sy;
    }
  }
}

int pixelList(int *list, int x0, int y0, int x1, int y1) {

  int dx = abs(x1 - x0);
  int sx = x0 < x1 ? 1 : -1;
  int dy = -abs(y1 - y0);
  int sy = y0 < y1 ? 1 : -1;
  int error = dx + dy;

  while (1) {
    if (y0 < SCREEN_HEIGHT && y0 >= 0)
      list[y0] = x0;
    if (x0 == x1 && y0 == y1)
      break;
    int e2 = 2 * error;
    if (e2 >= dy) {
      if (x0 == x1)
        break;
      error = error + dy;
      x0 = x0 + sx;
    }
    if (e2 <= dx) {
      if (y0 == y1)
        break;
      error = error + dx;
      y0 = y0 + sy;
    }
  }

  return 0;
}

void v2triangle(vector2_t u, vector2_t v, vector2_t w, int color) {
  triangle(u.x, u.y, v.x, v.y, w.x, w.y, color);
}

vector2_t v2NormalizedToScreen(vector2_t v) {
  return (vector2_t){
      (SCREEN_WIDTH + SCREEN_WIDTH * v.x) / 2,
      (SCREEN_HEIGHT - SCREEN_HEIGHT * v.y) / 2,
  };
}

void triangle(int x0, int y0, int x1, int y1, int x2, int y2, int color) {

#ifndef WIREFRAMEONLY
  if (y1 > y0) {
    triangle(x1, y1, x0, y0, x2, y2, color);
    return;
  }

  int linelist[SCREEN_HEIGHT];
  for (int i = 0; i < SCREEN_HEIGHT; i++)
    linelist[i] = -1;

  int linelist2[SCREEN_HEIGHT];
  for (int i = 0; i < SCREEN_HEIGHT; i++)
    linelist2[i] = -1;

  // need to draw 0 to 1 line and 0 to 2 line
  // and use linelist to check if im out of bounds

  int x0t = x0;
  int y0t = y0;
  int x1t = x1;
  int y1t = y1;
  int x2t = x2;
  int y2t = y2;

  pixelList(linelist, x1, y1, x2, y2);
  pixelList(linelist2, x0, y0, x2, y2);

  int dx = abs(x1 - x0);
  int sx = x0 < x1 ? 1 : -1;
  int dy = -abs(y1 - y0);
  int sy = y0 < y1 ? 1 : -1;
  int error = dx + dy;

  // draws from x0 to x2 line
  while (1) {
    if (linelist[y0] != -1)
      line(x0, y0, linelist[y0], y0, color);
    if (linelist2[y0] != -1)
      line(x0, y0, linelist2[y0], y0, color);
    if (x0 == x1 && y0 == y1)
      break;
    int e2 = 2 * error;
    if (e2 >= dy) {
      if (x0 == x1)
        break;
      error = error + dy;
      x0 = x0 + sx;
    }
    if (e2 <= dx) {
      if (y0 == y1)
        break;
      error = error + dx;
      y0 = y0 + sy;
    }
  }

  x0 = x0t;
  y0 = y0t;
  x1 = x2t;
  y1 = y2t;

  dx = abs(x1 - x0);
  sx = x0 < x1 ? 1 : -1;
  dy = -abs(y1 - y0);
  sy = y0 < y1 ? 1 : -1;
  error = dx + dy;

  // draws from x0 to x2 line
  while (1) {
    if (linelist[y0] != -1)
      line(x0, y0, linelist[y0], y0, color);
    if (linelist2[y0] != -1)
      line(x0, y0, linelist2[y0], y0, color);
    if (x0 == x1 && y0 == y1)
      break;
    int e2 = 2 * error;
    if (e2 >= dy) {
      if (x0 == x1)
        break;
      error = error + dy;
      x0 = x0 + sx;
    }
    if (e2 <= dx) {
      if (y0 == y1)
        break;
      error = error + dx;
      y0 = y0 + sy;
    }
  }

  x0 = x0t;
  x1 = x1t;
  x2 = x2t;
  y0 = y0t;
  y1 = y1t;
  y2 = y2t;

#endif
#ifdef DEBUG
  line(x0, y0, x1, y1, 0xFFFFFF);
  line(x1, y1, x2, y2, 0xFFFFFF);
  line(x0, y0, x2, y2, 0xFFFFFF);
#endif
}

void clear(int color) { SDL_FillRect(surface, NULL, color); }

void present() {
  texture = SDL_CreateTextureFromSurface(renderer, surface);
  SDL_RenderCopy(renderer, texture, NULL, NULL);
  SDL_RenderPresent(renderer);
  SDL_DestroyTexture(texture);
}
