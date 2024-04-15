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

typedef struct v3mesh_t {
  vector3_t **list;
  int len;
  int color;
} v3mesh_t;

SDL_Window *window = NULL;
SDL_Renderer *renderer = NULL;
SDL_Texture *texture = NULL;
SDL_Surface *surface = NULL;
double theta = PI / 2;
double phi = 0;
double sensitivity = 0.4;
vector3_t camera = {0, 0, 0};

int init();
int cleanup();
void pixel(int x, int y, int color);
void present();
void line(int x0, int y0, int x1, int y1, int color);
void orthographicProjection();
void matrixMultiply(int n, double matrix0[][n], double matrix1[][n],
                    double ret[n][n]);
void handleMouse();
void triangle(int x0, int y0, int x1, int y1, int x2, int y2, int color);
int pixelList(int *list, int x0, int y0, int x1, int y1);
void clear(int color);
void v2triangle(vector2_t u, vector2_t v, vector2_t w, int color);
void matrixstr(char *str, int n, double matrix[n][n]);
void v3cube(vector3_t p0, vector3_t p1, vector3_t p2, vector3_t p3,
            vector3_t p4, vector3_t p5, vector3_t p6, vector3_t p7, int color);
vector2_t v2NormalizedToScreen(vector2_t v);
void v2line(vector2_t u, vector2_t v, int color);
void normv2line(vector2_t u, vector2_t v, int color);
vector3_t project(vector3_t v);
void normv2triangle(vector2_t u, vector2_t v, vector2_t w, int color);
matrix4x4_t transformationMatrix();
void v3meshdraw(v3mesh_t mesh);
void v3meshfree(v3mesh_t mesh);
v3mesh_t v3meshfromlist(int len, vector3_t list[len][3], int color);
matrix4x4_t getCameraMatrix();
matrix4x4_t getProjectionMatrix();
matrix4x4_t getRotationZMatrix();
matrix4x4_t getRotationXMatrix();
matrix4x4_t getTranslationMatrix();

int main() {

  if (init())
    return 1;

  double r = 1;

  int x = 0;
  int y = 0;
  int z = 0;

  int quit = 0;
  SDL_Event e;

  vector3_t p0 = {0, 0, 0};
  vector3_t p1 = {1, 0, 0};
  vector3_t p2 = {0, 1, 0};
  vector3_t p3 = {1, 1, 0};
  vector3_t p4 = {0, 0, 1};
  vector3_t p5 = {1, 0, 1};
  vector3_t p6 = {0, 1, 1};
  vector3_t p7 = {1, 1, 1};

  vector3_t list[12][3] = {
      // front
      {p0, p3, p1},
      {p0, p2, p3},

      // bottom
      {p5, p4, p0},
      {p5, p0, p1},

      // back
      {p5, p7, p6},
      {p5, p6, p4},

      // top
      {p2, p6, p7},
      {p2, p7, p3},

      // left
      {p4, p6, p2},
      {p4, p2, p0},

      // right
      {p1, p3, p7},
      {p1, p7, p5},

  };

  v3mesh_t mesh = v3meshfromlist(12, list, 0xFFFFFF);

  int i = 0;
  while (!quit) {
    const unsigned char *keystates = SDL_GetKeyboardState(NULL);
    while (SDL_PollEvent(&e) != 0) {
      if (e.type == SDL_QUIT)
        quit = 1;
    }

    if (keystates[SDL_SCANCODE_W])

      y++;
    if (keystates[SDL_SCANCODE_S])
      y--;
    if (keystates[SDL_SCANCODE_A])
      x++;
    if (keystates[SDL_SCANCODE_D])
      x--;

    handleMouse();

    clear(0x0);

    // v3cube(p0, p1, p2, p3, p4, p5, p6, p7, 0xFFFFFF);
    v3meshdraw(mesh);

    //  m4x4print(A);

    // camera.z -= 0.101;
    // camera.x += 0.111;
    // camera.y -= 0.111;
    //  normv2triangle(u, v, w, 0xFFFFFF);

    // rotate(point1);

    present();

    SDL_Delay(1000 / 60);
  }

  v3meshfree(mesh);
  cleanup();

  SDL_Quit();

  return 0;
}

v3mesh_t v3meshfromlist(int len, vector3_t list[len][3], int color) {
  vector3_t **l = malloc(sizeof(vector3_t *) * len);
  for (int i = 0; i < len; i++)
    l[i] = malloc(sizeof(vector3_t) * 3);

  for (int i = 0; i < len; i++)
    for (int j = 0; j < 3; j++)
      l[i][j] = list[i][j];

  v3mesh_t ret;
  ret.list = l;
  ret.len = len;
  ret.color = color;
  return ret;
}

matrix4x4_t getCameraMatrix() {
  double cameraMatrixEntries[4][4] = {
      {1, 0, 0, camera.x},
      {0, 1, 0, camera.y},
      {0, 0, 1, camera.z},
      {0, 0, 0, 5},
  };
  return m4x4create(cameraMatrixEntries);
}

matrix4x4_t getProjectionMatrix() {
  double zfar = 100;
  double znear = 1;

  double a = (double)SCREEN_HEIGHT / (double)SCREEN_WIDTH;
  double f = 1 / tan(theta / 2);
  double q = zfar / (zfar - znear);
  double projectionMatrixEntries[4][4] = {
      {a * f, 0, 0, 0},
      {0, f, 0, 0},
      {0, 0, q, 1},
      {0, 0, -znear * q, 0},
  };

  return m4x4create(projectionMatrixEntries);
}

matrix4x4_t getRotationZMatrix() {
  double rotationz[4][4] = {
      {cos(phi), sin(phi), 0, 0},
      {-sin(phi), cos(phi), 0, 0},
      {0, 0, 1, 0},
      {0, 0, 0, 1},
  };
  return m4x4create(rotationz);
}

matrix4x4_t getRotationXMatrix() {
  double rotationx[4][4] = {
      {1, 0, 0, 0},
      {0, cos(phi), sin(phi), 0},
      {0, -sin(phi), cos(phi), 0},
      {0, 0, 0, 1},
  };
  return m4x4create(rotationx);
}

matrix4x4_t getTranslationMatrix() {
  double translationMatrixEntries[4][4] = {
      {1, 0, 0, 0},
      {0, 1, 0, 0},
      {0, 0, 1, 2},
      {0, 0, 0, 2},
  };
  return m4x4create(translationMatrixEntries);
}

void v3meshdraw(v3mesh_t mesh) {

  phi += 0.01;

  matrix4x4_t cameraMatrix = getCameraMatrix();
  matrix4x4_t projectionMatrix = getProjectionMatrix();
  matrix4x4_t rotationzMatrix = getRotationZMatrix();
  matrix4x4_t rotationxMatrix = getRotationXMatrix();
  matrix4x4_t translationMatrix = getTranslationMatrix();

  matrix4x4_t transform = m4x4mul(4, rotationzMatrix, rotationxMatrix,
                                  translationMatrix, projectionMatrix);

  for (int i = 0; i < mesh.len; i++) {
    vector4_t p0 = v3tov4(mesh.list[i][0]);
    vector4_t p1 = v3tov4(mesh.list[i][1]);
    vector4_t p2 = v3tov4(mesh.list[i][2]);

    p0 = v4m4x4mul(p0, transform);
    p1 = v4m4x4mul(p1, transform);
    p2 = v4m4x4mul(p2, transform);

    vector3_t line1 = SUB(p1, p0);
    vector3_t line2 = SUB(p2, p0);
    vector3_t norm = CROSS(line1, line2);
    norm = NORM(norm);

    if (DOT(norm, SUB(p0, camera)) > 0) {
      
      /*
       * NEED TO DO MORE ADVANCED COLOR CALCULATION WITH LUMINENCE VALUE
      */

      vector3_t light = {1000, 0, -1};
      light = NORM(light);
      double luminence = DOT(light, norm);

      double r = (double)((mesh.color & 0xFF0000) >> 16) / 255;
      double g = (double)((mesh.color & 0x00FF00) >> 8) / 255;
      double b = (double)((mesh.color & 0x0000FF)) / 255;

      double M = MAX(MAX(r, g), b);
      double m = MIN(MIN(r, g), b);
      double C = M - m;

      double h = 0;
      if (M == r)
        h = fmod((g - b) / C, 6);
      else if (M == g)
        h = (b - r) / C + 2;
      else if (M == b)
        h = (r - g) / C + 4;
      h *= PI / 3;

      double intensity = (r + g + b) / 3;

      // color = r << 16 | g << 8 | b;

      int color = mesh.color;
      color = (color & 0x7E7E7E) << (int)(luminence * 2);

      normv2triangle(v4tov2(p0), v4tov2(p1), v4tov2(p2), color);
    }
  }
}

void v3meshfree(v3mesh_t mesh) {
  for (int i = 0; i < mesh.len; i++)
    free(mesh.list[i]);
  free(mesh.list);
}

vector3_t project(vector3_t v) { return v; }

void handleMouse() {
  int mouseX;
  int mouseY;
  SDL_GetMouseState(&mouseX, &mouseY);
  int offsetX = mouseX - SCREEN_WIDTH / 2;
  int offsetY = mouseY - SCREEN_HEIGHT / 2;
  // SDL_WarpMouseInWindow(window, SCREEN_WIDTH / 2, SCREEN_HEIGHT / 2);
  //  theta = fmod(theta + offsetX * sensitivity / 100, (double)2 * PI);
  // if ((fabs(phi + offsetY * sensitivity / 100) <= 2 * PI))
  //   phi = phi + offsetY * sensitivity / 100;
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

#ifdef DEBUG
  line(x0, y0, x1, y1, 0);
  line(x1, y1, x2, y2, 0);
  line(x0, y0, x2, y2, 0);
#endif
}

void clear(int color) { SDL_FillRect(surface, NULL, color); }

void present() {
  texture = SDL_CreateTextureFromSurface(renderer, surface);
  SDL_RenderCopy(renderer, texture, NULL, NULL);
  SDL_RenderPresent(renderer);
  SDL_DestroyTexture(texture);
}
