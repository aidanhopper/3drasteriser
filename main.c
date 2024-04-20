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

// #define DEBUG

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
double theta = PI / 3;
double phi = 0;
double alpha = 0;
double sensitivity = 1;
vector3_t camera = {0, 0, 0};
vector3_t forward = {0, 0, -1};
vector3_t up = {0, 1, 0};
vector3_t right = {1, 0, 0};

double aspectratio = (double)SCREEN_HEIGHT / (double)(SCREEN_WIDTH);

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

int main() {

  if (init())
    return 1;

  double r = 1;

  int x = 0;
  int y = 0;
  int z = 0;

  int quit = 0;
  SDL_Event e;

  v3mesh_t mesh = loadobjfile("mountains.obj");
  mesh.color = 0xFFFFFF;

  unsigned int starttime = SDL_GetTicks();
  while (!quit) {

    const unsigned char *keystates = SDL_GetKeyboardState(NULL);
    while (SDL_PollEvent(&e) != 0) {
      if (e.type == SDL_QUIT)
        quit = 1;
    }
    unsigned int currenttime = SDL_GetTicks();
    double elapsedtime = (double)(currenttime - starttime) / 1000;
    double speed = 0.7;
    vector3_t velocity = {0, 0, 0};
    if (keystates[SDL_SCANCODE_UP]) {
      velocity.y += speed;
    }
    if (keystates[SDL_SCANCODE_DOWN]) {
      velocity.y -= speed;
    }
    if (keystates[SDL_SCANCODE_A]) {
      vector3_t velocityscaled = MUL(right, speed);
      velocity = ADD(velocity, velocityscaled);
    }
    if (keystates[SDL_SCANCODE_D]) {
      vector3_t velocityscaled = MUL(right, -speed);
      velocity = ADD(velocity, velocityscaled);
    }
    if (keystates[SDL_SCANCODE_W]) {
      vector3_t velocityscaled = MUL(forward, speed);
      velocity = ADD(velocity, velocityscaled);
    }
    if (keystates[SDL_SCANCODE_S]) {
      vector3_t velocityscaled = MUL(forward, -speed);
      velocity = ADD(velocity, velocityscaled);
    }

    camera = ADD(camera, velocity);

    if (keystates[SDL_SCANCODE_LEFT])
      phi += -0.01;
    if (keystates[SDL_SCANCODE_RIGHT])
      phi += 0.01;

    handleMouse();

    clear(0x0);

    vector2_t p0 = {-0.019036, -0.020741};
    vector2_t p1 = {1, -1};
    vector2_t p2 = {1, 1};
    // normv2triangle(p0, p1, p2, 0xFF0000);
    v3meshdraw(mesh);

    present();

    // SDL_Delay(1000 / 60);
  }

  v3meshfree(mesh);
  cleanup();

  SDL_Quit();

  return 0;
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
      double x, y, z;
      tok = strtok(NULL, " ");
      sscanf(tok, "%lf", &x);
      tok = strtok(NULL, " ");
      sscanf(tok, "%lf", &y);
      tok = strtok(NULL, " ");
      sscanf(tok, "%lf", &z);
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

void hextohsv(int color, double *hue, double *saturation, double *value) {
  double r = (double)((color & 0xFF0000) >> 16) / 255;
  double g = (double)((color & 0x00FF00) >> 8) / 255;
  double b = (double)((color & 0x0000FF)) / 255;

  double cmax = MAX(MAX(r, g), b);
  double cmin = MIN(MIN(r, g), b);
  double range = cmax - cmin;

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

int hsvtohex(double hue, double saturation, double value) {
  double H = hue, S = saturation, V = value, P, Q, T, fract;
  double r, g, b;

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

  double z00 = v4u[0].z;
  double z01 = v4u[1].z;
  double z02 = v4u[2].z;

  double z10 = v4v[0].z;
  double z11 = v4v[1].z;
  double z12 = v4v[2].z;

  double z0mid = (z00 + z01 + z02) / 3;
  double z1mid = (z10 + z11 + z12) / 3;

  if (z0mid < z1mid)
    return 1;
  if (z0mid > z1mid)
    return -1;
  return 0;
}

double lerp(double x, vector2_t p0, vector2_t p1) {
  double m =
      ((p1.y - p0.y) / (p1.x - p0.x));

  return m * (x - p0.x) + p0.y;
}

int cliptriangle(vector4_t p0, vector4_t p1, vector4_t p2,
                 vector4_t tribuf[3][3]) {

  int len = 0;
  tribuf[0][0] = p0;
  tribuf[0][1] = p1;
  tribuf[0][2] = p2;

  vector2_t in[5];
  int inlen = 0;

  vector2_t out[5];
  int outlen = 0;

  vector2_t p0projected = {p0.x / p0.w, p0.y / p0.w};
  vector2_t p1projected = {p1.x / p1.w, p1.y / p1.w};
  vector2_t p2projected = {p2.x / p2.w, p2.y / p2.w};

  // for each plane (x=-1, x=1, y=-1, y=1)
  //   clip

  if (p0projected.x < -1) {
    p0projected.y = lerp(-1, p0projected, p1projected);
    in[inlen++] = p0projected;
  }
  else {
    in[inlen++] = p0projected;
  }




  if (-p0.w <= p0.x && p0.x <= p0.w && -p0.w <= p0.y && p0.y <= p0.w &&
      -p1.w <= p1.x && p1.x <= p1.w && -p1.w <= p1.y && p1.y <= p1.w &&
      -p2.w <= p2.x && p2.x <= p2.w && -p2.w <= p2.y && p2.y <= p2.w)
    return 1;

  return 0;
}

void v3meshdraw(v3mesh_t mesh) {

  matrix4x4_t transform2 =
      createTranslationMatrix(-camera.x, -camera.y, -camera.z + 10);

  matrix4x4_t rotation = createRotationMatrix(alpha, phi, 0);

  matrix4x4_t transform5 = createProjectionMatrix(
      1, 100, (double)SCREEN_HEIGHT / (double)(SCREEN_WIDTH), theta);

  matrix4x4_t transform = m4x4mul(3, transform2, rotation, transform5);

  // forward vector for velocity calculation
  forward = (vector3_t){0, 0, 1};
  matrix4x4_t yrot = createRotationMatrix(0, -phi, 0);
  forward = v3m4x4mul(forward, yrot);
  forward = NORM(forward);

  // get right vector
  right = CROSS(forward, up);
  right = NORM(right);

  vector4_t draworder[mesh.tlen][3];
  int drawcount = 0;
  for (int i = 0; i < mesh.tlen; i++) {
    // for (int i = 0; i < 1; i++) {
    vector4_t p0 = v3tov4(mesh.vertices[mesh.triangles[i][0]]);
    vector4_t p1 = v3tov4(mesh.vertices[mesh.triangles[i][1]]);
    vector4_t p2 = v3tov4(mesh.vertices[mesh.triangles[i][2]]);

    p0 = v4m4x4mul(p0, transform);
    p1 = v4m4x4mul(p1, transform);
    p2 = v4m4x4mul(p2, transform);

    vector3_t line1 = SUB(p1, p0);
    vector3_t line2 = SUB(p2, p0);
    vector3_t norm = CROSS(line1, line2);
    norm = NORM(norm);

    if (DOT(norm, p0) < 0) {
      draworder[drawcount][0] = p0;
      draworder[drawcount][1] = p1;
      draworder[drawcount][2] = p2;
      drawcount++;
    }
  }

  qsort(draworder, drawcount, sizeof(vector4_t) * 3, v4cmpz);

  for (int i = 0; i < drawcount; i++) {
    // for (int i = 0; i < 1; i++) {
    vector4_t p0 = draworder[i][0];
    vector4_t p1 = draworder[i][1];
    vector4_t p2 = draworder[i][2];

    vector3_t line1 = SUB(p1, p0);
    vector3_t line2 = SUB(p2, p0);
    vector3_t norm = CROSS(line1, line2);
    norm = NORM(norm);

    vector3_t light = {0, 0, -1};
    vector3_t lightnorm = NORM(light);
    double luminence = DOT(lightnorm, norm) * 1 / LEN(light);

    double hue, saturation, value;
    hextohsv(mesh.color, &hue, &saturation, &value);

    value *= MAX(luminence, 0.4);

    int color = hsvtohex(hue, saturation, value);

    vector4_t tribuf[3][3];
    int len = cliptriangle(p0, p1, p2, tribuf);

    for (int i = 0; i < len; i++)
      normv2triangle(v4tov2(tribuf[i][0]), v4tov2(tribuf[i][1]),
                     v4tov2(tribuf[i][2]), color);
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
  //  theta = fmod(theta + offsetX * sensitivity / 100, (double)2 * PI);
  /*
  if (offsetX > 2)
    offsetX = 2;
  if (offsetX < -2)
    offsetX = -2;
  phi += ((double)offsetX / 100) * sensitivity;

  if (offsetY > 2)
    offsetY = 2;
  if (offsetY < -2)
    offsetY = -2;
  alpha += ((double)offsetY / 100) * sensitivity;

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
