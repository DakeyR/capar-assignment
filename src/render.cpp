#include "render.hpp"
#include <cstdint>
#include <cassert>

#define X_Low -2.5
#define X_High 1.0
#define Y_Low -1.0
#define Y_High 1.0

struct rgb8_t {
  std::uint8_t r;
  std::uint8_t g;
  std::uint8_t b;
};

rgb8_t heat_lut(float x)
{
  assert(0 <= x && x <= 1);
  float x0 = 1.f / 4.f;
  float x1 = 2.f / 4.f;
  float x2 = 3.f / 4.f;

  if (x < x0)
  {
    auto g = static_cast<std::uint8_t>(x / x0 * 255);
    return rgb8_t{0, g, 255};
  }
  else if (x < x1)
  {
    auto b = static_cast<std::uint8_t>((x1 - x) / x0 * 255);
    return rgb8_t{0, 255, b};
  }
  else if (x < x2)
  {
    auto r = static_cast<std::uint8_t>((x - x1) / x0 * 255);
    return rgb8_t{r, 255, 0};
  }
  else
  {
    auto b = static_cast<std::uint8_t>((1.f - x) / x0 * 255);
    return rgb8_t{255, b, 0};
  }
}

static float scale(float x, float c_low, float c_high, float mandel_low, float mandel_high)
{
  return mandel_low + ((mandel_high - mandel_low)/(c_high - c_low)) * (x - c_low);
}

void render(std::byte* buffer,
            int width,
            int height,
            std::ptrdiff_t stride,
            int n_iterations)
{
  for (int y = 0; y < height; ++y)
  {
    rgb8_t* lineptr = reinterpret_cast<rgb8_t*>(buffer);

    float y0 = scale((float)y, 0.0, (float)height, Y_Low, Y_High);
    for (int x = 0; x < width; ++x){
      float x0 = scale((float)x, 0.0, (float)width, X_Low, X_High);
      int iteration = 0;
      float nx = 0.0;
      float ny = 0.0;
      while (nx * nx + ny * ny < 4 && iteration < n_iterations) {
        auto xtemp = nx * nx - ny * ny + x0;
        ny = 2 * nx * ny + y0;
        nx = xtemp;
        iteration = iteration + 1;
      }
      lineptr[x] = heat_lut(iteration/n_iterations);
      //lineptr[x] = heat_lut((nx * nx + ny * ny) / float(width * width + height * height));
    }
    buffer += stride;
  }
}


void render_mt(std::byte* buffer,
               int width,
               int height,
               std::ptrdiff_t stride,
               int n_iterations)
{
  render(buffer, width, height, stride, n_iterations);
}
