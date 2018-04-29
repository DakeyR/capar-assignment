#include "render.hpp"
#include <cstdint>
#include <cassert>
#include <vector>

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

void render(std::byte* buffer,
            int width,
            int height,
            std::ptrdiff_t stride,
            int n_iterations)
{
  std::vector<int> histogram(n_iterations);
  int pixels[width * height];

  for (int i = 0; i < n_iterations; i++)
    histogram[i] = 0.0;

  //int total = 0;
  for (int j = 0; j < height; ++j)
  {

    float y0 = float(j) / float(height) * 2 - 1;

    for (int i = 0; i < width; ++i){

      float x0 = float(i) / float(width) * 3.5 - 2.5;
      int iteration = 0;
      float x = 0.0;
      float y = 0.0;

      while (x * x + y * y < 2 * 2 && iteration < n_iterations) {
        float xtemp = x * x - y * y + x0;
        y = 2 * x * y + y0;
        x = xtemp;
        iteration = iteration + 1;
      }

//      total += 1;//iteration;
      pixels[j * width + i] = iteration;
      histogram[iteration] += 1;
      //lineptr[x] = heat_lut((nx * nx + ny * ny) / float(width * width + height * height));
    }
  }
  int total = 0;
  for (int val: histogram)
    total += val;

  for (int j = 0; j < height; j++)
  {
    rgb8_t* lineptr = reinterpret_cast<rgb8_t*>(buffer);
    for (int i = 0; i < width; i++)
    {
      float hue = 0.0;
      int iter = pixels[j * width + i];
      if (iter == n_iterations)
      {
        lineptr[i] = rgb8_t{0,0,0};
        continue;
      }
      for (int k = 0; k <= iter; k++)
        hue += (float(histogram[k]) / float(total));
      lineptr[i] = heat_lut(hue);
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
