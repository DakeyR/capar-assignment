#include "render.hpp"
#include <cstdint>
#include <cassert>
#include <vector>
#include <iostream>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>

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
  std::vector<int> histogram(n_iterations + 1, 0);
  std::vector<int> pixels(width * height, 0);

  unsigned int total = 0;
  auto height_limit = height & 0x1 ? (height / 2) + 1 : height / 2;

  for (int j = 0; j < height_limit; ++j)
  {

    double y0 = double(j) / double(height - 1) * 2 - 1;

    for (int i = 0; i < width; ++i){

      double x0 = double(i) / double(width - 1) * 3.5 - 2.5;
      int iteration = 0;
      float x = 0.0;
      float y = 0.0;

      while (x * x + y * y < 2 * 2 && iteration < n_iterations) {
        float xtemp = x * x - y * y + x0;
        y = 2 * x * y + y0;
        x = xtemp;
        iteration = iteration + 1;
      }

      pixels[j * width + i] = iteration;
      //pixels[(height - 1 - j) * width + i] = iteration;
      total++;
      histogram[iteration] += 1;
    }
  }

  total -= histogram[histogram.size() - 1];

  std::byte *buff_start = buffer;
  for (int j = 0; j < height_limit; j++)
  {
    rgb8_t* lineptr = reinterpret_cast<rgb8_t*>(buffer);
    rgb8_t* lineptr_sym = reinterpret_cast<rgb8_t*>(buff_start + stride * (height - 1) - j * stride);
    for (int i = 0; i < width; i++)
    {
      double hue = 0.0;
      int iter = pixels[j * width + i];
      if (iter == n_iterations)
      {
        lineptr[i] = rgb8_t{0,0,0};
        continue;
      }
      for (int k = 0; k <= iter; k++)
        hue += (double(histogram[k]) / double(total));
      lineptr[i] = heat_lut(hue);
      lineptr_sym[i] = lineptr[i];
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
  std::vector<int> histogram(n_iterations + 1, 0);
  std::vector<int> pixels(width * height, 0);

  unsigned int total = 0;
  unsigned long height_limit = height & 0x1 ? (height / 2) + 1 : height / 2;

  auto f = [&](const tbb::blocked_range<size_t>& range){
    for (unsigned j = range.begin(); j < range.end(); ++j)
    {

      double y0 = double(j) / double(height - 1) * 2 - 1;

      for (int i = 0; i < width; ++i)
      {
        double x0 = double(i) / double(width - 1) * 3.5 - 2.5;
        int iteration = 0;
        float x = 0.0;
        float y = 0.0;

        while (x * x + y * y < 2 * 2 && iteration < n_iterations) {
          float xtemp = x * x - y * y + x0;
          y = 2 * x * y + y0;
          x = xtemp;
          iteration = iteration + 1;
        }

        pixels[j * width + i] = iteration;
        //pixels[(height - 1 - j) * width + i] = iteration;
      }
    }
  };

  const tbb::blocked_range<size_t>& half_range{0, height_limit};
  tbb::parallel_for(half_range, f);

  for (int k = 0; k < height_limit * width; k++)
  {
    histogram[pixels[k]]++;
  }

  total = height_limit * width  - histogram[histogram.size() - 1];


  //TODO: Fix range of this loop (computing twice the same lut)...

  std::byte *buff_start = buffer;
  for (int j = 0; j < height_limit; j++)
  {
    rgb8_t* lineptr = reinterpret_cast<rgb8_t*>(buffer);
    rgb8_t* lineptr_sym = reinterpret_cast<rgb8_t*>(buff_start + stride * (height - 1) - j * stride);
    for (int i = 0; i < width; i++)
    {
      double hue = 0.0;
      int iter = pixels[j * width + i];
      if (iter == n_iterations)
      {
        lineptr[i] = rgb8_t{0,0,0};
        continue;
      }
      for (int k = 0; k <= iter; k++)
        hue += (double(histogram[k]) / double(total));
      lineptr[i] = heat_lut(hue);
      lineptr_sym[i] = lineptr[i];
    }
    buffer += stride;
  }
}
