#include "render.hpp"
#include <cstdint>
#include <cassert>
#include <vector>
#include <iostream>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <immintrin.h>
#include <numeric>
#include <fstream>

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
  //assert(0 <= x && x <= 1);
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

  __m256 limit = _mm256_set1_ps(4);
  __m256 incr = _mm256_set1_ps(1);
  __m256 mask = _mm256_set1_ps(-1);
  unsigned int total = 0;
  auto height_limit = height & 0x1 ? (height / 2) + 1 : height / 2;

  for (int j = 0; j < height_limit; ++j)
  {

    double y0 = j / double(height - 1) * 2 - 1;
    __m256 y0s = _mm256_set1_ps(y0);

    for (int i = 0; i < width; i += 8){

      double x0 = i / double(width - 1) * 3.5 - 2.5;
      __m256 x0s = _mm256_set_ps((i + 7)/ double(width - 1) * 3.5 - 2.5,
                                (i + 6) / double(width - 1) * 3.5 - 2.5,
                                (i + 5) / double(width - 1) * 3.5 - 2.5,
                                (i + 4) / double(width - 1) * 3.5 - 2.5,
                                (i + 3) / double(width - 1) * 3.5 - 2.5,
                                (i + 2) / double(width - 1) * 3.5 - 2.5,
                                (i + 1) / double(width - 1) * 3.5 - 2.5,
                                i / double(width - 1) * 3.5 - 2.5);
      int iteration = 0;
      //float x = 0.0;
      //float y = 0.0;

      __m256 xs = _mm256_setzero_ps();
      __m256 ys = _mm256_setzero_ps();
      __m256 sx = _mm256_setzero_ps();
      __m256 sy = _mm256_setzero_ps();
      __m256 xy = _mm256_setzero_ps();
      __m256 xtmps = _mm256_setzero_ps();
      //while (x * x + y * y < 4.0 && iteration < n_iterations) {
      __m256 cond = _mm256_setzero_ps();

      __m256 iters = _mm256_setzero_ps();

      do {
        //float xtemp = x * x - y * y + x0;
        sx = xs * xs;
        sy = ys * ys;
        xy = xs * ys;

        xs = sx - sy + x0s;
        ys = xy + xy + y0s;
        //x = xtemp;
        cond = sx + sy < limit;
        iteration = iteration + 1;
        iters += _mm256_and_ps(cond, incr);
      } while (!_mm256_testz_ps(cond, mask) && iteration < n_iterations);

      //pixels[j * width + i] = iteration;
      //__m256i tmp_iters = _mm256_cvtps_epi32(iters);
      pixels[j * width + i] = (int)iters[0];
      pixels[j * width + i + 1] = (int)iters[1];
      pixels[j * width + i + 2] = (int)iters[2];
      pixels[j * width + i + 3] = (int)iters[3];
      pixels[j * width + i + 4] = (int)iters[4];
      pixels[j * width + i + 5] = (int)iters[5];
      pixels[j * width + i + 6] = (int)iters[6];
      pixels[j * width + i + 7] = (int)iters[7];
      //histogram[iteration] += 1;
    }
  }

  //std::ofstream of("simd.txt", std::ios::trunc);
  for (int k = 0; k < height_limit * width; k++)
  {
    histogram[pixels[k]]++;
    //of << pixels[k] << std::endl;
  }

  total = height_limit * width  - histogram[histogram.size() - 1];

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
        hue += (histogram[k] / double(total));
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

  __m256 limit = _mm256_set1_ps(4);
  __m256 incr = _mm256_set1_ps(1);
  __m256 mask = _mm256_set1_ps(-1);

  auto f = [&](const tbb::blocked_range<size_t>& range){
    for (int j = range.begin(); j < range.end(); ++j)
    {

      double y0 = j / double(height - 1) * 2 - 1;
      __m256 y0s = _mm256_set1_ps(y0);

      for (int i = 0; i < width; i += 8){

        double x0 = i / double(width - 1) * 3.5 - 2.5;
        __m256 x0s = _mm256_set_ps((i + 7)/ double(width - 1) * 3.5 - 2.5,
                                  (i + 6) / double(width - 1) * 3.5 - 2.5,
                                  (i + 5) / double(width - 1) * 3.5 - 2.5,
                                  (i + 4) / double(width - 1) * 3.5 - 2.5,
                                  (i + 3) / double(width - 1) * 3.5 - 2.5,
                                  (i + 2) / double(width - 1) * 3.5 - 2.5,
                                  (i + 1) / double(width - 1) * 3.5 - 2.5,
                                  i / double(width - 1) * 3.5 - 2.5);
        int iteration = 0;

        __m256 xs = _mm256_setzero_ps();
        __m256 ys = _mm256_setzero_ps();
        __m256 sx = _mm256_setzero_ps();
        __m256 sy = _mm256_setzero_ps();
        __m256 xy = _mm256_setzero_ps();
        __m256 xtmps = _mm256_setzero_ps();
        __m256 cond = _mm256_setzero_ps();
        __m256 iters = _mm256_setzero_ps();

        do {
          sx = xs * xs;
          sy = ys * ys;
          xy = xs * ys;

          xs = sx - sy + x0s;
          ys = xy + xy + y0s;
          cond = sx + sy < limit;
          iteration = iteration + 1;
          iters += _mm256_and_ps(cond, incr);
        } while (!_mm256_testz_ps(cond, mask) && iteration < n_iterations);

        pixels[j * width + i] = (int)iters[0];
        pixels[j * width + i + 1] = (int)iters[1];
        pixels[j * width + i + 2] = (int)iters[2];
        pixels[j * width + i + 3] = (int)iters[3];
        pixels[j * width + i + 4] = (int)iters[4];
        pixels[j * width + i + 5] = (int)iters[5];
        pixels[j * width + i + 6] = (int)iters[6];
        pixels[j * width + i + 7] = (int)iters[7];
      }
    }
  };

  const tbb::blocked_range<size_t> half_range{0, height_limit};
  tbb::parallel_for(half_range, f);

  for (int k = 0; k < height_limit * width; k++)
  {
    histogram[pixels[k]]++;
  }

  total = height_limit * width  - histogram[histogram.size() - 1];


  //TODO: Fix range of this loop (computing twice the same lut)...

  auto f2 = [&](const tbb::blocked_range<size_t>& range){
    for (int j = range.begin(); j < range.end(); j++)
    {
      rgb8_t* lineptr = reinterpret_cast<rgb8_t*>(buffer + j * stride);
      rgb8_t* lineptr_sym = reinterpret_cast<rgb8_t*>(buffer + stride * (height - 1) - j * stride);
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
    }
  };
  tbb::parallel_for(half_range, f2);
}
