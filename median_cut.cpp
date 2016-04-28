/*
 * median_cut
 *
 * Copyright (c) 2013 Tobias Alexander Franke
 * http://www.tobias-franke.eu
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 * To get this to run, you will need stb_image.c and stb_image_write.h:
 * http://nothings.org/stb_image.c
 * http://nothings.org/stb/stb_image_write.h
 */

#include <iostream>
#include <vector>
#include <cassert>

#define STB_IMAGE_WRITE_IMPLEMENTATION

#include "stb/stb_image.c"
#include "stb/stb_image_write.h"

struct float2
{
    float x, y;
};

struct float3
{
    float x, y, z;
};

template<typename T>
float luminance(T r, T g, T b)
{
    return r*0.2125f + g*0.7154f + b*0.0721f;
}

/**
 * Summed Area Table
 *
 * Create a luminance summed area table from an image.
 */
class summed_area_table
{
protected:
    int width_, height_;
    std::vector<float> sat_;

    float I(int x, int y) const
    {
        if (x < 0 || y < 0) return 0;
        size_t i = y*width_ + x;
        return sat_[i];
    }

public:
    template<typename T>
    void create_lum(T* rgb, int width, int height, int nc)
    {
        assert(nc > 2);

        width_ = width; height_ = height;

        sat_.clear();
        sat_.resize(width_ * height_);

        for (int y = 0; y < height_; ++y)
        for (int x = 0; x < width_;  ++x)
        {
            size_t i = y*width_ + x;

            T r = rgb[i*nc + 0];
            T g = rgb[i*nc + 1];
            T b = rgb[i*nc + 2];

            float ixy = luminance(r,g,b);

            sat_[i] = ixy + I(x-1, y) + I(x, y-1) - I(x-1, y-1);
        }
    }

    int width() const  { return width_;  }
    int height() const { return height_; }

    /**
     * Returns the sum of a region defined by A,B,C,D.
     *
     * A----B
     * |    |  sum = C+A-B-D
     * D----C
     */
    int sum(int ax, int ay, int bx, int by, int cx, int cy, int dx, int dy) const
    {
        return I(cx, cy) + I(ax, ay) - I(bx, by) - I(dx, dy);
    }
};

/**
 * A subregion in a summed_area_table.
 */
struct sat_region
{
    int x_, y_, w_, h_;
    float sum_;
    const summed_area_table* sat_;

    void create(int x, int y, int w, int h, const summed_area_table* sat, float init_sum = -1)
    {
        x_ = x; y_ = y; w_ = w; h_ = h; sum_ = init_sum; sat_ = sat;

        if (sum_ < 0)
            sum_ = sat_->sum(x,       y,
                             x+(w-1), y,
                             x+(w-1), y+(h-1),
                             x,       y+(h-1));
    }

    void split_w(sat_region& A) const
    {
        for (size_t w = 1; w <= w_; ++w)
        {
            A.create(x_, y_, w, h_, sat_);

            // if region left has approximately half the energy of the entire thing stahp
            if (A.sum_*2.f >= sum_)
                break;
        }
    }

    /**
     * Split region horizontally into subregions A and B.
     */
    void split_w(sat_region& A, sat_region& B) const
    {
        split_w(A);
        B.create(x_ + (A.w_-1), y_, w_ - A.w_, h_, sat_, sum_ - A.sum_);
    }

    void split_h(sat_region& A) const
    {
        for (size_t h = 1; h <= h_; ++h)
        {
            A.create(x_, y_, w_, h, sat_);

            // if region top has approximately half the energy of the entire thing stahp
            if (A.sum_*2.f >= sum_)
                break;
        }
    }

    /**
     * Split region vertically into subregions A and B.
     */
    void split_h(sat_region& A, sat_region& B) const
    {
        split_h(A);
        B.create(x_, y_ + (A.h_-1), w_, h_ - A.h_, sat_, sum_ - A.sum_);
    }

    float2 centroid() const
    {
        float2 c;

        sat_region A;

        split_w(A);
        c.x = A.x_ + (A.w_-1);

        split_h(A);
        c.y = A.y_ + (A.h_-1);

        return c;
    }
};

/**
 * Recursively split a region r and append new subregions
 * A and B to regions vector when at an end.
 */
void split_recursive(const sat_region& r, size_t n, std::vector<sat_region>& regions)
{
    // check: can't split any further?
    if (r.w_ < 2 || r.h_ < 2 || n == 0)
    {
        regions.push_back(r);
        return;
    }

    sat_region A, B;

    if (r.w_ > r.h_)
        r.split_w(A, B);
    else
        r.split_h(A, B);

    split_recursive(A, n-1, regions);
    split_recursive(B, n-1, regions);
}

/**
 * The median cut algorithm.
 *
 * img - Summed area table of an image
 * n - number of subdivision, yields 2^n cuts
 * regions - an empty vector that gets filled with generated regions
 */
void median_cut(const summed_area_table& img, size_t n, std::vector<sat_region>& regions)
{
    regions.clear();

    // insert entire image as start region
    sat_region r;
    r.create(0, 0, img.width(), img.height(), &img);

    // recursively split into subregions
    split_recursive(r, n, regions);
}

/**
 * Create a light source from each region by querying its centroid
 */
void red(float* d, int ci, int m)
{
    if (ci < 0) return;
    if (ci > m) return;

    d[ci + 0] = 1.f;
    d[ci + 1] = 0.f;
    d[ci + 2] = 0.f;
    d[ci + 3] = 1.f;
}

/**
 * Draw a cross at position l into image rgba
 */
void draw(float* rgba, int width, int height, float2 l)
{
    static int i = 0;

    int ci;

    int m = width*height*4;

    for (int x = -1; x < 2; ++x)
    {
        ci = std::min<int>((l.y*width + l.x+x)*4, m);
        red(rgba, ci, m);
    }

    ci = std::min<int>(((l.y+1)*width + l.x)*4, m);
    red(rgba, ci, m);

    ci = std::min<int>(((l.y-1)*width + l.x)*4, m);
    red(rgba, ci, m);
}

/**
 * Create a light source position from each region by querying its centroid
 */
void create_lights(const std::vector<sat_region>& regions, std::vector<float2>& lights)
{
    // set light at centroid
    for (size_t i = 0; i < regions.size(); ++i)
        lights.push_back(regions[i].centroid());
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cerr << "Use " << argv[0] << " filename" << std::endl;
        return 1;
    }

    // load image
    int width, height, nc;
    float* rgba = stbi_loadf(argv[1], &width, &height, &nc, 4);
    if (stbi_failure_reason())
    {
        std::cerr << "stbi: " << stbi_failure_reason() << std::endl;
        return 1;
    }

    // create summed area table of luminance image
    summed_area_table lum_sat;
    lum_sat.create_lum(rgba, width, height, 4);

    // apply median cut
    std::vector<sat_region> regions;
    median_cut(lum_sat, 9, regions); // max 2^n cuts

    // create 2d positions from regions
    std::vector<float2> lights;
    create_lights(regions, lights);

    // draw a marker into image for each position
    size_t i = 0;
    for (auto l = lights.begin(); l != lights.end(); ++l)
    {
        std::cout << "Light " << i++ << ": (" << l->x << ", " << l->y << ")" << std::endl;
        draw(rgba, width, height, *l);
    }

    // save image with marked samples
    std::vector<unsigned char> conv;
    conv.resize(width*height*4);

    for (size_t i = 0; i < width * height * 4; ++i)
        conv[i] = static_cast<unsigned char>(rgba[i]*255);

    stbi_write_bmp("test.bmp", width, height, 4, &conv[0]);

    return 0;
}
