#include <iostream>
#include <stdio.h>
#include <string>
#include <stdlib.h>
#include <math.h>
#include "ftipsy.hpp"
extern "C" {
#include "jpeglib.h"
}

using namespace std;

typedef struct
{
    char r,g,b;
} rgb_t;

typedef struct
{
    rgb_t **img;
    int **freq;
    float min, max;
    char *fname;
} slice_t;

int nSlices, width, height;

void color_ramp_hot2cold(float *r, float *g, float *b, float v)
{
    float cr, cg, cb;
    cr = cg = cb = 1.0F;

    if (v < 0.0F)
      v = 0;
    if (v > 1.0F)
      v = 1.0F;

    if (v < 0.25) {
      cr = 0;
      cg = 4 * v;
    } else if (v < (0.5)) {
      cr = 0;
      cb = 1 + 4 * (0.25 - v);
    } else if (v < (0.75)) {
      cr = 4 * (v - 0.5);
      cb = 0;
    } else {
      cg = 1 + 4 * (0.75 - v);
      cb = 0;
    }

    *r = cr;
    *g = cg;
    *b = cb;

}

void write_slice(slice_t *slice)
{
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_compress(&cinfo);

    FILE * outfile;

    if ((outfile = fopen(slice->fname, "wb")) == NULL) {
        cerr << "Can't open " << slice->fname << endl;
        exit(1);
    }

    //cerr << "Writing to " << slice->fname << endl;

    jpeg_stdio_dest(&cinfo, outfile);

    cinfo.image_width      = width;
    cinfo.image_height     = height;
    cinfo.input_components = 3;
    cinfo.in_color_space   = JCS_RGB;

    jpeg_set_defaults(&cinfo);
    jpeg_start_compress(&cinfo, TRUE);

    JSAMPROW row_pointer[1];              /* pointer to a single row */
    //int row_stride;                       /* physical row width in buffer */

    //row_stride = cinfo.image_width * 3;   /* JSAMPLEs per row in image_buffer */

    while (cinfo.next_scanline < cinfo.image_height) 
    {
        row_pointer[0] = (JSAMPROW) slice->img[cinfo.next_scanline]; //(env.frame_buffer[(cinfo.image_height-cinfo.next_scanline-1) * row_stride]);
        jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);

    fclose(outfile);
}

int main(int argc, char **argv)
{
    ifTipsy in;

    TipsyHeader h;
    TipsyDarkParticle d;

    if (argc < 5)
    {
        cerr << "Usage: tipmri <tipsy-file> <n-slices> <width> <height>" << endl;
        exit(2);
    }

    nSlices = atoi(argv[2]);

    if (nSlices < 1)
    {
        cerr << "Can only have a positive number of slices." << endl;
        exit(1);
    }

    width  = atoi(argv[3]);
    height = atoi(argv[4]);

    if (width < 0 || height < 0)
    {
        cerr << "Bad width/height" << endl;
        exit(1);
    }

    in.open(argv[1], "standard");
    if (!in.is_open())
    {
        cerr << "Can't open " << argv[1] << endl;
        exit(1);
    }

    in >> h;

    if (h.h_nDark == 0)
    {
        cerr << "No dark particles." << endl;
        exit(1);
    }

    slice_t **slices = new slice_t*[3];

    float sliceDelta = 1.0F / nSlices;


    for (int j=0; j < 3; j++)
    {
        slices[j] = new slice_t[nSlices];

        for (int i=0; i < nSlices; i++)
        {
            slices[j][i].img = (rgb_t **)malloc(height * sizeof(rgb_t *));
            slices[j][i].freq = (int **)malloc(width * sizeof(int *));
            for (int y=0; y < height; y++)
            {
                slices[j][i].img[y]  = (rgb_t *)calloc(width, sizeof(rgb_t));
                slices[j][i].freq[y] = (int *)calloc(width, sizeof(int));
            }

            slices[j][i].min  = i * sliceDelta;
            slices[j][i].min  = i * sliceDelta;
            slices[j][i].max  = slices[j][i].min + sliceDelta;
            slices[j][i].fname = new char[1024];
            snprintf(slices[j][i].fname, 1024, "%s_%.2f-%.2f_%03i.%c.jpg", 
                argv[1], 
                slices[j][i].min, slices[j][i].max,
                i,
                'x' + j);

            cerr << slices[j][i].fname << endl;
        }
    }

    int maxfreq = 0;
    int minfreq = h.h_nDark;

    cerr << "Reading input / creating images..." << endl;
    for (unsigned int i=0; i < h.h_nDark; i++)
    {
        in >> d;

        { /* x axis */
            int s = (int)((d.pos[0] + 0.5) * nSlices); if (s >= nSlices) s = nSlices - 1;

            int x = (int)((d.pos[2] + 0.5) * width);  if (x >= width)  x = width - 1;
            int y = (int)((d.pos[1] + 0.5) * height); if (y >= height) y = height - 1;

            //cerr << "s " << s << " x " << x << " y " << y << endl;

            slices[0][s].freq[y][x]++;
        }

        { /* y axis */
            int s = (int)((d.pos[1] + 0.5) * nSlices); if (s >= nSlices) s = nSlices - 1;

            int x = (int)((d.pos[0] + 0.5) * width);  if (x >= width)  x = width - 1;
            int y = (int)((d.pos[2] + 0.5) * height); if (y >= height) y = height - 1;

            //cerr << "s " << s << " x " << x << " y " << y << endl;

            slices[1][s].freq[y][x]++;
        }

        { /* z axis */
            int s = (int)((d.pos[2] + 0.5) * nSlices); if (s >= nSlices) s = nSlices - 1;

            int x = (int)((d.pos[0] + 0.5) * width);  if (x >= width)  x = width - 1;
            int y = (int)((d.pos[1] + 0.5) * height); if (y >= height) y = height - 1;

            //cerr << "s " << s << " x " << x << " y " << y << endl;

            slices[2][s].freq[y][x]++;
        }

        //if ((i % 1000000L) == 0)
        if ((i & ((1 << 23)-1)) == 0)
            cerr << "\r" << (100.0*i / h.h_nDark) << "%  ";
    }

    cerr << endl;

    in.close();

    cerr << "Finding min/max frequencies..." << endl;
    for (int j=0; j < 3; j++)
        for (int i=0; i < nSlices; i++)
            for (int y=0; y < height; y++)
                for (int x=0; x < width; x++)
                {
                    if (slices[j][i].freq[y][x] > maxfreq) maxfreq = slices[j][i].freq[y][x];
                    if (slices[j][i].freq[y][x] < minfreq) minfreq = slices[j][i].freq[y][x];
                }

    if (maxfreq == minfreq) { minfreq = 0; }

    for (int j=0; j < 3; j++)
        for (int i=0; i < nSlices; i++)
        {
            cerr << "Creating " << slices[j][i].fname << "...";
            cerr.flush();

            for (int y=0; y < height; y++)
                for (int x=0; x < width; x++)
                {
                    float c = 1.0 + 3.0 * slices[j][i].freq[y][x] / (maxfreq - minfreq);
                    c = 255.0 * log(c);

                    //if (j==0) printf("%f\n", c);
                    
                    float r,g,b;
                    //int c = (int)(127.0 * slices[j][i].freq[y][x] / maxfreq + 128.0);
                    color_ramp_hot2cold(&r, &g, &b, c);
                    slices[j][i].img[y][x].r = (int)(255.0 * r);
                    slices[j][i].img[y][x].g = (int)(255.0 * g);
                    slices[j][i].img[y][x].b = (int)(255.0 * b);
                }
            write_slice(&(slices[j][i]));
            cerr << "Done." << endl;
        }

    return 0;
}

