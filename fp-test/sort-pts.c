#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NPOINTS 10

typedef struct
{
    float r[3];
    unsigned long key;
} Point;

#define I(x,y) ((*((int *)&x)) & y)
#define make3key0(x) (I(x[0], 0x49249249) | I(x[1], 0x92492492) | I(x[2], 0x24924924))
#define make3key1(x) (I(x[1], 0x49249249) | I(x[2], 0x92492492) | I(x[0], 0x24924924))
#define make3key2(x) (I(x[2], 0x49249249) | I(x[0], 0x92492492) | I(x[1], 0x24924924))

#define make2key0(x) (I(x[1], 0x55555555) | I(x[0], 0xAAAAAAAA))
#define make2key1(x) (I(x[0], 0x55555555) | I(x[1], 0xAAAAAAAA))

unsigned long makekey(float r[3])
{
    unsigned long z = 0;
    unsigned int x, y;
    float f;

    f = r[0]; x = *((unsigned int *)&f);
    f = r[1]; y = *((unsigned int *)&f);

    for (int i = 0; i < sizeof(x) * 8; i++) // unroll for more speed...
          z |= (x & 1UL << i) << i | (y & 1UL << i) << (i + 1);

    //fprintf(stderr, "%x %x %lx\n", x, y, z);
    return z;
}

int compar(const void *a, const void *b)
{
    Point *a0 = (Point *)a;
    Point *b0 = (Point *)b;

    if (a0->key < b0->key) 
        return -1;
    else if (a0->key > b0->key) 
        return 1;
    else
        return 0;
}

char *tobinstr(long x, char *buf)
{
    unsigned long mask;
    char *ptr = buf;
    for (mask=0x8000000000000000UL; mask; mask >>= 1)
    {
        if (x & mask)
            *ptr++ = '1';
        else
            *ptr++ = '0';
    }

    *ptr = 0;

    return buf;
}

int main(int argc, char **argv)
{
    int i;
    float ran;
    char buf[65];
    Point p[NPOINTS];
    Point corner_min = { {1024, 1024, 1024}, 0 };
    Point corner_max = { {0, 0, 0}, 0 };

    srand(0);
    for (i=0; i < NPOINTS; i++)
    {
        ran = (rand() % 1024 + 1); p[i].r[0] = ran;
        ran = (rand() % 1024 + 1); p[i].r[1] = ran;
        ran = (rand() % 1024 + 1); p[i].r[2] = ran;


        corner_min.r[0] = fminf(corner_min.r[0], p[i].r[0]);
        corner_min.r[1] = fminf(corner_min.r[1], p[i].r[1]);
        corner_min.r[2] = fminf(corner_min.r[2], p[i].r[2]);

        corner_max.r[0] = fmaxf(corner_max.r[0], p[i].r[0]);
        corner_max.r[1] = fmaxf(corner_max.r[1], p[i].r[1]);
        corner_max.r[2] = fmaxf(corner_max.r[2], p[i].r[2]);
    }
    //printf("%f %f %f\n", corner.r[0], corner.r[1], corner.r[2]);

    int start_with_x 
        =  fabsf(corner_max.r[0] - corner_min.r[0]) 
        >= fabsf(corner_max.r[1] - corner_min.r[1]);

    fprintf(stderr, "start_with_x: %i\n", start_with_x);

    for (i=0; i < NPOINTS; i++)
    {
#if 0
        p[i].r[0] -= (corner_max.r[0] + corner_min.r[0]) / 2;
        p[i].r[1] -= (corner_max.r[1] + corner_min.r[1]) / 2;
        p[i].r[2] -= (corner_max.r[2] + corner_min.r[2]) / 2;
#else
        p[i].r[0] -= corner_min.r[0];
        p[i].r[1] -= corner_min.r[1];
        p[i].r[2] -= corner_min.r[2];
#endif

        //p[i].r[0] *= -1;
        //p[i].r[1] *= -1;
        //p[i].r[2] *= -1;

#if 0
        if (start_with_x)
            p[i].key = make2key0(p[i].r);
        else
            p[i].key = make2key1(p[i].r);
#endif
        p[i].key = makekey(p[i].r);

        printf("%f %f %f\n", p[i].r[0], p[i].r[1], p[i].r[2]);
        //fprintf(stderr, "%s\n", tobinstr(p[i].key, buf));
    }

    qsort(p, NPOINTS, sizeof(Point), compar);

    printf("\n\n");

        printf("%f\t(%8x)\t%f\t(%8x)\t%f\t(%8x)\t0x%lx\t%s\n", 
            p[0].r[0], I(p[0].r[0], 0xFFFFFFFF),
            p[0].r[1], I(p[0].r[1], 0xFFFFFFFF),
            p[0].r[2], I(p[0].r[2], 0xFFFFFFFF),
            p[0].key, 
            tobinstr(p[0].key, buf));
    for (i=1; i < NPOINTS; i++)
    {
        printf("%f\t(%8x)\t%f\t(%8x)\t%f\t(%8x)\t0x%lx\t%s\n", 
            p[i].r[0], I(p[i].r[0], 0xFFFFFFFF),
            p[i].r[1], I(p[i].r[1], 0xFFFFFFFF),
            p[i].r[2], I(p[i].r[2], 0xFFFFFFFF),
            p[i].key, 
            tobinstr(p[i].key, buf));
    }

    return 0;
}
