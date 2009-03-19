/* Just shows that if a and b are integers and a < b then
 * if we interpret the bits of a and b as floats it still
 * holds that a < b (as floats).
 */

#include <stdio.h>

char *tobinstr(int x, char *buf)
{
    unsigned int mask;
    char *ptr = buf;
    for (mask=0x80000000; mask; mask >>= 1)
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
    char buf[33];

    //for (i=0; i < 0xFFFFFFFF; i+= 1000)
    for (i=-2; i <= 2 ; i+= 4)
    {
        float f  = i; //*((float *)&i);
        int   fi = *((int *)&f);
        printf("%24.50f %s\n", f, tobinstr(fi, buf));
    }


    return 0;
}
