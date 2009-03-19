#ifndef TREE_JS
#define TREE_JS

struct bound 
{
    real r[3];
    real max[3];
};

struct tree 
{
    real r[3];
    real mass;

    struct bound bnd;

    struct tree *left;
    struct tree *right;

    int iUpper, iLower;
};

void tree_build_js(struct env *env);
void tree_free_js(struct env *env);

#endif

