#include <math.h>
#include <R_ext/Applic.h>
#include <R_ext/Arith.h>
#include <R_ext/Memory.h>

typedef void mode_function(const int *, const int *, const int *,
                           const double *, double *, double *, int *);

extern mode_function fg_long_for_vmmin;
extern mode_function fg_eta_nested_for_vmmin;

typedef struct {
    const int *subject;
    const int *offset;
    const int *dimension;
    mode_function *evaluate;
    double *gradient;
} mode_context;

static double mode_value(int n, double *x, void *data)
{
    mode_context *context = data;
    double value;
    int ok;
    context->evaluate(context->subject, context->offset, context->dimension,
                      x, &value, context->gradient, &ok);
    return ok && R_FINITE(value) ? value : 1e100;
}

static void mode_gradient(int n, double *x, double *gradient, void *data)
{
    mode_context *context = data;
    double value;
    int ok;
    context->evaluate(context->subject, context->offset, context->dimension,
                      x, &value, gradient, &ok);
    for (int k = 0; k < n; ++k) {
        if (!ok || !R_FINITE(gradient[k]))
            gradient[k] = 0.0;
    }
}

static void find_mode(const int *i, const int *sm, const int *dim, double *x,
                      double *value, int *fail, int *fn, int *gr,
                      const int *maxit, const double *abstol, const double *reltol,
                      const int *trace, mode_function *evaluate)
{
    void *stack = vmaxget();
    int *mask = (int *) R_alloc(*dim, sizeof(int));
    double *gradient = (double *) R_alloc(*dim, sizeof(double));
    for (int k = 0; k < *dim; ++k)
        mask[k] = 1;
    mode_context context = {i, sm, dim, evaluate, gradient};
    vmmin(*dim, x, value, mode_value, mode_gradient, *maxit, *trace,
          mask, *abstol, *reltol, 10, &context, fn, gr, fail);
    vmaxset(stack);
}

void vmmin_long_mode(const int *i, const int *sm, const int *dim, double *x,
                     double *value, int *fail, int *fn, int *gr,
                     const int *maxit, const double *abstol, const double *reltol,
                     const int *trace)
{
    find_mode(i, sm, dim, x, value, fail, fn, gr, maxit, abstol, reltol,
              trace, fg_long_for_vmmin);
}

void vmmin_eta_nested_mode(const int *i, const int *sm, const int *dim, double *x,
                           double *value, int *fail, int *fn, int *gr,
                           const int *maxit, const double *abstol, const double *reltol,
                           const int *trace)
{
    find_mode(i, sm, dim, x, value, fail, fn, gr, maxit, abstol, reltol,
              trace, fg_eta_nested_for_vmmin);
}
