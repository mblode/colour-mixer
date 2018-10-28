#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MAX(a, b) (((a) > (b)) ? (a) : (b))
#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define ROUND(x) ((int)((x) + 0.5))
#define SIGN(x) ((x) > 0 ? 1 : (-1))
#define SQR(x) ((x) * (x))
#define CLAMP(x, low, high) (((x) > (high)) ? (high) : (((x) < (low)) ? (low) : (x)))
#define MAX3(a, b, c) ((a) > (b) ? MAX((a), (c)) : MAX((b), (c)))
#define MIN3(a, b, c) ((a) < (b) ? MIN((a), (c)) : MIN((b), (c)))
#define WGM_EPSILON 0.0001

float T_MATRIX[3][36] = { { 5.47813E-05, 0.000184722, 0.000935514, 0.003096265, 0.009507714, 0.017351596, 0.022073595, 0.016353161, 0.002002407, -0.016177731, -0.033929391, -0.046158952, -0.06381706, -0.083911194, -0.091832385, -0.08258148, -0.052950086, -0.012727224, 0.037413037, 0.091701812, 0.147964686, 0.181542886, 0.210684154, 0.210058081, 0.181312094, 0.132064724, 0.093723787, 0.057159281, 0.033469657, 0.018235464, 0.009298756, 0.004023687, 0.002068643, 0.00109484, 0.000454231, 0.000255925 },
    { -4.65552E-05, -0.000157894, -0.000806935, -0.002707449, -0.008477628, -0.016058258, -0.02200529, -0.020027434, -0.011137726, 0.003784809, 0.022138944, 0.038965605, 0.063361718, 0.095981626, 0.126280277, 0.148575844, 0.149044804, 0.14239936, 0.122084916, 0.09544734, 0.067421931, 0.035691251, 0.01313278, -0.002384996, -0.009409573, -0.009888983, -0.008379513, -0.005606153, -0.003444663, -0.001921041, -0.000995333, -0.000435322, -0.000224537, -0.000118838, -4.93038E-05, -2.77789E-05 },
    { 0.00032594, 0.001107914, 0.005677477, 0.01918448, 0.060978641, 0.121348231, 0.184875618, 0.208804428, 0.197318551, 0.147233899, 0.091819086, 0.046485543, 0.022982618, 0.00665036, -0.005816014, -0.012450334, -0.015524259, -0.016712927, -0.01570093, -0.013647887, -0.011317812, -0.008077223, -0.005863171, -0.003943485, -0.002490472, -0.001440876, -0.000852895, -0.000458929, -0.000248389, -0.000129773, -6.41985E-05, -2.71982E-05, -1.38913E-05, -7.35203E-06, -3.05024E-06, -1.71858E-06 } };

float spectral_r[36] = { 0.022963, 0.022958, 0.022965, 0.022919973, 0.022752, 0.022217, 0.021023, 0.019115, 0.016784, 0.014467, 0.012473, 0.010941, 0.009881, 0.009311, 0.009313, 0.010067, 0.011976, 0.015936, 0.024301, 0.043798, 0.09737, 0.279537, 0.903191, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };

float spectral_g[36] = { 0.023091, 0.023094, 0.02311, 0.023186201, 0.023473, 0.02446, 0.02704, 0.032978, 0.045879, 0.075263, 0.148204, 0.344509, 0.810966, 1, 1, 1, 1, 1, 1, 1, 0.644441, 0.332202, 0.19032, 0.127354, 0.097532, 0.082594, 0.074711, 0.070621, 0.068506, 0.067451, 0.066952, 0.066725, 0.066613, 0.066559, 0.066531, 0.066521 };

float spectral_b[36] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.870246, 0.582997, 0.344759, 0.198566, 0.116401, 0.071107, 0.045856, 0.031371, 0.022678, 0.017256, 0.013751, 0.011428, 0.009878, 0.008831, 0.008138, 0.007697, 0.00743, 0.007271, 0.007182, 0.007136, 0.007111, 0.007099, 0.007094, 0.007092, 0.007091, 0.007089, 0.007089 };

void rgb_to_spectral(float r, float g, float b, float* spectral_)
{

    r = MAX(r, WGM_EPSILON);
    g = MAX(g, WGM_EPSILON);
    b = MAX(b, WGM_EPSILON);
    //upsample rgb to spectral primaries
    float spec_r[36] = { 0 };
    for (int i = 0; i < 36; i++) {
        spec_r[i] = spectral_r[i] * r;
    }
    float spec_g[36] = { 0 };
    for (int i = 0; i < 36; i++) {
        spec_g[i] = spectral_g[i] * g;
    }
    float spec_b[36] = { 0 };
    for (int i = 0; i < 36; i++) {
        spec_b[i] = spectral_b[i] * b;
    }
    //collapse into one spd
    for (int i = 0; i < 36; i++) {
        spectral_[i] += spec_r[i] + spec_g[i] + spec_b[i];
    }
}

void spectral_to_rgb(float* spectral, float* rgb_)
{

    for (int i = 0; i < 36; i++) {
        rgb_[0] += T_MATRIX[0][i] * spectral[i];
        rgb_[1] += T_MATRIX[1][i] * spectral[i];
        rgb_[2] += T_MATRIX[2][i] * spectral[i];
    }
    rgb_[0] = CLAMP(rgb_[0], 0.0f, 1.0f);
    rgb_[1] = CLAMP(rgb_[1], 0.0f, 1.0f);
    rgb_[2] = CLAMP(rgb_[2], 0.0f, 1.0f);
}

//function to make it easy to blend normal and subtractive color blending modes
//a is the first color b is the other color
float* mix_colors(float* a, float* b, float fac, float normsub,
    float smudge_darken, float smudge_desat, float spectral, float multiply)
{
    float rgbmixnorm[4] = { 0 };
    float rgbmixsub[4] = { 0 };
    float rgbmixsub_wgm[4] = { 0 };
    float rgbmixsub_mult[4] = { 0 };
    float spectralmixnorm[4] = { 0 };
    float spectralmixsub[4] = { 0 };
    float spectralmixsub_wgm[4] = { 0 };
    float spectralmixsub_mult[4] = { 0 };
    float normmix[4] = { 0 };
    float submix[4] = { 0 };
    multiply = powf(multiply, 64.0);
    static float result[4] = { 0 };
    #pragma omp threadprivate(result)
    //normsub is the ratio of normal to subtractive
    normsub = CLAMP(normsub, 0.0f, 1.0f);

    //do the mix
    //normal
    if (normsub < 1.0) {

        //RGB normal mixing
        if (spectral < 1.0) {

            rgbmixnorm[0] = fac * a[0] + (1 - fac) * b[0];
            rgbmixnorm[1] = fac * a[1] + (1 - fac) * b[1];
            rgbmixnorm[2] = fac * a[2] + (1 - fac) * b[2];
        }

        //spectral normal mixing
        if (spectral > 0.0) {

            float spec_a[36] = { 0 };
            float spec_b[36] = { 0 };

            rgb_to_spectral(a[0], a[1], a[2], spec_a);
            rgb_to_spectral(b[0], b[1], b[2], spec_b);

            //blend spectral primaries normal mode
            float spectralmix[36] = { 0 };
            for (int i = 0; i < 36; i++) {
                spectralmix[i] = fac * spec_a[i] + (1 - fac) * spec_b[i];
            }
            //convert to RGB
            spectral_to_rgb(spectralmix, spectralmixnorm);
        }

        //mix rgb and spectral normal modes:
        for (int i = 0; i < 3; i++) {
            normmix[i] = (((1 - spectral) * rgbmixnorm[i]) + (spectral * spectralmixnorm[i]));
        }
    }

    //subtractive blend modes (spectral and rgb)
    if (normsub > 0.0) {

        //calculate a different smudge ratio for subtractive modes
        //may not be coherent but it looks good
        //this is a weighted geometric mean method inspired
        //by Scott Allen Burns
        //we want the sum of the ratios to be 1.0
        //but we're using alpha for that ratio so need to scale it
        float paint_a_ratio, paint_b_ratio;
        float subfac = fac;
        float alpha_sum = a[3] + b[3];
        if (alpha_sum > 0.0) {
            paint_a_ratio = (a[3] / alpha_sum) * fac;
            paint_b_ratio = (b[3] / alpha_sum) * (1 - fac);
            paint_a_ratio /= (paint_a_ratio + paint_b_ratio);
            paint_b_ratio /= (paint_a_ratio + paint_b_ratio);
        }

        float paint_ratio_sum = paint_a_ratio + paint_b_ratio;
        if (paint_ratio_sum > 0.0) {
            subfac = paint_a_ratio / (paint_ratio_sum);
        }

        if (spectral < 1.0) {
            //mix with Weighted Geometric Mean
            //don't allow absolute zero because it has infinite
            //darkening power
            if (multiply < 1.0) {
                rgbmixsub_wgm[0] = powf(MAX(a[0], WGM_EPSILON), subfac) * powf(MAX(b[0], WGM_EPSILON), (1 - subfac));
                rgbmixsub_wgm[1] = powf(MAX(a[1], WGM_EPSILON), subfac) * powf(MAX(b[1], WGM_EPSILON), (1 - subfac));
                rgbmixsub_wgm[2] = powf(MAX(a[2], WGM_EPSILON), subfac) * powf(MAX(b[2], WGM_EPSILON), (1 - subfac));
            }
            if (multiply > 0.0) {
                rgbmixsub_mult[0] = MAX(a[0], WGM_EPSILON) * MAX(b[0], WGM_EPSILON);
                rgbmixsub_mult[1] = MAX(a[1], WGM_EPSILON) * MAX(b[1], WGM_EPSILON);
                rgbmixsub_mult[2] = MAX(a[2], WGM_EPSILON) * MAX(b[2], WGM_EPSILON);
            }
            rgbmixsub[0] = rgbmixsub_wgm[0] * (1 - multiply) + rgbmixsub_mult[0] * multiply;
            rgbmixsub[1] = rgbmixsub_wgm[1] * (1 - multiply) + rgbmixsub_mult[1] * multiply;
            rgbmixsub[2] = rgbmixsub_wgm[2] * (1 - multiply) + rgbmixsub_mult[2] * multiply;
        }

        if (spectral > 0.0) {

            float spec_a[36] = { 0 };
            float spec_b[36] = { 0 };

            rgb_to_spectral(a[0], a[1], a[2], spec_a);
            rgb_to_spectral(b[0], b[1], b[2], spec_b);
            //blend spectral primaries subtractive WGM and/or mult
            float spectralmix_wgm[36] = { 0 };
            float spectralmix_mult[36] = { 0 };
            float spectralmix[36] = { 0 };
            for (int i = 0; i < 36; i++) {
                if (multiply < 1.0) {
                    spectralmix_wgm[i] = powf(MAX(spec_a[i], WGM_EPSILON), subfac) * powf(MAX(spec_b[i], WGM_EPSILON), (1 - subfac));
                }
                if (multiply > 0.0) {
                    spectralmix_mult[i] = MAX(spec_a[i], WGM_EPSILON) * MAX(spec_b[i], WGM_EPSILON);
                }
                spectralmix[i] = spectralmix_wgm[i] * (1.0 - multiply) + spectralmix_mult[i] * multiply;
            }
            //convert to RGB
            spectral_to_rgb(spectralmix, spectralmixsub);
        }
        //mix rgb and spectral sub modes:
        for (int i = 0; i < 3; i++) {
            submix[i] = (((1 - spectral) * rgbmixsub[i]) + (spectral * spectralmixsub[i]));
        }
    }

    //combine normal and subtractive RGB modes

    for (int i = 0; i < 3; i++) {
        result[i] = ((1 - normsub) * normmix[i]) + (normsub * submix[i]);
    }

    //alpha is simple
    result[0] = CLAMP(result[0], 0.0f, 1.0f);
    result[1] = CLAMP(result[1], 0.0f, 1.0f);
    result[2] = CLAMP(result[2], 0.0f, 1.0f);
    result[3] = CLAMP(fac * a[3] + (1 - fac) * b[3], 0.0f, 1.0f);

    for (int i = 0; i < 3; i++) {
        if (isnan(result[i])) {
            result[i] = 0.0;
        }
    }

    //Chroma and Luminosity tweak
    //compare result to the smudge_state (a) and tweak the chroma and luma
    //based on hue angle difference in linear HCY
    if ((smudge_desat != 0 || smudge_darken != 0)) {
        float smudge_h = a[0];
        float smudge_c = a[1];
        float smudge_y = a[2];

        float result_h = result[0];
        float result_c = result[1];
        float result_y = result[2];

        rgb_to_hcy_float(&smudge_h, &smudge_c, &smudge_y);
        rgb_to_hcy_float(&result_h, &result_c, &result_y);

        //set our Brightness of the mix according to mode result. and also process saturation
        //Don't process achromatic colors (HCY has 3 achromatic states C=0 or Y=1 or Y=0)
        if (result_c != 0 && smudge_c != 0 && result_y != 0 && smudge_y != 0 && result_y != 1 && smudge_y != 1) {

            //determine hue diff, proportional to smudge ratio
            //if fac is .5 the hueratio should be 1. When fac closer to 0 and closer to 1 should decrease towards zero
            //why- because if smudge is 0 or 1, only 100% of one of the brush or smudge color will be used so there is no comparison to make.
            //when fac is 0.5 the smudge and brush are mixed 50/50, so the huedifference should be respected 100%
            float huediff_sat;
            float huediff_bright;
            float hueratio = (0.5 - fabs(0.5 - fac)) / 0.5;
            float anglediff = fabs(smallest_angular_difference(result_h * 360, smudge_h * 360) / 360);

            //calculate the adjusted hue difference and apply that to the saturation level and/or brightness
            //if smudge_desaturation setting is zero, the huediff will be zero.  Likewise when smudge (fac) is 0 or 1, the huediff will be zero.
            huediff_sat = anglediff * smudge_desat * hueratio;
            //do the same for brightness
            huediff_bright = anglediff * smudge_darken * hueratio;

            //do the desaturation.  More strongly saturated colors will reduce S more than less saturated colors
            result_c = CLAMP(result_c * (1 - huediff_sat), 0.0f, 1.0f);

            //attempt to simulate subtractive mode by darkening colors if they are different hues
            result_y = CLAMP(result_y * (1 - huediff_bright), 0.0f, 1.0f);

            //convert back to sRGB
            hcy_to_rgb_float(&result_h, &result_c, &result_y);

            for (int i = 0; i < 3; i++) {
                if (isnan(result[i])) {
                    result[i] = 0.0;
                }
            }

            result[0] = CLAMP(result_h, 0.0f, 1.0f);
            result[1] = CLAMP(result_c, 0.0f, 1.0f);
            result[2] = CLAMP(result_y, 0.0f, 1.0f);
        }
    }
    return result;
}

int main(int argc, char* argv[])
{

    return 0;
}
