import numpy as np
import math

# weighted geometric mean must avoid absolute zero
_WGM_EPSILON = 0.000030517
# _WGM_EPSILON = 0.000030517

def RGB_to_Spectral(rgb):
    """Converts RGB to 36 segments spectral power distribution curve.
    Upsamples to spectral primaries and sums them together into one SPD
    Based on work by Scott Allen Burns.
    """

    r, g, b = rgb

    r = r / 255
    g = g / 255
    b = b / 255

    r = max(r, _WGM_EPSILON)
    g = max(g, _WGM_EPSILON)
    b = max(b, _WGM_EPSILON)

    # Spectral primaries derived by an optimization routine devised by
    # Allen Burns. Smooth curves <= 1.0 to match XYZ
    spectral_r = r * np.array(
        [0.022963, 0.022958, 0.022965, 0.022919973, 0.022752, 0.022217,
         0.021023, 0.019115, 0.016784, 0.014467, 0.012473, 0.010941,
         0.009881, 0.009311, 0.009313, 0.010067, 0.011976, 0.015936,
         0.024301, 0.043798, 0.09737, 0.279537, 0.903191, 1, 1, 1, 1,
         1, 1, 1, 1, 1, 1, 1, 1, 1
         ]
    )

    spectral_g = g * np.array(
        [0.023091, 0.023094, 0.02311, 0.023186201, 0.023473, 0.02446,
         0.02704, 0.032978, 0.045879, 0.075263, 0.148204, 0.344509,
         0.810966, 1, 1, 1, 1, 1,
         1, 1, 0.644441, 0.332202, 0.19032, 0.127354,
         0.097532, 0.082594, 0.074711, 0.070621, 0.068506, 0.067451,
         0.066952, 0.066725, 0.066613, 0.066559, 0.066531, 0.066521
         ]
    )

    spectral_b = b * np.array(
        [1, 1, 1, 1, 1, 1,
         1, 1, 1, 0.870246, 0.582997, 0.344759,
         0.198566, 0.116401, 0.071107, 0.045856, 0.031371, 0.022678,
         0.017256, 0.013751, 0.011428, 0.009878, 0.008831, 0.008138,
         0.007697, 0.00743, 0.007271, 0.007182, 0.007136, 0.007111,
         0.007099, 0.007094, 0.007092, 0.007091, 0.007089, 0.007089
         ]
    )

    x = np.sum([spectral_r, spectral_g, spectral_b], axis=0)
    return x

def Spectral_to_RGB(spd):
    """Converts 36 segments spectral power distribution curve to RGB.
    Based on work by Scott Allen Burns.
    """

    # Spectral_to_XYZ CIE matrix weighted w/ diagonal D65 matrix.  sRGB
    T_MATRIX = (
        np.array(
            [[5.47813E-05, 0.000184722, 0.000935514, 0.003096265,
              0.009507714, 0.017351596, 0.022073595, 0.016353161,
              0.002002407, -0.016177731, -0.033929391, -0.046158952,
              -0.06381706, -0.083911194, -0.091832385, -0.08258148,
              -0.052950086, -0.012727224, 0.037413037, 0.091701812,
              0.147964686, 0.181542886, 0.210684154, 0.210058081,
              0.181312094, 0.132064724, 0.093723787, 0.057159281,
              0.033469657, 0.018235464, 0.009298756, 0.004023687,
              0.002068643, 0.00109484, 0.000454231, 0.000255925],
             [-4.65552E-05, -0.000157894, -0.000806935, -0.002707449,
              -0.008477628, -0.016058258, -0.02200529, -0.020027434,
              -0.011137726, 0.003784809, 0.022138944, 0.038965605,
              0.063361718, 0.095981626, 0.126280277, 0.148575844,
              0.149044804, 0.14239936, 0.122084916, 0.09544734,
              0.067421931, 0.035691251, 0.01313278, -0.002384996,
              -0.009409573, -0.009888983, -0.008379513, -0.005606153,
              -0.003444663, -0.001921041, -0.000995333, -0.000435322,
              -0.000224537, -0.000118838, -4.93038E-05, -2.77789E-05],
             [0.00032594, 0.001107914, 0.005677477, 0.01918448,
              0.060978641, 0.121348231, 0.184875618, 0.208804428,
              0.197318551, 0.147233899, 0.091819086, 0.046485543,
              0.022982618, 0.00665036, -0.005816014, -0.012450334,
              -0.015524259, -0.016712927, -0.01570093, -0.013647887,
              -0.011317812, -0.008077223, -0.005863171, -0.003943485,
              -0.002490472, -0.001440876, -0.000852895, -0.000458929,
              -0.000248389, -0.000129773, -6.41985E-05, -2.71982E-05,
              -1.38913E-05, -7.35203E-06, -3.05024E-06, -1.71858E-06]]
        )
    )
    r, g, b = np.sum(spd*T_MATRIX, axis=1)


    r = r * 255
    g = g * 255
    b = b * 255

    print(r)

    return int(r), int(g), int(b)

def Spectral_Mix_WGM(spd_a, spd_b, ratio):
    """Mixes two SPDs via weighted geomtric mean and returns an SPD.
    Based on work by Scott Allen Burns.
    """
    print(spd_a**(1.0 - ratio) * spd_b**ratio)
    return spd_a**(1.0 - ratio) * spd_b**ratio

# rgb_color_1 = input('Enter hex 1: ').lstrip('#')
# rgb_color_2 = input('Enter hex 2: ').lstrip('#')

# rgb_color_1 = tuple(int(rgb_color_1[i:i+2], 16) for i in (0, 2, 4))
# rgb_color_2 = tuple(int(rgb_color_2[i:i+2], 16) for i in (0, 2, 4))

rgb_color_1 = (255, 0, 0)
rgb_color_2 = (255, 255, 0)

spectral_color_1 = RGB_to_Spectral(rgb_color_1)
spectral_color_2 = RGB_to_Spectral(rgb_color_2)

# blend = float(input('Enter blend (0 to 1): '))
blend = 0

spd_blended = Spectral_Mix_WGM(spectral_color_1, spectral_color_2, blend)

rgb_blended = Spectral_to_RGB(spd_blended)
# rgb_blended = '#%02x%02x%02x' % rgb_blended

print(rgb_blended)
