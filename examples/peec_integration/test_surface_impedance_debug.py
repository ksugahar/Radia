"""
Debug surface impedance calculation.
"""

import numpy as np


def dowell_formula(xi):
    """Standard Dowell's formula."""
    if xi < 0.01:
        return 1.0 + xi**4 / 45

    sh2 = np.sinh(2*xi)
    sn2 = np.sin(2*xi)
    ch2 = np.cosh(2*xi)
    cs2 = np.cos(2*xi)

    return xi * (sh2 + sn2) / (ch2 - cs2)


def surface_impedance_debug(xi):
    """Debug version of surface impedance method."""
    gamma_a = complex(1, 1) * xi

    tanh_ga = np.tanh(gamma_a)

    # Z_s / Z_dc = (gamma * a) * tanh(gamma * a) / (gamma * a) * a
    #            = tanh(gamma * a) * a * gamma / (rho/a) * rho
    # Hmm, let me redo this

    # Z_s = rho * gamma * tanh(gamma*a)
    # Z_dc = rho / a
    # R_ac/R_dc = Re(Z_s) / Z_dc = Re(rho * gamma * tanh(gamma*a)) / (rho/a)
    #           = a * Re(gamma * tanh(gamma*a))

    # gamma = (1+j) / delta
    # gamma * a = (1+j) * a / delta = (1+j) * xi

    # So: R_ac/R_dc = Re((gamma*a) * tanh(gamma*a)) = Re((1+j)*xi * tanh((1+j)*xi))

    product = gamma_a * tanh_ga
    result = product.real

    print(f"xi = {xi:.2f}")
    print(f"  gamma*a = {gamma_a}")
    print(f"  tanh(gamma*a) = {tanh_ga}")
    print(f"  product = {product}")
    print(f"  Re(product) = {result}")
    print(f"  Dowell = {dowell_formula(xi):.4f}")
    print()

    return result


def main():
    print("=" * 60)
    print("Surface Impedance Debug")
    print("=" * 60)

    xi_values = [0.1, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0]

    for xi in xi_values:
        surface_impedance_debug(xi)


if __name__ == '__main__':
    main()
