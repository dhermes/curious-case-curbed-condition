# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import fractions

import matplotlib.pyplot as plt
import numpy as np

import de_casteljau
import utils


F = fractions.Fraction
GAMMA15 = utils.gamma(15)
ALPHA = 0.5


def bound4(N):
    r"""Compute an a priori error bound.

    Will be :math:`\gamma_{15} \widetilde{p}(s) / p(s)` where
    :math:`p(s) = (1 - 4s)^5 = [(1 - s) - 3s]^n` and
    :math:`\widetilde{p}(s) = [(1 - s) + 3s]^n = (1 + 2s)^5`

    Evaluations at the point :math:`s = a_N = 1/4 + 6/(16N)`.
    """
    aN = F(1, 4) + F(6, 16 * N)
    p_exact = (1 - 4 * aN) ** 5
    ptilde_exact = (1 + 2 * aN) ** 5
    bound = GAMMA15 * ptilde_exact / p_exact
    return float(abs(bound))


def error4(N):
    aN = F(1, 4) + F(6, 16 * N)
    coeffs = (1, -3, 9, -27, 81, -243)
    p_exact = de_casteljau.basic(aN, coeffs)
    p_computed = de_casteljau.basic(float(aN), coeffs)
    err_exact = abs((F(p_computed) - p_exact) / p_exact)
    return float(err_exact)


def bound5(N):
    r"""Compute an a priori error bound.

    Will be :math:`\gamma_{15} \widetilde{p}(s) / p(s)` where
    :math:`p(s) = (1 - 5s)^5 = [(1 - s) - 4s]^n` and
    :math:`\widetilde{p}(s) = [(1 - s) + 4s]^n = (1 + 3s)^5`

    Evaluations at the point :math:`s = b_N = 1/5 + 8/(25N)`.
    """
    bN = F(1, 5) + F(8, 25 * N)
    p_exact = (1 - 5 * bN) ** 5
    ptilde_exact = (1 + 3 * bN) ** 5
    bound = GAMMA15 * ptilde_exact / p_exact
    return float(abs(bound))


def error5(N):
    bN = F(1, 5) + F(8, 25 * N)
    coeffs = (1, -4, 16, -64, 256, -1024)
    p_exact = de_casteljau.basic(bN, coeffs)
    p_computed = de_casteljau.basic(float(bN), coeffs)
    err_exact = abs((F(p_computed) - p_exact) / p_exact)
    return float(err_exact)


def bound6(N):
    r"""Compute an a priori error bound.

    Will be :math:`\gamma_{15} \widetilde{p}(s) / p(s)` where
    :math:`p(s) = (1 - 6s)^5 = [(1 - s) - 5s]^n` and
    :math:`\widetilde{p}(s) = [(1 - s) + 5s]^n = (1 + 4s)^5`

    Evaluations at the point :math:`s = c_N = 1/6 + 10/(36N)`.
    """
    cN = F(1, 6) + F(10, 36 * N)
    p_exact = (1 - 6 * cN) ** 5
    ptilde_exact = (1 + 4 * cN) ** 5
    bound = GAMMA15 * ptilde_exact / p_exact
    return float(abs(bound))


def error6(N):
    cN = F(1, 6) + F(10, 36 * N)
    coeffs = (1, -5, 25, -125, 625, -3125)
    p_exact = de_casteljau.basic(cN, coeffs)
    p_computed = de_casteljau.basic(float(cN), coeffs)
    err_exact = abs((F(p_computed) - p_exact) / p_exact)
    return float(err_exact)


def main(filename=None):
    figure = plt.figure()
    ax = figure.gca()

    bounds = []
    for exponent in range(1, 45 + 1):
        N = 2.1 ** exponent
        fN = F(N)
        bounds.append(
            (
                N,
                bound4(fN),
                bound5(fN),
                bound6(fN),
                error4(fN),
                error5(fN),
                error6(fN),
            )
        )

    bounds = np.array(bounds)
    ax.loglog(
        bounds[:, 0], bounds[:, 1], alpha=ALPHA, color="black", label="Bound"
    )
    # NOTE: We intentionally omit bound5() and bound6() since they are
    #       essentially identical.
    ax.loglog(
        bounds[:, 0],
        bounds[:, 4],
        marker="o",
        linestyle="none",
        markersize=7,
        markeredgewidth=1,
        markerfacecolor="none",
        label="$u(s)$",
    )
    ax.loglog(
        bounds[:, 0],
        bounds[:, 5],
        marker="d",
        linestyle="none",
        label="$v(s)$",
    )
    ax.loglog(
        bounds[:, 0],
        bounds[:, 6],
        marker="o",
        linestyle="none",
        markersize=4,
        color="black",
        zorder=2,
        label="$w(s)$",
    )
    # Label the axes.
    ax.set_xlabel("$N$")
    ax.set_ylabel("Relative Forward Error")
    # Add the legend.
    ax.legend(framealpha=1.0, frameon=True)
    # Set the major x- and y-ticks.
    ax.set_xticks([1e0, 1e3, 1e6, 1e9, 1e12])
    ax.set_yticks([1e-15, 1e0, 1e15, 1e30, 1e45, 1e60])

    if filename is None:
        plt.show()
    else:
        path = utils.get_path(filename)
        figure.savefig(path, bbox_inches="tight")
        print("Saved {}".format(filename))
        plt.close(figure)


if __name__ == "__main__":
    utils.set_styles()
    main(filename="against_a_priori.pdf")
