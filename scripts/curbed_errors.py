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
GAMMA3 = utils.gamma(3)
ALPHA = 0.5


def custom_de_casteljau(s, n):
    """Evaluate :math:`(1 - 5s)^n` in a "special" way.

    Does so via ``p{n} = 1``, ``p{k} = (1 - s) p{k+1} - 4s p{k + 1}``.
    """
    p = 1
    r = 1 - s
    scaled = 4 * s
    for _ in range(n):
        p = r * p - scaled * p
    return p


def abs_phi(s):
    s = F(s)
    return abs((1 + 3 * s) / (1 - 5 * s))


def bounds_curbed(s, n):
    phi = abs_phi(s)
    bound1 = utils.gamma(3 * n) * phi ** n
    bound2 = (1 + phi * GAMMA3) ** n - 1
    if bound1 <= 0:
        raise ValueError(bound1, float(bound1))
    if bound2 <= 0:
        raise ValueError(bound2, float(bound2))

    computed_p = custom_de_casteljau(s, n)
    exact_p = custom_de_casteljau(F(s), n)
    if not isinstance(exact_p, F):
        raise TypeError(exact_p)
    observed_err = abs((computed_p - exact_p) / exact_p)

    return float(bound1), float(bound2), float(observed_err)


def main(filename=None):
    bound_vals = []

    for exponent in range(1, 45 + 1):
        N = 2.1 ** exponent
        bN = F(1, 5) + F(8, 25 * F(N))
        bound1, bound2, observed_err = bounds_curbed(float(bN), 5)
        bound_vals.append((N, bound1, bound2, observed_err))

    bound_vals = np.array(bound_vals)

    figure = plt.figure()
    ax = figure.gca()
    # Add the "curbed" plot.
    ax.loglog(
        bound_vals[:, 0],
        bound_vals[:, 1],
        color="black",
        alpha=ALPHA,
        linestyle=":",
        label=r"Na\"ive Bound",
    )
    ax.loglog(
        bound_vals[:, 0],
        bound_vals[:, 2],
        color="black",
        alpha=ALPHA,
        label="Improved Bound",
    )
    ax.loglog(
        bound_vals[:, 0],
        bound_vals[:, 3],
        marker="o",
        linestyle="none",
        markersize=5,
        label="Observed Error",
    )
    ax.legend(framealpha=1.0, frameon=True)
    # Label the axes.
    ax.set_xlabel("$N$")
    ax.set_ylabel("Relative Forward Error")
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
    main(filename="curbed_condition.pdf")
