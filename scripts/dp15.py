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

import de_casteljau
import vs_method
import utils


F = fractions.Fraction
# f(s) = (s - 1/20)(s - 2/20) ... (s - 19/20)(s - 20/20)
WILKINSON1 = (
    float.fromhex("0x1.8e9b4e661311ep-26"),
    float.fromhex("-0x1.02de7c6051641p-24"),
    float.fromhex("0x1.1e770f8cd0529p-23"),
    float.fromhex("-0x1.1423646be98bbp-22"),
    float.fromhex("0x1.d65c5a1952526p-22"),
    float.fromhex("-0x1.654ea80cca589p-21"),
    float.fromhex("0x1.e74f07afab5fdp-21"),
    float.fromhex("-0x1.2b909657bbd25p-20"),
    float.fromhex("0x1.4cd9eaa240078p-20"),
    float.fromhex("-0x1.4e9588fcf799ep-20"),
    float.fromhex("0x1.302ad9a026e8fp-20"),
    float.fromhex("-0x1.f346dff3600b4p-21"),
    float.fromhex("0x1.70b1f41d35ef2p-21"),
    float.fromhex("-0x1.e74f07afab5fdp-22"),
    float.fromhex("0x1.1dd88670a1e07p-22"),
    float.fromhex("-0x1.25f9b84fd3738p-23"),
    float.fromhex("0x1.03e5133863565p-24"),
    float.fromhex("-0x1.7df414bbc06e1p-26"),
    float.fromhex("0x1.b3fd7328f4de7p-28"),
    float.fromhex("-0x1.3ee2a51e75a7fp-30"),
    float.fromhex("0x0.0p+0"),
)
# g(s) = (s - 2/2)(s - 2/4) ... (s - 2/2^{19})(s - 2/2^{20})
WILKINSON2 = (
    float.fromhex("0x1.0000000000000p-190"),
    float.fromhex("-0x1.9997800000000p-175"),
    float.fromhex("0x1.cbe01cc2d7943p-160"),
    float.fromhex("-0x1.5e5b5d5f563cfp-145"),
    float.fromhex("0x1.5faf54aece13bp-131"),
    float.fromhex("-0x1.c5acfc21c2483p-118"),
    float.fromhex("0x1.70884aaeef9a9p-105"),
    float.fromhex("-0x1.731dccc8a7da7p-93"),
    float.fromhex("0x1.c9d249cd6378bp-82"),
    float.fromhex("-0x1.57076b1c416fcp-71"),
    float.fromhex("0x1.3678531b90eb4p-61"),
    float.fromhex("-0x1.525677d4ef30bp-52"),
    float.fromhex("0x1.bb4180365b8a2p-44"),
    float.fromhex("-0x1.5cd503454aebdp-36"),
    float.fromhex("0x1.49772c71e764ap-29"),
    float.fromhex("-0x1.741a90baf536fp-23"),
    float.fromhex("0x1.f16d53866846bp-18"),
    float.fromhex("-0x1.7f99def2a0b19p-13"),
    float.fromhex("0x1.401687e02e12fp-9"),
    float.fromhex("-0x1.d926bcbd9b881p-7"),
    float.fromhex("0x0.0p+0"),
)
# h(s) = (s - 1/2)^{20} = [-1/2(1 - s) + 1/2s]^{20}
MULTIPLE_ROOT = (
    0.5 ** 20,
    -0.5 ** 20,
    0.5 ** 20,
    -0.5 ** 20,
    0.5 ** 20,
    -0.5 ** 20,
    0.5 ** 20,
    -0.5 ** 20,
    0.5 ** 20,
    -0.5 ** 20,
    0.5 ** 20,
    -0.5 ** 20,
    0.5 ** 20,
    -0.5 ** 20,
    0.5 ** 20,
    -0.5 ** 20,
    0.5 ** 20,
    -0.5 ** 20,
    0.5 ** 20,
    -0.5 ** 20,
    0.5 ** 20,
)


def do_plot(ax, s_vals, coeffs, title, add_legend=False, add_ylabel=False):
    bounds = []
    rel_errors1 = []
    rel_errors2 = []
    n = len(coeffs) - 1
    gamma3n = utils.gamma(3 * n)
    exact_coeffs = tuple(map(F, coeffs))
    abs_coeffs = tuple(map(abs, exact_coeffs))

    for s in s_vals:
        exact_s = F(s)
        exact_p = de_casteljau.basic(exact_s, exact_coeffs)
        if not isinstance(exact_p, F):
            raise TypeError(exact_p)

        exact_p_tilde = de_casteljau.basic(exact_s, abs_coeffs)
        a_priori_bound = gamma3n * exact_p_tilde / abs(exact_p)
        bounds.append(float(a_priori_bound))

        error1 = F(de_casteljau.basic(s, coeffs)) - exact_p
        rel_errors1.append(float(abs(error1 / exact_p)))

        error2 = F(vs_method.basic(s, coeffs)) - exact_p
        rel_errors2.append(float(abs(error2 / exact_p)))

    size = 5
    ax.semilogy(
        s_vals,
        bounds,
        marker="o",
        markersize=size,
        linestyle=":",
        alpha=0.5,
        color="black",
        label="Bound",
    )
    ax.semilogy(
        s_vals,
        rel_errors1,
        marker="o",
        markersize=size,
        label=r"$\mathtt{DeCasteljau}$",
    )
    ax.semilogy(
        s_vals,
        rel_errors2,
        marker="s",
        markersize=size,
        label=r"$\mathtt{VS}$",
    )
    if add_legend:
        ax.legend(
            loc="lower center",
            framealpha=1.0,
            frameon=True,
            markerscale=1.25,
            fontsize=16,
        )
    # Label the axes.
    ax.set_xlabel("$s$", fontsize=20)
    if add_ylabel:
        ax.set_ylabel("Relative Forward Error", fontsize=20)
    # Set the axis title.
    ax.set_title(title, fontsize=20)


def main(filename=None):
    figure, (ax1, ax2, ax3) = plt.subplots(1, 3)

    s_vals1 = [(2 * i + 1) / 72.0 for i in range(35 + 1)]
    do_plot(
        ax1, s_vals1, WILKINSON1, "$f(s)$", add_legend=True, add_ylabel=True
    )
    s_vals2 = [i / 39.0 for i in range(1, 38 + 1)]
    do_plot(ax2, s_vals2, WILKINSON2, "$g(s)$")
    s_vals3 = [4 * i / 100.0 for i in range(1, 24 + 1)]
    do_plot(ax3, s_vals3, MULTIPLE_ROOT, "$h(s)$")

    if filename is None:
        plt.show()
    else:
        figure.set_size_inches(12.99, 5.33)
        figure.subplots_adjust(
            left=0.06,
            bottom=0.09,
            right=0.98,
            top=0.96,
            wspace=0.15,
            hspace=0.20,
        )
        path = utils.get_path(filename)
        figure.savefig(path, bbox_inches="tight")
        print("Saved {}".format(filename))
        plt.close(figure)


if __name__ == "__main__":
    utils.set_styles()
    main(filename="compare_dp15.pdf")
