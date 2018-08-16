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

r"""Evaluates a polynomial in Bernstein form via the `VS method`_.

.. _VS method: https://doi.org/10.1016/0167-8396(86)90018-X
.. _survey: http://dx.doi.org/10.1016/j.amc.2015.08.086

This method has been described in the 1D case in a `survey`_
paper that followed the original formulation.
"""

import utils


def basic(s, coeffs):
    n = len(coeffs) - 1
    r = 1.0 - s
    if s >= 0.5:
        sigma = r / s
        multiplier = s
    else:
        sigma = s / r
        multiplier = r
        coeffs = coeffs[::-1]

    result = coeffs[0]
    for j in range(1, n + 1):
        binom_val = utils.binomial(n, j)
        modified_coeff = binom_val * coeffs[j]
        result = result * sigma + modified_coeff

    for _ in range(n):
        result = multiplier * result

    return result
