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
import math
import os


U = fractions.Fraction(1, 2 ** 53)


def binomial(n, k):
    numerator = math.factorial(n)
    denominator = math.factorial(k) * math.factorial(n - k)
    result = fractions.Fraction(numerator, denominator)
    if float(result) != result:
        raise ValueError("Cannot be represented exactly")
    return float(result)


def set_styles():
    """Set the styles used for plotting."""
    import seaborn

    seaborn.set(style="white")


def get_path(filename):
    """Get a file path in the ``images/`` directory.

    This assumes the script is currently in the ``src/``
    directory.
    """
    curr_dir = os.path.abspath(os.path.dirname(__file__))
    root_dir = os.path.dirname(curr_dir)
    images_dir = os.path.join(root_dir, "images")
    return os.path.join(images_dir, filename)


def gamma(n):
    return (n * U) / (1 - n * U)
