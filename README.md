# A Curious Case of Curbed Condition

This is the repository for a [paper][1] I wrote and posted on the arXiv.

This repository is laid out in a manner described in
[Good Enough Practices in Scientific Computing][2].

## Abstract

In computer aided geometric design a polynomial is usually represented in
Bernstein form. The de Casteljau algorithm is the most well-known algorithm
for evaluating a polynomial in this form. Evaluation via the de Casteljau
algorithm has relative forward error proportional to the condition number
of evaluation. However, for a particular family of polynomials, a curious
phenomenon occurs: the observed error is much smaller than the expected
error bound. We examine this family and prove a much stronger error bound
than the one that applies to the general case. Then we provide a few examples
to demonstrate the difference in rounding.

## Installation

The code used to build the manuscript is written in Python. To run the
code, Python 3.6 should be installed, along with ``nox-automation``:

```
python -m pip install --upgrade nox-automation
```

Once installed, the various build jobs can be listed.

```
$ nox --list-sessions
Available sessions:
* build_tex
* make_images
* update_requirements
```

To run ``nox -s build_tex`` (i.e. to build the PDF), ``pdflatex`` is required.

[1]: https://arxiv.org/abs/1806.05145
[2]: https://arxiv.org/abs/1609.00037
