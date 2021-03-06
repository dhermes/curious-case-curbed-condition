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

"""This is a configuration file for running ``nox`` on this project.

To determine the supported actions run ``nox --list-sessions`` from the
project root.
"""

import os

import nox
import py.path


NOX_DIR = os.path.abspath(os.path.dirname(__file__))
SINGLE_INTERP = "python3.6"


def get_path(*names):
    return os.path.join(NOX_DIR, *names)


class Remove(object):
    def __init__(self, prefix, extensions):
        self.prefix = prefix
        self.extensions = extensions

    def __call__(self):
        for extension in self.extensions:
            path = "{}.{}".format(self.prefix, extension)
            os.remove(path)


def build_tex_file(session, base, new_id, extensions=()):
    # NOTE: This assumes that ``session.chdir(get_path('doc'))``
    #       has been called.
    modify_id = get_path("scripts", "modify_pdf_id.py")

    session.run("pdflatex", base)
    session.run("bibtex", base)
    session.run("pdflatex", base)
    session.run("bibtex", base)
    session.run("pdflatex", base)
    session.run("pdflatex", base)

    path = get_path("doc", base)
    remove = Remove(path, extensions)
    session.run(remove)
    session.run("python", modify_id, "--base", path, "--id", new_id)


@nox.session
def build_tex(session):
    session.interpreter = SINGLE_INTERP

    if py.path.local.sysfind("pdflatex") is None:
        session.skip("`pdflatex` must be installed")

    if py.path.local.sysfind("bibtex") is None:
        session.skip("`bibtex` must be installed")

    # No need to create a virtualenv.
    session.virtualenv = False

    session.chdir(get_path("doc"))

    build_tex_file(
        session,
        "curious-case",
        "663FB1F131563E4B4B4208A6D18A69EC",
        extensions=("aux", "bbl", "blg", "log", "out", "toc"),
    )


@nox.session
def make_images(session):
    session.interpreter = SINGLE_INTERP
    # Install all dependencies.
    session.install("--requirement", "make-images-requirements.txt")
    # Run the script(s).
    # Make sure
    # - Custom ``matplotlibrc`` is used
    # - Code in ``src/`` is importable
    # - PDFs have deterministic ``CreationDate``
    env = {
        "MATPLOTLIBRC": get_path("images"),
        "PYTHONPATH": get_path("src"),
        "SOURCE_DATE_EPOCH": "0",
    }
    names = (
        "dp15.py",
        "curbed_errors.py",
        "curious_intro.py",
    )
    for name in names:
        script = get_path("scripts", name)
        session.run("python", script, env=env)


@nox.session
def update_requirements(session):
    session.interpreter = SINGLE_INTERP

    if py.path.local.sysfind("git") is None:
        session.skip("`git` must be installed")

    # Install all dependencies.
    session.install("pip-tools")

    # Update all of the requirements file(s).
    names = ("make-images",)
    for name in names:
        in_name = "{}-requirements.in".format(name)
        txt_name = "{}-requirements.txt".format(name)
        session.run(
            "pip-compile", "--upgrade", "--output-file", txt_name, in_name
        )
        session.run("git", "add", txt_name)
