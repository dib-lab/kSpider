
<p align="center">
  <img src="https://i.ibb.co/r66VhYc/6373059048-001abe61-1a3c-48c7-af51-0fd327b9c18a.png" alt="Logo"/>
</p>

<h1 align="center"> @dib-lab/kSpider </h1>
<p align="center">
<a href="https://github.com/dib-lab/kSpider/actions/workflows/ubuntu.yml"><img alt="Ubuntu" src="https://github.com/dib-lab/kSpider/actions/workflows/ubuntu.yml/badge.svg"></a>
<a href=""><img alt="Open Issues" src="https://img.shields.io/github/issues-raw/dib-lab/kSpider" height="20"/></a> <a href="https://github.com/dib-lab/kSpider/blob/master/LICENSE"><img alt="GitHub" src="https://img.shields.io/github/license/dib-lab/kSpider"></a> <a href="https://pypi.org/project/kSpider/#files"><img alt="PyPI - Wheel" src="https://img.shields.io/pypi/wheel/kSpider"></a> <a href=""><img alt="Maintained" src="https://img.shields.io/badge/Maintained%3F-yes-green.svg" height="20"/></a> <a href="https://pypi.org/project/kSpider"><img alt="PyPI - Python Version" src="https://img.shields.io/pypi/pyversions/kSpider"></a> 
</p>

<details>
<summary>📖 Table of Contents</summary>
<br />

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#table-of-contents)

## ➤ Table of Contents

- [➤ Table of Contents](#-table-of-contents)
- [➤ Introduction](#-introduction)
- [➤ Quick Installation (pip)](#-quick_installation)
- [➤ Build from source](#-build_source)
- [➤ Authors](#-authors)
- [➤ License](#-license)

</details>

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#introduction)

## ➤ Introduction

**kSpider** is a user-friendly command line interface program to perform sequence clustering. First, it creates an index using kProcessor for the source sequences. Second, it constructs a pairwise containment matrix through a single iteration over the index. Finally, it builds a graph from the pairwise matrix and applies a connected-components graph algorithm to extract the clusters with a user-defined containment threshold.

Documentations are hosted at [https://dib-lab.github.io/kSpider](https://dib-lab.github.io/kSpider)


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#quick_installation)

## ➤ Quick Installation (pip)

```bash
pip install kSpider
```

[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#build_source)

## ➤ Manual build / Development

### Install dependencies

```bash
sudo apt-get install g++ swig cmake python3-dev zlib1g-dev libghc-bzlib-dev python3-distutils libboost-all-dev
```

```bash
git clone https://github.com/dib-lab/kSpider.git
cd kSpider
git submodule update --init --recursive
cmake -Bbuild
cmake --build build
bash build_wrapper.sh
```


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#authors)

## ➤ Authors

| [<img alt="You?" src="https://avatars2.githubusercontent.com/u/7165864?s=460&&v=4" width="100">](https://github.com/mr-eyes) | [<img alt="Tamer Mansour" src="https://avatars3.githubusercontent.com/u/6537740?s=400&&v=4" width="100">](https://github.com/drtamermansour) |
|------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------|
| [Mohamed Abuelanin](https://github.com/mr-eyes)                                                                              | [Tamer Manosur](https://github.com/drtamermansour)                                                                                           |


[![-----------------------------------------------------](https://raw.githubusercontent.com/andreasbm/readme/master/assets/lines/colored.png)](#license)

## ➤ License

Licensed under [MIT License](https://opensource.org/licenses/MIT).
