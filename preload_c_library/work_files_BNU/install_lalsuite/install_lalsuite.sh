#!/bin/bash

# 默认下载路径为当前用户目录下的/Downloads文件夹，可自行修改
# LALSUITE_PATH="${HOME}/Downloads"
LALSUITE_PATH=/workspace/chendawww

echo "正在进入lalsuite下载目录..."
cd ${LALSUITE_PATH}
echo "正在从git拉取lalsuite..."
git clone https://github.com/lscsoft/lalsuite.git
# git@github.com:lscsoft/lalsuite.git
# git clone https://gitlab.com/odidev/lalsuite.git

echo "正在进入lalsuite目录..."
cd lalsuite
echo "正在创建lalsuite的conda环境..."
conda env create -f ./conda-dev-env.yml
conda activate lalsuite-dev
echo "lalsuite-dev环境创建完毕 ( conda activate lalsuite-dev )"

PYTHON_PATH="${CONDA_PREFIX}/bin/python"
SWIG_PATH="${CONDA_PREFIX}/bin/swig"
PKG-CONFIG_PATH="${CONDA_PREFIX}/bin/pkg-config"

LALSUITE_INST_PATH="${CONDA_PREFIX}"

echo "正在安装lalsuite..."
echo "正在执行：./00boot"
./00boot
echo "正在执行：./configure [参数略]"
./configure --enable-lalsimulation --disable-swig --with-hdf5 --with-cfitsio \
    PYTHON=${PYTHON_PATH} \
    SWIG=${SWIG_PATH} \
    PKG_CONFIG=${PKG-CONFIG_PATH} \
    prefix=${LALSUITE_INST_PATH}
echo "正在执行：make"
make
echo "正在执行：make install"
make install
echo "安装完毕。"



./configure --enable-lalsimulation --disable-swig --with-hdf5 --with-cfitsio \
    PYTHON=/workspace/minconda3/envs/lalsuite-dev/bin/python \
    SWIG=/workspace/minconda3/envs/lalsuite-dev/bin/swig \
    PKG_CONFIG=/workspace/minconda3/envs/lalsuite-dev/bin/pkg-config \
    prefix=${CONDA_PREFIX}