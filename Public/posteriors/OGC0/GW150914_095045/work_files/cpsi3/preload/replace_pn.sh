#!/bin/bash

WORK_PATH=$(pwd)

# CONFIG_FILE="$WORK_PATH/preload/replace_pn.conf"

GCC_PATH="/usr/bin/gcc"
PYTHON_PATH="/usr/bin/python"

LALSUITE_c_PATH="/usr/include/lalsuite"
LALSUITE_c_INST_PATH="${LALSUITE_c_PATH}/_inst"

PRELOAD_PATH="${WORK_PATH}/preload"
PRELOAD_FILE_NAME_c="my_dchii.c"
PRELOAD_FILE_NAME_so="my_dchii.so"

TEST_PATH="${WORK_PATH}/test"
TEST_FILE_NAME_c="test_preload.c"
TEST_FILE_NAME_py="test_preload.py"


# 初始化所有配置变量
# for key in "${!DEFAULTS[@]}"; do
#     declare -g "$key"="${DEFAULTS[$key]}"
# done

export LD_LIBRARY_PATH="$(echo -n "$LD_LIBRARY_PATH:$LALSUITE_c_INST_PATH/lib" | tr ':' '\n' | awk '!x[$0]++' | tr '\n' ':' | sed 's/:$//')"

######################################################

# 加载配置文件
load_config() {
    local CONFIG_FILE="$1"

    # CONFIG_FILE=${CONFIG_FILE:-"$2"}

    # 检查配置文件是否存在
    # if [ ! -f "$CONFIG_FILE" ]; then
    #     echo "警告: 配置文件不存在 '$CONFIG_FILE'" >&2
    #     return 1
    # fi
    
    while IFS='=' read -r key value; do
        # 跳过注释和空行
        if [[ "$key" =~ ^[^#] && -n "$key" ]]; then
            # 安全赋值（避免代码注入）
            case $key in
                WORK_PATH| \
                GCC_PATH|PYTHON_PATH| \
                LALSUITE_c_PATH|LALSUITE_c_INST_PATH| \
                PRELOAD_PATH|PRELOAD_FILE_NAME_c|PRELOAD_FILE_NAME_so| \
                TEST_PATH|TEST_FILE_NAME_c|TEST_FILE_NAME_py)
                    declare -g "$key"="$value"
                    ;;
                *)
                    echo "警告: 忽略未知配置项 '$key'"
                    ;;
            esac
        fi
    done < "$CONFIG_FILE"
}

# 解析命令行参数
while [[ $# -gt 0 ]]; do
    case "$1" in
        --work-path)
            WORK_PATH="$2"
            shift 2
            ;;
        --gcc-path)
            GCC_PATH="$2"
            shift 2
            ;;
        --python-path)
            PYTHON_PATH="$2"
            shift 2
            ;;
        --lalsuite-path)
            LALSUITE_c_PATH="$2"
            LALSUITE_c_INST_PATH="${LALSUITE_c_PATH}/_inst"
            shift 2
            ;;
        --preload-path)
            PRELOAD_PATH="$2"
            shift 2
            ;;
        --test-path)
            TEST_PATH="$2"
            shift 2
            ;;
        --config)
            load_config "$2"
            shift 2
            ;;
        *)
            echo "错误: 未知选项 $1" >&2
            exit 1
            ;;
    esac
done

# 验证路径是否存在
check_path() {
    if [ ! -e "$1" ]; then
        echo "错误: 路径 '$1' 不存在" >&2
        exit 1
    fi
}

check_path "$GCC_PATH"
check_path "$PYTHON_PATH"
check_path "$LALSUITE_c_PATH"
check_path "$PRELOAD_PATH"
check_path "$TEST_PATH"

######################################################

# Go into work directory where there should be preloading .c file and test .py file.
cd "$WORK_PATH" || { echo "错误: 无法进入目录 $WORK_PATH" >&2; exit 1; }

echo "正在编译 preload 预加载库..."


${GCC_PATH} -shared -fPIC -o "${PRELOAD_PATH}/${PRELOAD_FILE_NAME_so}" \
    "${PRELOAD_PATH}/${PRELOAD_FILE_NAME_c}" \
    -I${LALSUITE_c_PATH}/lalsimulation/include \
    -I${LALSUITE_c_PATH}/lalsimulation/lib \
    -I${LALSUITE_c_PATH}/lalsimulation/include/lal \
    -I${LALSUITE_c_PATH}/lal/include \
    -I${LALSUITE_c_PATH}/lal/include/lal \
    -I${LALSUITE_c_PATH}/lal/lib \
    -I${LALSUITE_c_INST_PATH}/include \
    -L/usr/lib -L${LALSUITE_c_INST_PATH}/lib \
    -Wl,-rpath,/root/Downloads/lalsuite/_inst/lib \
    -lm -ldl \
    -llalsimulation -llalsupport -llal \
    -D_GNU_SOURCE

# ${GCC_PATH} -shared -fPIC -o "${PRELOAD_PATH}/${PRELOAD_FILE_NAME_so}" \
#     "${PRELOAD_PATH}/${PRELOAD_FILE_NAME_c}" \
#     -I/usr/include \
#     -I${LALSUITE_c_PATH}/lalsimulation/include \
#     -I${LALSUITE_c_PATH}/lalsimulation/lib \
#     -I${LALSUITE_c_PATH}/lalsimulation/include/lal \
#     -I${LALSUITE_c_PATH}/lal/include \
#     -I${LALSUITE_c_PATH}/lal/include/lal \
#     -I${LALSUITE_c_PATH}/lal/lib \
#     -I${LALSUITE_c_INST_PATH}/include \
#     -L/usr/lib -L${LALSUITE_c_INST_PATH}/lib \
#     -Wl,-rpath,/root/Downloads/lalsuite/_inst/lib \
#     -lm -ldl \
#     -llalsimulation -llalsupport -llal


# ${GCC_PATH} -shared -fPIC -o "${PRELOAD_PATH}/LALSimIMRPhenomX_inspiral_my.so" "${PRELOAD_PATH}/preloaded_libraries/LALSimIMRPhenomX_inspiral_my.c" \
#     -I/usr/include \
#     -I${LALSUITE_c_PATH}/lalsimulation/include \
#     -I${LALSUITE_c_PATH}/lalsimulation/lib \
#     -I${LALSUITE_c_PATH}/lalsimulation/include/lal \
#     -I${LALSUITE_c_PATH}/lal/include \
#     -I${LALSUITE_c_PATH}/lal/include/lal \
#     -I${LALSUITE_c_PATH}/lal/lib \
#     -I/workspace/minconda3/envs/preload/lib \
#     -I/workspace/minconda3/envs/preload/include \
#     -I${LALSUITE_c_INST_PATH}/include \
#     -I${LALSUITE_c_INST_PATH}/include/lal \
#     -L/usr/lib -L${LALSUITE_c_INST_PATH}/lib \
#     -lm -ldl \
#     -llalsimulation -llalsupport -llal

#     # -I${LALSUITE_c_INST_PATH}/include \
#     # -I${LALSUITE_c_INST_PATH}/include/lal \

# ${GCC_PATH} -shared -fPIC -o "${PRELOAD_PATH}/LALSimIMRPhenomX_internals_my.so" "${PRELOAD_PATH}/preloaded_libraries/LALSimIMRPhenomX_internals_my.c" \
#     -I/usr/include \
#     -I${LALSUITE_c_PATH}/lalsimulation/include \
#     -I${LALSUITE_c_PATH}/lalsimulation/lib \
#     -I${LALSUITE_c_PATH}/lalsimulation/include/lal \
#     -I${LALSUITE_c_PATH}/lal/include \
#     -I${LALSUITE_c_PATH}/lal/include/lal \
#     -I${LALSUITE_c_PATH}/lal/lib \
#     -I/workspace/minconda3/envs/preload/lib \
#     -I/workspace/minconda3/envs/preload/include \
#     -I${LALSUITE_c_INST_PATH}/include \
#     -I${LALSUITE_c_INST_PATH}/include/lal \
#     -L/usr/lib -L${LALSUITE_c_INST_PATH}/lib \
#     -lm -ldl \
#     -llalsimulation -llalsupport -llal

    
# ${GCC_PATH} -shared -fPIC -o "${PRELOAD_PATH}/LALSimIMRPhenomX_my.so" "${PRELOAD_PATH}/preloaded_libraries/LALSimIMRPhenomX_my.c" \
#     -I/usr/include \
#     -I${LALSUITE_c_PATH}/lalsimulation/include \
#     -I${LALSUITE_c_PATH}/lalsimulation/lib \
#     -I${LALSUITE_c_PATH}/lalsimulation/include/lal \
#     -I${LALSUITE_c_PATH}/lal/include \
#     -I${LALSUITE_c_PATH}/lal/include/lal \
#     -I${LALSUITE_c_PATH}/lal/lib \
#     -I/workspace/minconda3/envs/preload/lib \
#     -I/workspace/minconda3/envs/preload/include \
#     -I${LALSUITE_c_INST_PATH}/include \
#     -I${LALSUITE_c_INST_PATH}/include/lal \
#     -L/usr/lib  -L${LALSUITE_c_INST_PATH}/lib \
#     -lm -ldl \
#     -llalsimulation -llalsupport -llal

if [ $? -ne 0 ]; then
    echo "错误: 编译失败" >&2
    exit 1
fi


# 设置 LD_PRELOAD（自动去重）
# export LD_PRELOAD="$(echo -n "$LD_PRELOAD:$WORK_PATH/$PRELOAD_FILE_NAME_so" | tr ':' '\n' | awk '!x[$0]++' | tr '\n' ':' | sed 's/:$//')"


echo "-------------------------"
echo "1. 配置完成！请先查看环境变量："
echo "echo \$LD_LIBRARY_PATH"
echo "echo \$LD_PRELOAD"
echo "并完成环境变量配置："
echo "export LD_LIBRARY_PATH=${LALSUITE_c_INST_PATH}/lib:\$LD_LIBRARY_PATH"
echo "export LD_PRELOAD=${PRELOAD_PATH}/${PRELOAD_FILE_NAME_so}:\$LD_PRELOAD"
echo "-------------------------"
echo "2. 通过以下命令运行gcc测试："
echo "$GCC_PATH $TEST_PATH/$TEST_FILE_NAME_c"
echo "以检查lalsuite的c语言库是否安装成功、环境变量是否配置成功、${PRELOAD_PATH}/${PRELOAD_FILE_NAME_so}文件是否成功preload；"
echo "-------------------------"
echo "3. 通过以下命令运行Python测试："
echo "$PYTHON_PATH $TEST_PATH/$TEST_FILE_NAME_py"
echo "以检查lalsuite的c语言库与PyCBC的lalsuite动态链接so库是否头文件版本一致、环境变量是否配置成功、${PRELOAD_PATH}/${PRELOAD_FILE_NAME_so}文件preload是否成功应用于PyCBC的lalsuite动态链接so库；"
echo "-------------------------"
echo "若出现问题，请使用 ldd xxx.so 与 nm -D xxx.so 等命令检查动态链接库链接情况、以及版本的适配。"

source ~/.bashrc

export C_INCLUDE_PATH=/usr/include/lalsuite/lal/include:$C_INCLUDE_PATH
export C_INCLUDE_PATH=/usr/include/lalsuite/lalsimulation/include/lal:$C_INCLUDE_PATH
export C_INCLUDE_PATH=/usr/include/lalsuite/lalsimulation/include:$C_INCLUDE_PATH
export C_INCLUDE_PATH=/usr/include/lalsuite/lalsimulation/lib:$C_INCLUDE_PATH