#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
用法：
  download_core_nt.sh [输出目录]

默认输出目录：
  项目下的 02_resources/blast/core_nt

说明：
  下载并解压 NCBI 预格式化 core_nt BLAST 数据库，随后检查数据库完整性。
  再次运行会由 update_blastdb.pl 检查并下载有更新的数据库文件。
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if (( $# > 1 )); then
    usage >&2
    exit 2
fi

for command_name in update_blastdb.pl blastdbcheck; do
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        printf '错误：未找到命令 %s。\n' "${command_name}" >&2
        printf '请先安装 NCBI BLAST+：sudo apt update && sudo apt install ncbi-blast+\n' >&2
        exit 1
    fi
done

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -- "${script_dir}/../../.." && pwd)
output_dir=${1:-"${project_dir}/02_resources/blast/core_nt"}

mkdir -p -- "${output_dir}"
output_dir=$(cd -- "${output_dir}" && pwd)

printf '下载目录：%s\n' "${output_dir}"
printf '可用空间：'
df -h --output=avail "${output_dir}" | awk 'NR == 2 {print $1}'

(
    cd -- "${output_dir}"
    update_blastdb.pl --decompress core_nt
    blastdbcheck -db core_nt -dbtype nucl
)

printf '完成：core_nt 已下载并通过完整性检查：%s\n' "${output_dir}"
