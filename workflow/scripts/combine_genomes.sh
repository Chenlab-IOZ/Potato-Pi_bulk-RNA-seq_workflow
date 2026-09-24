#!/usr/bin/env bash

set -euo pipefail

usage() {
    cat <<'EOF'
用法：
  combine_genomes.sh [资源目录] [输出 FASTA]

默认值：
  资源目录：项目下的 02_resources
  输出 FASTA：<资源目录>/combined_genome.fa

序列名前缀：
  Solanum_tuberosum.SolTub_3.0.dna.toplevel.fa       Stu_
  Phytophthora_infestans.ASM14294v1.dna.toplevel.fa Pinf_
  potato_organelle.rename.fa                         Stu_

脚本仅保留 FASTA 标题行中的首个序列 ID，删除其后的描述。
已有正确前缀的 ID 不会被重复添加。
EOF
}

if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
fi

if (( $# > 2 )); then
    usage >&2
    exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
project_dir=$(cd -- "${script_dir}/../../.." && pwd)
resource_dir=${1:-"${project_dir}/02_resources"}
output_fasta=${2:-"${resource_dir}/combined_genome.fa"}

input_names=(
    "Solanum_tuberosum.SolTub_3.0.dna.toplevel.fa"
    "potato_organelle.rename.fa"
    "Phytophthora_infestans.ASM14294v1.dna.toplevel.fa"
)
prefixes=(
    "Stu_"
    "Stu_"
    "Pinf_"
)

if [[ ! -d "${resource_dir}" ]]; then
    printf '错误：资源目录不存在：%s\n' "${resource_dir}" >&2
    exit 1
fi

output_dir=$(dirname -- "${output_fasta}")
mkdir -p -- "${output_dir}"
tmp_fasta=$(mktemp --tmpdir="${output_dir}" '.combined_genome.XXXXXX.fa')
trap 'rm -f -- "${tmp_fasta}"' EXIT

expected_records=0

for i in "${!input_names[@]}"; do
    input_fasta="${resource_dir}/${input_names[$i]}"
    prefix=${prefixes[$i]}

    if [[ ! -r "${input_fasta}" ]]; then
        printf '错误：FASTA 不存在或不可读：%s\n' "${input_fasta}" >&2
        exit 1
    fi

    record_count=$(awk '/^>/{n++} END{print n+0}' "${input_fasta}")
    if (( record_count == 0 )); then
        printf '错误：FASTA 中没有序列记录：%s\n' "${input_fasta}" >&2
        exit 1
    fi
    expected_records=$((expected_records + record_count))

    awk -v prefix="${prefix}" '
        /^>/ {
            header = substr($0, 2)
            separator = match(header, /[[:space:]]/)
            if (separator == 0) {
                id = header
            } else {
                id = substr(header, 1, separator - 1)
            }
            if (index(id, prefix) != 1) {
                id = prefix id
            }
            print ">" id
            next
        }
        { print }
    ' "${input_fasta}" >> "${tmp_fasta}"

    printf '已处理：%s（%d 条序列，前缀 %s）\n' \
        "${input_names[$i]}" "${record_count}" "${prefix}"
done

actual_records=$(awk '/^>/{n++} END{print n+0}' "${tmp_fasta}")
if (( actual_records != expected_records )); then
    printf '错误：合并后序列数异常（预期 %d，实际 %d）\n' \
        "${expected_records}" "${actual_records}" >&2
    exit 1
fi

if ! awk '
    /^>/ {
        id = substr($1, 2)
        if (seen[id]++) {
            printf "错误：合并后序列 ID 重复：%s\n", id > "/dev/stderr"
            duplicate = 1
        }
    }
    END { exit duplicate }
' "${tmp_fasta}"; then
    exit 1
fi

chmod 0644 -- "${tmp_fasta}"
mv -f -- "${tmp_fasta}" "${output_fasta}"
trap - EXIT

printf '完成：%s（共 %d 条序列）\n' "${output_fasta}" "${actual_records}"
