#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
テスト: .cir → .asc 自動変換
各テストケースのネットリストから.ascファイルを生成する
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', 'src'))

from netlist_to_asc import NetlistToAsc

# テスト用ネットリスト
TEST_CASES = {
    'test01_rc_lowpass': """\
* RC Lowpass Filter
V1 in 0 AC 1
R1 in out 1k
C1 out 0 1u
.ac dec 100 1 100k
.end
""",

    'test02_rlc_series': """\
* RLC Series Circuit
I1 0 in AC 1
R1 in mid 100
L1 mid out 10m
C1 out 0 1u
.ac dec 100 10 100k
.end
""",

    'test03_parallel_lc': """\
* Parallel LC Tank
I1 0 top AC 1
L1 top 0 1m
C1 top 0 100n
.ac dec 100 1k 100k
.end
""",

    'test04_pi_filter': """\
* Pi Filter (C-L-C)
I1 0 in AC 1
C1 in 0 100n
L1 in out 1m
C2 out 0 100n
R1 out 0 50
.ac dec 100 1k 1Meg
.end
""",

    'test05_voltage_rc': """\
* Voltage Source RC
V1 in 0 AC 1
R1 in out 10k
C1 out 0 10n
.ac dec 100 100 1Meg
.end
""",
}


def main():
    converter = NetlistToAsc()
    output_dir = os.path.dirname(__file__)

    for name, netlist in TEST_CASES.items():
        print(f"\n{'='*60}")
        print(f"Converting: {name}")
        print(f"{'='*60}")

        # .cirファイルを書き出し
        cir_path = os.path.join(output_dir, f'{name}.cir')
        with open(cir_path, 'w') as f:
            f.write(netlist)

        # .ascに変換
        asc_path = os.path.join(output_dir, f'{name}.asc')
        asc = converter.convert_file(cir_path, asc_path)

        print(f"\nGenerated ASC:")
        print(asc)
        print()


if __name__ == '__main__':
    main()
