#!/usr/bin/env python3
import pandas as pd
import sys

main_file = sys.argv[1]
src_file = sys.argv[2]
output_file = sys.argv[3]

df_main = pd.read_csv(main_file)
df_src = pd.read_csv(src_file)

df_merged = pd.concat([df_main[df_main['kinase'] != 'SRC'], df_src], ignore_index=True)
df_merged.to_csv(output_file, index=False)

print(f"Done. Output: {output_file}")
