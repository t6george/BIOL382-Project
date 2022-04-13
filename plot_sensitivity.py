# Usage: python3 plot_sensitivity.py <path_to_csv> <output_fpath>

from sys import argv
import pandas as pd

if __name__ == "__main__":
    assert(len(argv) == 3)

    print("Reading csv...")

    df = pd.read_csv(argv[1])

    print("Computing numerical derivatives...")

    for i in range(len(df["TSH"]) - 1):
        df["TSH"][i] = (df["TSH"][i + 1] - df["TSH"][i]) / (df["GT"][i + 1] - df["GT"][i])

    df.drop(df.tail(1).index, inplace = True)

    print("Generating figure...")

    fig = df.plot(x = "GT", y = ["TSH"],
        kind = "line", title = "TSH Sensitivity to Changes in GT", legend = False)

    fig.set_xlabel("GT [mol/s]")
    fig.set_ylabel("Absolute Sensitivity of TSH")

    print("Saving figure...")

    fig.get_figure().savefig(argv[2])

    print("Done!")

