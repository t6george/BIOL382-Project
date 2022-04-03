# Usage: python3 plot.py <path_to_csv> <output_directory for images>
from os.path import join
from sys import argv
import pandas as pd


if __name__ == "__main__":
    assert(len(argv) == 3)

    print("Reading csv...")

    df = pd.read_csv(argv[1])

    print("Converting time units...")
    # Convert to hours, we might want to make this days later
    # df["t"] = df["t"].div(60 * 60)#.round(2)

    print("Generating figure...")

    for i, col in enumerate(df):
        if i > 0:
            fig = df.plot(x = "t", y = [col],
                kind = "line", title = "Hormone Conservations vs. Time")

            fig.set_xlabel("Time [h]")
            fig.set_ylabel("Concentration [units vary]")
        
            print(f"Saving figure {i} / {len(df.columns) - 1}")
            fig.get_figure().savefig(join(argv[2], f"{col}.png"))

    print("Done!")
