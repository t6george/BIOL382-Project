# Usage: python3 plot.py <path_to_csv> <output_image_path>
from sys import argv
import pandas as pd


if __name__ == "__main__":
    assert(len(argv) == 3)

    print("Reading csv...")

    df = pd.read_csv(argv[1])

    print("Converting time units...")
    # Convert to hours, we might want to make this days later
    df["t"] = df["t"].div(60 * 60)#.round(2)

    print("Generating figure...")
    fig = df.plot(x = "t", y = ["T4", "T3P", "T3c", "TSH", "TSHz", "T4th", "FT3", "FT4", "T3N", "T3R"],
        kind = "line", title = "Hormone Conservations vs. Time")

    fig.set_xlabel("Time [h]")
    fig.set_ylabel("Concentration")
    
    print("Saving figure...")
    fig.get_figure().savefig(argv[2])

    print("Done!")
