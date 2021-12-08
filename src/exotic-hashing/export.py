import json
import math
import pandas as pd
import sys
import git
from pathlib import Path
from inspect import cleandoc

from plotly.subplots import make_subplots
import plotly.express as px
import plotly.graph_objects as go
import plotly
from pretty_html_table import build_table

# plot colors
pal = px.colors.qualitative.Plotly
color_sequence = ["#BBB", "#777", "#111", pal[9], pal[4], pal[6], pal[1], pal[0], "#58a2c4", pal[5], pal[2], pal[7], pal[8], pal[3]]

# plot labels
plot_labels = dict(
    cpu_time_per_key='ns per key',
    dataset_elem_count='dataset size',
    elem_magnitude='dataset size',
    hashfn_bits_per_key='bits per key',
    throughput='keys per second')


file = "results.json" if len(sys.argv) < 2 else sys.argv[1]
with open(file) as data_file:
    data = json.load(data_file)

    # convert json results to dataframe
    df = pd.json_normalize(data, 'benchmarks')

    # augment additional computed columns
    df["hashfn"] = df["label"].apply(lambda x : x.split(":")[0])
    df["dataset"] = df["label"].apply(lambda x : x.split(":")[1])
    df["probe_distribution"] = df["label"].apply(lambda x : x.split(":")[2] if len(x.split(":")) > 2 else "-")

    # order data (important for legend & colors)
    def order(x):
        x = x.lower()
        if x == "donothinghash":
            return 1
        if x == "rankhash":
            return 2
        if x == "recsplit_leaf12_bucket9":
            return 3
        if x == "compacttrie":
            return 4
        if x == "fastsuccincttrie":
            return 5
        if x == "simplehollowtrie":
            return 6
        if x == "hollowtrie":
            return 7
        if x == "mwhc":
            return 8
        if x == "compressedmwhc":
            return 9
        if x == "compactedmwhc":
            return 10
        if x == "rmirank":
            return 11
        if x == "compressedrmirank":
            return 12
        if x == "learnedlinear":
            return 13
        if x == "adaptivelearnedmmphf":
            return 14
        if x == "mapomphf":
            return 15
        return 0
    df["order"] = df.apply(lambda x : order(x["hashfn"]), axis=1)
    df = df.sort_values(by=["order", "dataset_elem_count"])

    # augment plotting datasets
    def magnitude(x):
        l = math.log(x, 10)
        rem = round(x/pow(10, l), 2)
        exp = int(round(l, 0))
        #return f'${rem} \cdot 10^{{{exp}}}$'
        return f'{rem}e-{exp}'
    df["elem_magnitude"] = df.apply(lambda x : magnitude(x["dataset_elem_count"]), axis=1)

    # prepare datasets for plotting & augment dataset specific columns
    lt_df = df[df["name"].str.lower().str.contains("lookuptime")].copy(deep=True)
    bt_df = df[df["name"].str.lower().str.contains("buildtime")].copy(deep=True)

    lt_df["cpu_time_per_key"] = lt_df['cpu_time']
    lt_df["throughput"] = lt_df.apply(lambda x : 10**9 / x["cpu_time_per_key"], axis=1)

    bt_df["cpu_time_per_key"] = bt_df.apply(lambda x : x["cpu_time"] / x["dataset_elem_count"], axis=1)
    bt_df["throughput"] = bt_df.apply(lambda x : 10**9 / x["cpu_time_per_key"], axis=1)
    bt_df["sorted"] = bt_df.apply(lambda x : x["name"].lower().startswith("presorted"), axis=1)

    # ensure export output folder exists
    results_path = "docs" if len(sys.argv) < 3 else sys.argv[2]
    Path(results_path).mkdir(parents=True, exist_ok=True)

    def convert_to_html(fig):
        #fig.show()
        return fig.to_html(full_html=False, include_plotlyjs=False)

    def plot_lookup_times():
        name = "lookup_time"
        fig = px.line(
            lt_df,
            x="dataset_elem_count",
            y="cpu_time_per_key",
            color="hashfn",
            facet_row="probe_distribution",
            facet_col="dataset",
            category_orders={"dataset": ["seq", "gap_10", "uniform", "normal", "wiki", "osm", "fb"]},
            markers=True,
            log_x=True,
            labels=plot_labels,
            color_discrete_sequence=color_sequence,
            height=600,
            title="Lookup - nanoseconds per key"
            )
        return convert_to_html(fig)

    def plot_hashfn_bits_per_key():
        name = "bits_per_key"
        fig = px.line(
            lt_df,
            x="dataset_elem_count",
            y="hashfn_bits_per_key",
            color="hashfn",
            facet_col="dataset",
            facet_col_wrap=3,
            category_orders={"dataset": ["seq", "gap_10", "uniform", "normal", "wiki", "osm", "fb"]},
            log_x=True,
            markers=True,
            labels=plot_labels,
            color_discrete_sequence=color_sequence,
            height=600,
            title="Space - total bits per key"
            )
        fig.update_yaxes(range=[-50, 700])
        return convert_to_html(fig)

    def plot_build_time():
        # copy to enable value changes
        f_bt_df = bt_df.copy(deep=True)
        #f_bt_df = f_bt_df[f_bt_df["dataset_elem_count"].isin([10**6, 10**8])]
        f_bt_df = f_bt_df[f_bt_df["dataset_elem_count"] > 9 * 10**7]
        name = "build_time"
        fig = px.bar(
            f_bt_df,
            x="elem_magnitude",
            y="throughput",
            color="hashfn",
            barmode="group",
            facet_col="dataset",
            facet_row="sorted",
            category_orders={"dataset": ["seq", "gap_10", "uniform", "normal", "wiki", "osm", "fb"]},
            labels=plot_labels,
            color_discrete_sequence=color_sequence,
            height=600,
            title="Build - throughput in keys per second"
            )
        fig.update_traces(
                patch={'visible': 'legendonly'},
                selector=lambda go : go.legendgroup.lower() in ["donothinghash"])
        return convert_to_html(fig)

    def plot_raw_data():
        raw_data = df.sort_values(by=["name"])
        raw_data = raw_data.rename({"cpu_time": "ns", 'hashfn': 'function', "probe_distribution": "probe distribution", 'dataset_elem_count': 'keys', 'hashfn_bits_per_key': 'bits per key'}, axis='columns')
        raw_data = raw_data[["name", "function", "probe distribution", "dataset", "keys", "bits per key", "ns"]]
        raw_data["ns"] = raw_data.apply(lambda x : str(int(float(x["ns"]))), axis=1)
        raw_data["keys"] = raw_data.apply(lambda x : str(int(x["keys"])), axis=1)

        return cleandoc(f"""
        <div style="width: 100%; height: 500px; overflow-y: scroll;">
            {build_table(raw_data, 'blue_light', width="100%")}
        </div>
        """)

    with open(f'{results_path}/index.html', 'w') as readme:
        readme.write(cleandoc(f"""
        <!doctype html>
        <html>
          <head>
              <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
          </head>

          <body style="display: grid; grid-template-columns: repeat(auto-fit, minmax(1200px, 1fr))">
            <embed src="functions.html" style="width: 100%; height: 500px;"/>
            {plot_lookup_times()}
            {plot_hashfn_bits_per_key()}
            {plot_build_time()}
            <div style="margin: 15px">
                <h3 style="color: rgb(42, 63, 95)">Raw Data</h2>
                {plot_raw_data()}
            </div>
          </body>
        </html>
        """))
