import yaml

with open("config.yaml") as f:
    config = yaml.safe_load(f)

data_repo = "../data".strip()
cell_line_sra_pairs = [(cell_line, sra) for cell_line, sras in config["metadata"].items() for sra in sras]
cell_line_single = config["metadata"].keys()