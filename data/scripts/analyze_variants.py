from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt

project_dir = Path(__file__).resolve().parents[2]
data_dir = project_dir / "data"
vep_files = sorted(data_dir.glob("*_VEP.tsv"))
if len(vep_files) != 1:
    raise FileNotFoundError(
        f"Se esperaba un único TSV de VEP en {data_dir}; encontrados: {len(vep_files)}"
    )

results_dir = project_dir / "results"
results_dir.mkdir(exist_ok=True)

df = pd.read_csv(vep_files[0], sep="\t", na_values="-", low_memory=False)
df.columns = df.columns.str.lstrip("#")
required_columns = {
    "SYMBOL",
    "CANONICAL",
    "Uploaded_variation",
    "Location",
    "Ref",
    "Allele",
    "HGVSc",
    "ClinVar_CLNSIG",
    "gnomAD_exomes_AF",
    "REVEL",
    "am_pathogenicity",
}
missing_columns = required_columns.difference(df.columns)
if missing_columns:
    raise ValueError(f"Faltan columnas requeridas en el TSV: {', '.join(sorted(missing_columns))}")

df = df[(df["SYMBOL"] == "MYH7") & (df["CANONICAL"] == "YES")].copy()
df = df.drop_duplicates(subset=["Location", "Ref", "Allele"])
if df.empty:
    raise ValueError("No se encontraron variantes MYH7 canónicas en el TSV.")

vcf_path = data_dir / "clinvar_labels_grch38.vcf"
if vcf_path.exists():
    hgvs_by_variant = {}
    with vcf_path.open(encoding="utf-8") as vcf_file:
        for line in vcf_file:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 8:
                raise ValueError(f"Registro VCF no válido en {vcf_path}: {line.rstrip()}")
            info = dict(
                item.split("=", 1)
                for item in fields[7].split(";")
                if "=" in item
            )
            if "HGVS" in info:
                hgvs_by_variant[(fields[0], fields[1], fields[3], fields[4])] = info["HGVS"]

    def find_input_hgvs(row):
        chromosome, _, location = str(row["Location"]).partition(":")
        position = location.split("-", 1)[0]
        key = (chromosome, position, str(row["Ref"]), str(row["Allele"]))
        return hgvs_by_variant.get(key)

    input_hgvs = df.apply(find_input_hgvs, axis=1).astype("string")
    df["HGVSc"] = df["HGVSc"].astype("string").fillna(input_hgvs)

for col in ["REVEL", "am_pathogenicity", "gnomAD_exomes_AF"]:
    df[col] = pd.to_numeric(df[col], errors="coerce")

df["Rare"] = df["gnomAD_exomes_AF"].isna() | (df["gnomAD_exomes_AF"] < 0.001)

cols = [
    "Uploaded_variation",
    "Location",
    "Ref",
    "Allele",
    "HGVSc",
    "ClinVar_CLNSIG",
    "gnomAD_exomes_AF",
    "Rare",
    "REVEL",
    "am_pathogenicity",
]
df[cols].to_csv(results_dir / "vep_clean.csv", index=False)

df.boxplot(column=["REVEL", "am_pathogenicity"], by="ClinVar_CLNSIG", rot=45, figsize=(10, 4))
plt.suptitle("REVEL y AlphaMissense vs ClinVar")
plt.tight_layout()
plt.savefig(results_dir / "benchmark.png", dpi=300)
plt.close()

print(f"¡Listo! {len(df)} variantes; resultados guardados en {results_dir}")
